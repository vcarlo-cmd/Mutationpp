#!/usr/bin/env python
"""
Génère et trace la table B' du carbure de silicium (SiC) dans l'air.
Utilise le binaire C++ MutationPP `bprime` (bilan de masse généralisé,
options -char/-char-elem) comme moteur de calcul.

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique) et une plage de températures de 300 à 5000 K.

Deux particularités par rapport au carbone et à la silice :

  1. Élément de suivi du bilan = N (et non Si).
     Le SiC ne se consomme pas toujours de façon congruente : selon (Tw, P),
     il laisse à la paroi de la silice (oxydation passive), du silicium
     liquide ou un résidu de carbone (décomposition SiC -> Si(g) + C(gr)).
     Suivre Si ou C donne alors un B'c faux (nul, voire négatif puis écrêté,
     alors que le SiC se décompose). L'azote de l'air est inerte et ne
     figure dans aucune phase condensée : son bilan donne directement la
     masse nette gazéifiée par la paroi,
         B'c = y_e,N / y_w,N - 1
     Pour le carbone et la silice, où le char est la seule phase condensée,
     cette définition est identique à celle des autres tables.

  2. Calcul point par point, avec contrôle de stabilité des phases.
     Le solveur multiphase de Mutation++ boucle indéfiniment sur ~7 % des
     points (frontières de phases : SiO2/SiO, SiC/C(gr), SiC/Si(L)).
     Chaque point est donc calculé séparément avec un délai maximal, puis
     vérifié : on déduit les potentiels élémentaires lambda de la phase gaz,
     et toute phase condensée absente doit avoir une force motrice
     g°/RT - sum(a.lambda) >= 0 (polynômes NASA-9 lus dans nasa9.dat).
     Si le solveur bloque ou si le résultat n'est pas stable, on énumère
     les assemblages de phases condensées possibles (bprime sur un mélange
     restreint) et on retient celui qui satisfait le critère de stabilité :
     c'est l'équilibre vrai (minimum de l'enthalpie libre).
     Aux points où deux ou trois solides coexistent (SiC + C(gr),
     SiC + SiO2(L) + Si(L)...), le solveur de Mutation++ échoue même sur le
     mélange restreint. L'équilibre à assemblage fixé est alors résolu ici
     directement en potentiels élémentaires (mêmes données NASA-9, mêmes
     masses atomiques, même excès de char que bprime) ; ce solveur reproduit
     bprime à 1e-5 près là où celui-ci converge.

Usage :
    python sic_bprime.py

Prérequis :
    - Le binaire `bprime` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/sic-air.xml doit exister
    - MPP_DATA_DIRECTORY pointe sur data/ (sinon ../data est utilisé)
    - matplotlib, numpy et scipy installés
"""

import subprocess
import shutil
import sys
import os
import csv
import math
import itertools
import tempfile
import re
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.optimize import root, nnls
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paramètres du calcul
# ---------------------------------------------------------------------------
BPRIME_CMD = "bprime"           # nom du binaire (ou chemin absolu)
T_LIST     = np.arange(300.0, 5000.0 + 1.0e-6, 25.0)   # Tw de 300 à 5000 K, pas de 25 K
BG         = "0.0"              # débit pyrolyse adimensionné (0 = SiC pur)
MIXTURE    = "sic-air"
BL_COMP    = "air"
CHAR_COMP  = "sic"              # composition élémentaire du char (Si:1, C:1)
CHAR_ELEM  = "N"                # élément de suivi pour B'c (voir en-tête)

ONEATM = 101325.0   # Pa (aussi pression de référence des polynômes NASA-9)

# Pressions : 0.001 -> 1000 atm, espacement logarithmique (25 isobares)
PRESSURES_ATM = np.logspace(-3, 3, 25)

TIMEOUT_S  = 2.0     # délai maximal d'un appel à bprime pour un point
STAB_TOL   = 1.0e-4  # tolérance sur la force motrice des phases absentes
N_WORKERS  = os.cpu_count() or 1

# Familles de phases condensées : à une température donnée, un seul
# polymorphe de chaque famille est valide (plages NASA-9 disjointes).
CONDENSED_FAMILIES = [
    ["C(gr)"],
    ["Si(cr)", "Si(L)"],
    ["SiC(b)", "SiC(L)"],
    ["SiO2(a-qz)", "SiO2(b-qz)", "SiO2(b-crt)", "SiO2(L)"],
]
ELEMENTS = ["Si", "C", "N", "O"]
# Masses atomiques de data/thermo/elements.xml (celles de Mutation++)
ATOMIC_MW = {"Si": 28.085, "C": 12.011, "N": 14.0067, "O": 15.9994}
RU = 6.0221415e23 * 1.3806503e-23   # J/mol/K, src/general/Constants.h
AIR_X = {"N": 0.79, "O": 0.21}      # composition "air" de sic-air.xml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.environ.get("MPP_DATA_DIRECTORY",
                          os.path.join(SCRIPT_DIR, "..", "data"))


# ---------------------------------------------------------------------------
# Fonctions utilitaires
# ---------------------------------------------------------------------------

def find_bprime():
    """Localise le binaire bprime (PATH ou répertoires build courants)."""
    cmd = shutil.which(BPRIME_CMD)
    if cmd:
        return cmd
    candidates = [
        os.path.join(SCRIPT_DIR, "../build/src/apps/bprime"),
        os.path.join(SCRIPT_DIR, "../../build/src/apps/bprime"),
        "build/src/apps/bprime",
        "../build/src/apps/bprime",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def mixture_species(name):
    """Liste des espèces du fichier data/mixtures/<name>.xml."""
    txt = open(os.path.join(DATA_DIR, "mixtures", name + ".xml")).read()
    txt = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
    return re.search(r"<species>(.*?)</species>", txt, re.S).group(1).split()


def write_variant(tmpdir, name, gas, condensed):
    """Écrit un mélange sic-air restreint à un assemblage de phases donné."""
    with open(os.path.join(tmpdir, name + ".xml"), "w") as f:
        f.write('<mixture thermo_db="NASA-9">\n'
                f'    <species> {" ".join(gas + list(condensed))} </species>\n'
                '    <element_compositions default="air">\n'
                '        <composition name="air">N:0.79, O:0.21</composition>\n'
                '        <composition name="sic">Si:1.0, C:1.0</composition>\n'
                '    </element_compositions>\n'
                '</mixture>\n')


def run_bprime_sic(bprime_path, T, pressure_pa, mixture=MIXTURE, cwd=None):
    """
    Exécute bprime (char SiC) en un point (T, P).
    Retourne (header, row) ou None si le solveur ne rend pas la main.
    """
    cmd = [
        bprime_path,
        "-T", f"{T:g}",
        "-P", repr(float(pressure_pa)),
        "-b", BG,
        "-m", mixture,
        "-bl", BL_COMP,
        "-char", CHAR_COMP,
        "-char-elem", CHAR_ELEM,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=TIMEOUT_S, cwd=cwd)
    except subprocess.TimeoutExpired:
        return None
    if result.returncode != 0:
        print(f"\nERREUR bprime à T = {T} K, P = {pressure_pa:.2f} Pa :")
        print(result.stderr)
        sys.exit(1)
    header, data = parse_output(result.stdout)
    if len(data) != 1:
        return None
    return header, data[0]


def parse_output(output):
    """
    Parse la sortie de bprime.
    La première ligne est l'en-tête avec les noms entre guillemets.
    Retourne (header: list[str], data: np.ndarray).
    """
    lines = [l.strip() for l in output.strip().splitlines() if l.strip()]
    header = [h.strip('"') for h in lines[0].split()]
    data = []
    for line in lines[1:]:
        try:
            data.append([float(v) for v in line.split()])
        except ValueError:
            continue
    return header, np.array(data)


# ---------------------------------------------------------------------------
# Thermodynamique NASA-9 et critère de stabilité des phases
# ---------------------------------------------------------------------------

def load_nasa9(path, names):
    """Lit les polynômes NASA-9 (format de nasa9.dat) des espèces demandées."""
    want = set(names)
    lines = open(path).read().splitlines()
    db = {}
    i = 0
    while i < len(lines):
        ln = lines[i]
        name = ln[:24].split(" ")[0] if ln and ln[0] not in " !#" else None
        if name in want and i + 1 < len(lines) and lines[i + 1][:2].strip().isdigit():
            l2 = lines[i + 1]
            nr = int(l2[:2])
            form = {}
            for k in range(5):
                el = l2[10 + 8 * k:12 + 8 * k].strip().capitalize()
                n = float(l2[12 + 8 * k:18 + 8 * k] or 0.0)
                if el and n:
                    form[el] = form.get(el, 0.0) + n
            tb, co = [], []
            for r in range(nr):
                h = lines[i + 2 + 3 * r]
                a1 = lines[i + 3 + 3 * r].replace("D", "E")
                a2 = lines[i + 4 + 3 * r].replace("D", "E")
                if r == 0:
                    tb.append(float(h[1:11]))
                tb.append(float(h[11:21]))
                co.append([float(a1[16 * k:16 * k + 16]) for k in range(5)]
                          + [float(a2[0:16]), float(a2[16:32]),
                             float(a2[48:64]), float(a2[64:80])])
            db[name] = dict(tb=tb, co=co, cond=int(l2[50:52] or 0) != 0,
                            form=form)
            i += 2 + 3 * nr
            continue
        i += 1
    missing = want - set(db)
    if missing:
        sys.exit(f"Espèces absentes de {path} : {sorted(missing)}")
    return db


def g_over_rt(sp, T):
    """g°/RT à la pression de référence, formule de Nasa9Polynomial::gibbs."""
    tb, co = sp["tb"], sp["co"]
    tr = len(co) - 1
    for k in range(1, len(co)):
        if T < tb[k]:
            tr = k - 1
            break
    c = co[tr]
    lt = math.log(T)
    p = [-0.5 / T**2, (lt + 1.0) / T, 1.0 - lt, -0.5 * T,
         -T**2 / 6.0, -T**3 / 12.0, -T**4 / 20.0, 1.0 / T]
    return -c[8] + sum(ci * pi for ci, pi in zip(c[:8], p))


def valid_at(sp, T):
    """Même critère que NasaDB::speciesThermoValidAtT."""
    return sp["tb"][0] < T <= sp["tb"][-1]


def h_over_rt(sp, T):
    """h/RT, formule de Nasa9Polynomial::enthalpy."""
    tb, co = sp["tb"], sp["co"]
    tr = len(co) - 1
    for k in range(1, len(co)):
        if T < tb[k]:
            tr = k - 1
            break
    c = co[tr]
    p = [-1.0 / T**2, math.log(T) / T, 1.0, 0.5 * T,
         T**2 / 3.0, 0.25 * T**3, T**4 / 5.0, 1.0 / T]
    return sum(ci * pi for ci, pi in zip(c[:8], p))


def element_potentials(db, T, P, X):
    """
    Potentiels élémentaires lambda déduits de la phase gaz (moindres carrés
    sur mu_i/RT = g°_i/RT + ln(x_i P/P°) = sum_j a_ij lambda_j).
    """
    A, b = [], []
    for s, x in X.items():
        if not db[s]["cond"] and x > 1.0e-250:
            A.append([db[s]["form"].get(e, 0.0) for e in ELEMENTS])
            b.append(g_over_rt(db[s], T) + math.log(x) + math.log(P / ONEATM))
    return np.linalg.lstsq(np.array(A), np.array(b), rcond=None)[0]


def driving_forces(db, T, P, X):
    """
    Force motrice d_k = g°_k/RT - sum_j a_kj lambda_j de chaque phase
    condensée valide à T. À l'équilibre : d_k = 0 si k est présente,
    d_k >= 0 sinon.
    """
    lam = element_potentials(db, T, P, X)
    return {s: g_over_rt(db[s], T)
               - sum(db[s]["form"].get(e, 0.0) * l for e, l in zip(ELEMENTS, lam))
            for s in X if db[s]["cond"] and valid_at(db[s], T)}


def stability_violation(db, T, P, X):
    """Plus forte force motrice négative d'une phase condensée absente."""
    d = driving_forces(db, T, P, X)
    return max([0.0] + [-v for s, v in d.items() if X[s] <= 0.0])


def solve_fixed_assemblage(db, species, T, P, cond, lam0):
    """
    Équilibre à assemblage de phases condensées FIXÉ, pour le système de
    bprime : 1 kg d'air + 200 kg de SiC (surfaceMassBalance, B'g = 0).
    Inconnues : lambda (4), ln n_gaz, ln n_k (k dans cond).
    Équations : g°_k/RT = a_k.lambda (k dans cond), sum x_i = 1,
                bilans des 4 éléments.
    Retourne (row, X) au format de bprime, ou None si pas de solution
    (typiquement : une phase de l'assemblage devrait avoir n_k <= 0).
    """
    gas = [s for s in species if not db[s]["cond"]]
    Ag = np.array([[db[s]["form"].get(e, 0.0) for e in ELEMENTS] for s in gas])
    gg = np.array([g_over_rt(db[s], T) for s in gas]) + math.log(P / ONEATM)
    Ac = np.array([[db[s]["form"].get(e, 0.0) for e in ELEMENTS]
                   for s in cond]).reshape(len(cond), len(ELEMENTS))
    gc = np.array([g_over_rt(db[s], T) for s in cond])
    m_air = sum(AIR_X[e] * ATOMIC_MW[e] for e in AIR_X)
    ye = {e: AIR_X.get(e, 0.0) * ATOMIC_MW[e] / m_air for e in ELEMENTS}
    m_sic = ATOMIC_MW["Si"] + ATOMIC_MW["C"]
    yc = {"Si": ATOMIC_MW["Si"] / m_sic, "C": ATOMIC_MW["C"] / m_sic}
    char_amount = 200.0
    b = np.array([(ye[e] + char_amount * yc.get(e, 0.0)) / ATOMIC_MW[e]
                  for e in ELEMENTS])
    iN = ELEMENTS.index("N")

    def xg(lam):
        return np.exp(np.clip(Ag @ lam - gg, -700.0, 50.0))

    def F(u):
        lam, lng, lnc = u[:4], u[4], u[5:]
        x = xg(lam)
        el = (math.exp(min(lng, 700.0)) * (Ag.T @ x)
              + Ac.T @ np.exp(np.minimum(lnc, 700.0)))
        return np.r_[Ac @ lam - gc, x.sum() - 1.0, el / b - 1.0]

    # Estimations initiales : gaz de la graine, solides par NNLS sur les bilans
    x0 = xg(lam0)
    x0 = x0 / x0.sum()
    ng0 = b[iN] / max(Ag[:, iN] @ x0, 1.0e-300)
    guesses = [np.full(len(cond), math.log(b.max() * f)) for f in (0.5, 1e-3, 1e-8)]
    if cond:
        nc = nnls(Ac.T, np.maximum(b - ng0 * (Ag.T @ x0), 0.0))[0]
        guesses.insert(0, np.log(np.maximum(nc, 1.0e-10)))
    best = None
    with np.errstate(all="ignore"):
        for lc0 in guesses:
            u0 = np.r_[lam0, math.log(max(ng0, 1.0e-30)), lc0]
            sol = root(F, u0, method="hybr", options={"xtol": 1e-13, "maxfev": 20000})
            res = np.abs(F(sol.x)).max()
            if best is None or res < best[1]:
                best = (sol.x, res)
            if res < 1.0e-10:
                break
    if best[1] > 1.0e-9:
        return None

    x = xg(best[0][:4])
    x = x / x.sum()
    mw = np.array([sum(db[s]["form"].get(e, 0.0) * ATOMIC_MW[e] for e in ELEMENTS)
                   for s in gas])
    mwg = x @ mw
    ywN = ATOMIC_MW["N"] * (Ag[:, iN] @ x) / mwg
    Bc = max(ye["N"] / ywN - 1.0, 0.0)                       # -char-elem N
    hw = RU * T / (mwg * 1.0e-3) * (x @ np.array([h_over_rt(db[s], T) for s in gas]))
    X = {s: 0.0 for s in species}
    X.update(zip(gas, x))
    for s in cond:
        X[s] = 1.0                                           # IN_PHASE
    row = np.array([T, Bc, hw / 1.0e6] + [X[s] for s in species])
    return row, X


# ---------------------------------------------------------------------------
# Résolution d'un point (Tw, P)
# ---------------------------------------------------------------------------

def solve_point(ctx, T, P_pa):
    """
    Retourne (row, méthode, phases condensées présentes, violation).
    méthode : 'direct'     bprime sur sic-air complet
              'enum'       bprime sur un mélange restreint à un assemblage
              'assemblage' équilibre à assemblage fixé résolu en Python
              'approx'     aucun assemblage strictement stable trouvé : le
                           moins instable est retenu (ne se produit pas sur
                           la grille actuelle)
    """
    bprime_path, db, species, gas, tmpdir = ctx
    seeds = []
    res = run_bprime_sic(bprime_path, T, P_pa)
    if res is not None:
        header, row = res
        X = dict(zip(header[3:], row[3:]))
        viol = stability_violation(db, T, P_pa, X)
        if viol <= STAB_TOL:
            return row, "direct", [s for s in X if db[s]["cond"] and X[s] > 0], viol
        seeds.append(element_potentials(db, T, P_pa, X))

    # Énumération des assemblages : un polymorphe valide par famille,
    # de 0 à 4 familles présentes.
    fams = [[s for s in f if valid_at(db[s], T)] for f in CONDENSED_FAMILIES]
    fams = [f[0] for f in fams if f]
    subsets = [sub for n in range(len(fams) + 1)
               for sub in itertools.combinations(fams, n)]
    best = None
    for subset in subsets:
        if len(subset) == len(fams):
            continue       # identique au mélange complet, déjà essayé
        name = "sic_var_" + "_".join(
            re.sub(r"[()]", "", s) for s in subset) if subset else "sic_var_gaz"
        out = run_bprime_sic(bprime_path, T, P_pa, mixture=name, cwd=tmpdir)
        if out is None:
            continue
        h, r = out
        vals = dict(zip(h, r))
        row = np.array([vals.get(c, 0.0) for c in ["Tw[K]", "B'c", "hw[MJ/kg]"] + species])
        X = dict(zip(species, row[3:]))
        viol = stability_violation(db, T, P_pa, X)
        seeds.append(element_potentials(db, T, P_pa, X))
        if best is None or viol < best[3]:
            best = (row, "enum", [s for s in subset if X[s] > 0], viol)
        if viol <= STAB_TOL:
            return best

    # Coexistence de phases que Mutation++ ne sait pas converger : équilibre
    # à assemblage fixé, graines lambda tirées des calculs précédents.
    for subset in subsets:
        for lam0 in seeds:
            out = solve_fixed_assemblage(db, species, T, P_pa, list(subset), lam0)
            if out is None:
                continue
            row, X = out
            viol = stability_violation(db, T, P_pa, X)
            if best is None or viol < best[3]:
                best = (row, "assemblage", list(subset), viol)
            if viol <= STAB_TOL:
                return best
            break
    if best is None:
        sys.exit(f"Aucun assemblage calculable à T = {T} K, P = {P_pa} Pa")
    return (best[0], "approx", best[2], best[3])


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

# Pressions à tracer : uniquement les puissances de 10
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)   # 0.001, 0.01, 0.1, 1, 10, 100, 1000 atm


def plot_bprime_table(all_data, pressures_atm):
    """
    Trace la table B' SiC :
      - Gauche  : B'c vs Tw en échelle log10, une courbe par puissance de 10
      - Droite  : h_w vs Tw, une courbe par puissance de 10
    Seules les 7 pressions puissances de 10 sont tracées.
    """
    data_map = {P: d for P, d in zip(pressures_atm, all_data)}

    n = len(PLOT_PRESSURES_ATM)
    colors = plt.get_cmap("plasma", n + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        r"Table B' — Carbure de silicium (SiC) dans l'air  "
        r"($B'_g = 0$,  $P \in [10^{-3},\,10^3]$ atm)",
        fontsize=13
    )

    for idx, P_atm in enumerate(PLOT_PRESSURES_ATM):
        closest = min(pressures_atm, key=lambda p: abs(np.log10(p) - np.log10(P_atm)))
        header, data = data_map[closest]

        Tw = data[:, 0]
        Bc = data[:, 1]
        hw = data[:, 2]

        exp = int(round(np.log10(P_atm)))
        lbl = rf"$10^{{{exp}}}$ atm" if exp != 0 else "1 atm"

        color = colors(idx)
        ax1.plot(Tw, Bc, color=color, lw=2, label=lbl)
        ax2.plot(Tw, hw, color=color, lw=2, label=lbl)

    ax1.set_yscale("log")
    ax1.set_ylim(1e-5, 1e3)
    ax1.set_xlabel("Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(r"Taux d'ablation SiC $B'_c$ (masse nette gazéifiée)")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    ax2.set_xlabel("Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"Enthalpie de paroi $h_w$")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    plt.tight_layout()
    out_png = "sic_bprime_table.png"
    plt.savefig(out_png, dpi=150)
    print(f"Figure sauvegardée : {out_png}")
    plt.show()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # 1. Localiser le binaire
    bprime_path = find_bprime()
    if bprime_path is None:
        print(
            f"Binaire '{BPRIME_CMD}' introuvable.\n"
            "Compilez MutationPP :\n"
            "  cmake -B build -DCMAKE_BUILD_TYPE=Release .\n"
            "  cmake --build build --target bprime\n"
            "puis ajoutez build/src/apps/ au PATH."
        )
        sys.exit(1)
    print(f"Binaire trouvé : {bprime_path}\n")

    species = mixture_species(MIXTURE)
    db = load_nasa9(os.path.join(DATA_DIR, "thermo", "nasa9.dat"), species)
    gas = [s for s in species if not db[s]["cond"]]
    header = ["Tw[K]", "B'c", "hw[MJ/kg]"] + species

    # 2. Calcul point par point pour chaque pression
    all_data = []
    diag = []
    with tempfile.TemporaryDirectory() as tmpdir:
        # mélanges restreints (un par assemblage) pour l'énumération
        for fams in itertools.product(*[[None] + f for f in CONDENSED_FAMILIES]):
            subset = [s for s in fams if s]
            name = "sic_var_" + "_".join(
                re.sub(r"[()]", "", s) for s in subset) if subset else "sic_var_gaz"
            write_variant(tmpdir, name, gas, subset)
        ctx = (bprime_path, db, species, gas, tmpdir)

        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            for P_atm in PRESSURES_ATM:
                P_pa = P_atm * ONEATM
                print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
                      end=" ", flush=True)
                results = list(pool.map(lambda T: solve_point(ctx, T, P_pa), T_LIST))
                data = np.array([r[0] for r in results])
                all_data.append((header, data))
                for T, (row, meth, present, viol) in zip(T_LIST, results):
                    diag.append((T, P_atm, meth, present, viol))
                n_enum = sum(r[1] != "direct" for r in results)
                print(f"{len(data)} points ({n_enum} par énumération des phases)")

    counts = {m: sum(d[2] == m for d in diag)
              for m in ("direct", "enum", "assemblage", "approx")}
    print(f"\nPoints : {counts['direct']} directs, {counts['enum']} par "
          f"énumération bprime, {counts['assemblage']} par assemblage fixé, "
          f"{counts['approx']} approchés (violation max "
          f"{max([d[4] for d in diag] + [0.0]):.2e})")

    # 3. Sauvegarde CSV globale
    out_csv = "sic_bprime_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["P_atm"] + header)
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            for row in data:
                writer.writerow(
                    [f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row]
                )
    print(f"\nTable complète sauvegardée : {out_csv}")

    # 3bis. Table B'c au format long, unites SI : Tw_K, P_bar, Bc
    #       (nTw x nP ; B'g = 0 fixe, pas de colonne Bg)
    out_csv = "sic_bprime_bc_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Tw_K", "P_bar", "Bc"])
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            P_bar = P_atm * ONEATM / 1.0e5
            for row in data:
                writer.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}",
                                  f"{row[1]:.6e}"])
    print(f"Table B'c (nTw x nP) sauvegardée : {out_csv}")

    # 3ter. Table h_w au format long, unites SI : Tw_K, P_bar, hw_Jkg (nTw x nP)
    out_csv = "sic_bprime_hw_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Tw_K", "P_bar", "hw_Jkg"])
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            P_bar = P_atm * ONEATM / 1.0e5
            for row in data:
                writer.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}",
                                  f"{row[2] * 1.0e6:.6e}"])
    print(f"Table h_w (nTw x nP) sauvegardée : {out_csv}")

    # 3quater. Assemblage de phases condensées retenu en chaque point
    out_csv = "sic_bprime_phases.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Tw_K", "P_atm", "phases_condensees", "methode",
                         "violation_stabilite"])
        for T, P_atm, meth, present, viol in sorted(diag, key=lambda d: (d[1], d[0])):
            writer.writerow([f"{T:g}", f"{P_atm:.6g}", " ".join(present) or "-",
                             meth, f"{viol:.2e}"])
    print(f"Assemblages de phases sauvegardés : {out_csv}")

    # 4. Visualisation
    plot_bprime_table(all_data, PRESSURES_ATM)
