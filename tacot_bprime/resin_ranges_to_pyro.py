#!/usr/bin/env python
"""
Des fractions MASSIQUES de la résine (données en plages) à la composition
du gaz de pyrolyse à renseigner dans le XML.

Deux étapes, souvent confondues :

  1. massique -> molaire        x_E = (y_E/M_E) / SUM_j (y_j/M_j)
     (facultatif : on peut aussi déclarer type="mass" dans le XML et
      laisser Mutation++ convertir)

  2. résine -> GAZ DE PYROLYSE  n_E(gaz) = n_E(résine) - n_E(char)
     C'est cette étape qui compte : le XML attend le GAZ, pas la résine.

Quand la composition de la résine est connue par plages, ce script
explore le domaine admissible (les fractions doivent sommer à 1), propage
la fermeture, et donne l'enveloppe résultante sur le gaz. Il évalue
ensuite l'impact réel sur B'c et h_w en appelant `bprime`, de façon à
savoir si l'incertitude est gênante ou négligeable.

Hypothèse sur l'azote : à haute température l'azote de la résine part
essentiellement avec les gaz (NH3, HCN, N2). Le char est donc supposé
carbone pur, et TOUT l'azote passe dans le gaz. Option --n-in-char pour
en retenir une fraction dans le char.

Usage :
    python resin_ranges_to_pyro.py
    python resin_ranges_to_pyro.py --char-yield 0.55 --n-in-char 0.3
"""

import argparse
import itertools
import os
import shutil
import subprocess
import sys

import numpy as np

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}
ELEMENTS = ["C", "H", "O", "N"]

# Plages de fractions MASSIQUES de la résine (données utilisateur)
RANGES = {
    "C": (0.71, 0.78),
    "O": (0.16, 0.18),
    "H": (0.05, 0.06),
    "N": (0.00, 0.07),
}

ONEATM = 101325.0


# ---------------------------------------------------------------------------
# Domaine admissible
# ---------------------------------------------------------------------------

def feasible_set(ranges, n=41):
    """
    Grille sur C, H, O ; N est imposé par la fermeture SUM y = 1.
    Ne garde que les points dont les quatre fractions sont dans leurs plages.
    Retourne un tableau (m, 4) de fractions massiques.
    """
    axes = [np.linspace(*ranges[e], n) for e in ("C", "H", "O")]
    keep = []
    lo_n, hi_n = ranges["N"]
    for c, h, o in itertools.product(*axes):
        nn = 1.0 - c - h - o
        if lo_n - 1e-12 <= nn <= hi_n + 1e-12:
            keep.append((c, h, o, max(nn, 0.0)))
    return np.array(keep)


def mass_to_mole(y):
    """Fractions massiques -> fractions molaires (élémentaires)."""
    n = np.array([y[i] / M[e] for i, e in enumerate(ELEMENTS)])
    return n / n.sum()


# ---------------------------------------------------------------------------
# Fermeture résine -> gaz
# ---------------------------------------------------------------------------

def pyro_gas(y_resin, char_yield, n_in_char=0.0):
    """
    y_resin    : fractions MASSIQUES de la résine, ordre ELEMENTS
    char_yield : fraction massique de résine restant en char
    n_in_char  : fraction de l'azote de la résine retenue dans le char

    Retourne les fractions MOLAIRES élémentaires du gaz, ou None si la
    fermeture est impossible (char plus riche en C que la résine).
    """
    basis = 100.0                                   # g de résine
    m = {e: basis * y_resin[i] for i, e in enumerate(ELEMENTS)}
    n_res = {e: m[e] / M[e] for e in ELEMENTS}

    m_char = char_yield * basis
    # azote retenu dans le char
    m_n_char = min(n_in_char * m["N"], m_char)
    m_c_char = m_char - m_n_char
    if m_c_char > m["C"] + 1e-9:
        return None                                  # pas assez de carbone

    n_gas = {
        "C": n_res["C"] - m_c_char / M["C"],
        "H": n_res["H"],
        "O": n_res["O"],
        "N": n_res["N"] - m_n_char / M["N"],
    }
    if any(v < -1e-9 for v in n_gas.values()):
        return None
    tot = sum(max(v, 0.0) for v in n_gas.values())
    if tot <= 0:
        return None
    return np.array([max(n_gas[e], 0.0) / tot for e in ELEMENTS])


# ---------------------------------------------------------------------------
# bprime
# ---------------------------------------------------------------------------

def find_bprime():
    c = shutil.which("bprime")
    if c:
        return c
    here = os.path.dirname(os.path.abspath(__file__))
    for p in (os.path.join(here, "../build/src/apps/bprime"),
              "build/src/apps/bprime"):
        if os.path.isfile(p):
            return os.path.abspath(p)
    return None


def make_env():
    env = os.environ.copy()
    env.setdefault("MPP_DATA_DIRECTORY", os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../data")))
    return env


def run_case(path, mixdir, gas_mole, bg, pressure_pa, temps):
    """
    Écrit un mélange temporaire portant la composition `gas_mole` puis
    appelle bprime. Retourne {T: (B'c, hw)}.
    """
    name = "tmp_resin_uq"
    comp = ", ".join(f"{e}:{v:.8f}"
                     for e, v in zip(ELEMENTS, gas_mole) if v > 1e-12)
    xml = f"""<mixture thermo_db="NASA-9">
    <species>
       C H O N CH4 CN CO CO2 C2 C2H C2H2,acetylene C3 C4 C4H2,butadiyne C5
       HCN H2 H2O N2 CH2OH CNN CNC CNCOCN C6H6 HNC C(gr)
    </species>
    <element_compositions default="air">
        <composition name="air">N:0.79, O:0.21</composition>
        <composition name="pyro">{comp}</composition>
        <composition name="char">C:1.0</composition>
    </element_compositions>
</mixture>
"""
    with open(os.path.join(mixdir, "mixtures", f"{name}.xml"), "w") as f:
        f.write(xml)

    env = make_env()
    env["MPP_DATA_DIRECTORY"] = mixdir
    lo, hi = min(temps), max(temps)
    cmd = [path, "-T", f"{lo}:25:{hi}", "-P", str(pressure_pa), "-b", str(bg),
           "-m", name, "-bl", "air", "-py", "pyro",
           "-char", "char", "-char-elem", "C"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"\nERREUR bprime :\n{r.stderr}")
        sys.exit(1)
    out = {}
    for line in r.stdout.strip().splitlines()[1:]:
        p = line.split()
        try:
            out[round(float(p[0]))] = (float(p[1]), float(p[2]))
        except (ValueError, IndexError):
            continue
    return out


# ---------------------------------------------------------------------------
# Programme
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--char-yield", type=float, default=0.50)
    ap.add_argument("--n-in-char", type=float, default=0.0,
                    help="fraction de l'azote retenue dans le char")
    ap.add_argument("--no-bprime", action="store_true")
    args = ap.parse_args()

    print("=" * 78)
    print("1. DOMAINE ADMISSIBLE DE LA RÉSINE")
    print("=" * 78)
    print("  plages fournies (fractions massiques) :")
    for e in ELEMENTS:
        print(f"     {e} : {100*RANGES[e][0]:5.1f} – {100*RANGES[e][1]:5.1f} %")
    smin = sum(RANGES[e][0] for e in ELEMENTS)
    smax = sum(RANGES[e][1] for e in ELEMENTS)
    print(f"\n  somme des bornes basses : {100*smin:.0f} %")
    print(f"  somme des bornes hautes : {100*smax:.0f} %")
    print("  -> les plages ne sont donc PAS indépendantes : toute composition")
    print("     réelle doit sommer à 100 %. On explore l'intersection.")

    Y = feasible_set(RANGES)
    print(f"\n  points admissibles trouvés : {len(Y)}")
    print("  plages effectives une fois la contrainte appliquée :")
    for i, e in enumerate(ELEMENTS):
        print(f"     {e} : {100*Y[:,i].min():5.2f} – {100*Y[:,i].max():5.2f} %")

    # --- 2. massique -> molaire ------------------------------------------
    print("\n" + "=" * 78)
    print("2. RÉSINE : MASSIQUE -> MOLAIRE")
    print("=" * 78)
    ymid = np.array([Y[:, i].mean() for i in range(4)])
    ymid = ymid / ymid.sum()
    xmid = mass_to_mole(ymid)
    print("  composition médiane du domaine :")
    print("     % masse  : " + "  ".join(
        f"{e} {100*ymid[i]:5.2f}" for i, e in enumerate(ELEMENTS)))
    print("     x molaire: " + "  ".join(
        f"{e} {xmid[i]:6.4f}" for i, e in enumerate(ELEMENTS)))
    print("\n  (rappel : l'hydrogène pèse peu mais compte beaucoup en moles —")
    print("   ~5.5 % de la masse devient ~30 % des atomes)")

    # --- 3. fermeture -> gaz ---------------------------------------------
    print("\n" + "=" * 78)
    print(f"3. GAZ DE PYROLYSE — fermeture, rendement char "
          f"{100*args.char_yield:.0f} %, azote au char {100*args.n_in_char:.0f} %")
    print("=" * 78)
    G = [g for g in (pyro_gas(y, args.char_yield, args.n_in_char) for y in Y)
         if g is not None]
    if not G:
        print("  Aucune fermeture possible : rendement en char trop élevé "
              "pour ce domaine.")
        sys.exit(1)
    G = np.array(G)
    gmid = pyro_gas(ymid, args.char_yield, args.n_in_char)

    print(f"  fermetures réussies : {len(G)} / {len(Y)}")
    print("\n  ENVELOPPE du gaz (fractions molaires élémentaires) :")
    print(f"  {'':3s} {'min':>9} {'médiane':>9} {'max':>9}   {'étendue':>9}")
    for i, e in enumerate(ELEMENTS):
        lo, hi = G[:, i].min(), G[:, i].max()
        span = (hi - lo) / gmid[i] * 100 if gmid[i] > 1e-9 else float("nan")
        print(f"  {e:3s} {lo:9.4f} {gmid[i]:9.4f} {hi:9.4f}   {span:8.1f} %")

    print("\n  >>> À RENSEIGNER DANS LE XML (valeur médiane) :")
    comp = ", ".join(f"{e}:{gmid[i]:.4f}"
                     for i, e in enumerate(ELEMENTS) if gmid[i] > 1e-6)
    print(f'  <composition name="x_pyro">{comp}</composition>')

    # --- 4. impact sur la table B' ---------------------------------------
    if args.no_bprime:
        return
    path = find_bprime()
    if path is None:
        print("\n(binaire bprime introuvable : impact non évalué)")
        return

    print("\n" + "=" * 78)
    print("4. IMPACT RÉEL SUR LA TABLE B'  (1 atm, B'g = 0.5)")
    print("=" * 78)

    here = os.path.dirname(os.path.abspath(__file__))
    mixdir = os.path.abspath(os.path.join(here, "../data"))
    temps = [1000, 2000, 3000, 3500, 3800]

    # cas extrêmes : min et max de C dans le gaz, plus la médiane
    cases = {
        "médiane":  gmid,
        "C min":    G[int(np.argmin(G[:, 0]))],
        "C max":    G[int(np.argmax(G[:, 0]))],
        "N max":    G[int(np.argmax(G[:, 3]))],
    }
    res = {}
    tmpxml = os.path.join(mixdir, "mixtures", "tmp_resin_uq.xml")
    try:
        for label, g in cases.items():
            res[label] = run_case(path, mixdir, g, 0.5, ONEATM,
                                  [min(temps), max(temps)])
    finally:
        if os.path.exists(tmpxml):
            os.remove(tmpxml)

    print(f"  {'T [K]':>7} | " + " | ".join(f"{k:>19}" for k in cases))
    print(f"  {'':7s} | " + " | ".join(f"{'B''c':>9}{'hw':>10}"
                                       for _ in cases))
    for T in temps:
        row = f"  {T:>7d} | "
        row += " | ".join(
            f"{res[k][T][0]:>9.4f}{res[k][T][1]:>10.3f}"
            if T in res[k] else f"{'--':>19}" for k in cases)
        print(row)

    print("\n  écart maximal sur B'c entre les cas extrêmes :")
    for T in temps:
        vals = [res[k][T][0] for k in cases if T in res[k]]
        hws = [res[k][T][1] for k in cases if T in res[k]]
        if not vals or max(vals) < 1e-9:
            continue
        sp = (max(vals) - min(vals)) / max(np.mean(vals), 1e-12) * 100
        sph = (max(hws) - min(hws)) / max(abs(np.mean(hws)), 1e-12) * 100
        print(f"     T = {T:5d} K :  B'c {sp:6.2f} %      hw {sph:6.2f} %")

    # --- 5. quel élément pilote l'incertitude ? --------------------------
    print("\n" + "=" * 78)
    print("5. SENSIBILITÉ ÉLÉMENT PAR ÉLÉMENT  (un à la fois, autres au centre)")
    print("=" * 78)
    print("  Chaque élément est porté à sa borne basse puis haute, les trois")
    print("  autres gardant leurs proportions médianes (renormalisation à 1).")
    print(f"\n  {'élément':>8} | " +
          " | ".join(f"T={T} K".rjust(16) for T in (1000, 3000, 3500)))
    print(f"  {'':8s} | " + " | ".join(
        f"{'B''c min':>8}{'B''c max':>8}" for _ in (1, 2, 3)))

    ranking = {}
    try:
        for i, e in enumerate(ELEMENTS):
            vals = {}
            for bound in (0, 1):
                y = ymid.copy()
                y[i] = RANGES[e][bound]
                others = [j for j in range(4) if j != i]
                rest = 1.0 - y[i]
                s = ymid[others].sum()
                for j in others:
                    y[j] = ymid[j] / s * rest
                if any(v < 0 for v in y):
                    vals[bound] = None
                    continue
                g = pyro_gas(y, args.char_yield, args.n_in_char)
                vals[bound] = (run_case(path, mixdir, g, 0.5, ONEATM,
                                        [1000, 3500]) if g is not None else None)
            cells, spans = [], []
            for T in (1000, 3000, 3500):
                a = vals[0][T][0] if vals[0] and T in vals[0] else None
                b = vals[1][T][0] if vals[1] and T in vals[1] else None
                if a is None or b is None:
                    cells.append(f"{'--':>16}")
                    continue
                cells.append(f"{a:>8.4f}{b:>8.4f}")
                # métrique dégénérée si l'une des bornes annule B'c :
                # on ne classe que sur les points bien conditionnés
                if min(a, b) > 1e-6:
                    spans.append(abs(b - a) / np.mean([a, b]) * 100)
            print(f"  {e:>8s} | " + " | ".join(cells))
            ranking[e] = max(spans) if spans else float("nan")
    finally:
        if os.path.exists(tmpxml):
            os.remove(tmpxml)

    print("\n  Amplitude induite sur B'c (points où B'c ne s'annule pas) :")
    for e, v in sorted(ranking.items(),
                       key=lambda kv: -(kv[1] if kv[1] == kv[1] else -1)):
        if v != v:
            print(f"     {e} :      n/a  (B'c s'annule sur une borne)")
            continue
        print(f"     {e} : {v:7.1f} %  " + "#" * int(min(v, 100) / 2))
    valid = {e: v for e, v in ranking.items() if v == v}
    if valid:
        top = max(valid, key=valid.get)
        print(f"\n  -> l'élément à mesurer en priorité est {top}.")
    print("\n  Note : à 1000 K, porter C à sa borne haute annule B'c — le gaz")
    print("  de pyrolyse apporte alors assez de carbone pour que le char ne")
    print("  soit plus consommé. Ce basculement est le plus sensible de tous.")


if __name__ == "__main__":
    main()
