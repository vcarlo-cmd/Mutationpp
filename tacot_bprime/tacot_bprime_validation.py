#!/usr/bin/env python
"""
Génère la table B' du TACOT sur la grille EXACTE de la table de référence
distribuée avec le classeur TACOT 3.0, puis compare les deux.

Référence : feuille "B' table" du classeur TACOT_3.0_1.xls
    - calculée initialement avec Mutation (base CEA), puis avec TARGET
      [Bianchi, CRAS-TTG-1001, Sapienza, 2013]
    - hypothèses : coefficients de diffusion égaux, équilibre chimique à la
      paroi, air O2/N2 = 0.21/0.79 molaire
    - grille : 4 pressions x 25 valeurs de B'g x 151 températures = 15100 points

Grille de référence :
    p    : 0.001, 0.01, 0.1, 1.0 atm   (stockées en bar dans le classeur)
    B'g  : 0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4, 0.5, 0.6,
           0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0, 4.0, 5.5, 7.5, 10.0
    T    : 250 -> 4000 K, pas 25 K

Usage :
    python tacot_bprime_validation.py [--mixture tacot-air_25] [--ref chemin.xls]

Prérequis :
    - binaire `bprime` dans le PATH ou dans build/src/apps/
    - xlrd, numpy, matplotlib
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Grille de référence (extraite de la feuille "B' table")
# ---------------------------------------------------------------------------
ONEATM = 101325.0  # Pa

REF_PRESSURES_ATM = [0.001, 0.01, 0.1, 1.0]

REF_BG = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4,
          0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4,
          3.0, 4.0, 5.5, 7.5, 10.0]

T_RANGE = "250:25:4000"

BL_COMP   = "air"
PYRO_COMP = "tacot_pyro"
CHAR_COMP = "tacot_char"
CHAR_ELEM = "C"

DEFAULT_XLS = os.path.expanduser(
    "~/.claude/uploads/214c7635-1e0d-5852-968d-a13e109a5f11/f011c908-TACOT_3.0_1.xls")


# ---------------------------------------------------------------------------
# Binaire bprime
# ---------------------------------------------------------------------------

def find_bprime():
    """Localise le binaire bprime (PATH ou répertoires build courants)."""
    cmd = shutil.which("bprime")
    if cmd:
        return cmd
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, "../build/src/apps/bprime"),
              os.path.join(here, "../../build/src/apps/bprime"),
              "build/src/apps/bprime"):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def make_env():
    """Environnement d'exécution avec MPP_DATA_DIRECTORY défini."""
    env = os.environ.copy()
    env.setdefault("MPP_DATA_DIRECTORY",
                   os.path.abspath(os.path.join(
                       os.path.dirname(os.path.abspath(__file__)), "../data")))
    return env


def run_bprime(path, mixture, pressure_pa, bg):
    """Exécute bprime et retourne {T: (B'c, hw_MJ)}."""
    cmd = [path, "-T", T_RANGE, "-P", str(pressure_pa), "-b", str(bg),
           "-m", mixture, "-bl", BL_COMP, "-py", PYRO_COMP,
           "-char", CHAR_COMP, "-char-elem", CHAR_ELEM]
    res = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if res.returncode != 0:
        print(f"\nERREUR bprime (P={pressure_pa:.2f} Pa, B'g={bg}) :\n{res.stderr}")
        sys.exit(1)
    out = {}
    for line in res.stdout.strip().splitlines()[1:]:
        parts = line.split()
        try:
            T, bc, hw = float(parts[0]), float(parts[1]), float(parts[2])
        except (ValueError, IndexError):
            continue
        out[round(T)] = (bc, hw)
    return out


# ---------------------------------------------------------------------------
# Lecture de la table de référence
# ---------------------------------------------------------------------------

def load_reference(xls_path):
    """
    Lit la feuille "B' table" du classeur TACOT.
    Retourne {(P_atm_arrondi, B'g, T): (B'c, hw_MJ)}.
    Le classeur stocke la pression en bar et hw en kJ/kg.
    """
    import xlrd
    wb = xlrd.open_workbook(xls_path, on_demand=True)
    sh = wb.sheet_by_name("B' table")

    def num(v):
        # Le classeur utilise la notation Fortran '1.01325D-03'
        if isinstance(v, str):
            v = v.strip().replace("D", "E").replace("d", "E")
            if not v:
                return None
            try:
                return float(v)
            except ValueError:
                return None
        return v if isinstance(v, float) else None

    ref = {}
    for r in range(9, sh.nrows):
        p, bg, T, bc, hw = (num(sh.cell_value(r, c)) for c in range(5))
        if None in (p, bg, T, bc, hw):
            continue
        p_atm = round(p / 1.01325, 6)          # bar -> atm
        ref[(p_atm, round(bg, 6), round(T))] = (bc, hw / 1000.0)  # kJ -> MJ
    return ref


# ---------------------------------------------------------------------------
# Statistiques d'écart
# ---------------------------------------------------------------------------

def relative_error(computed, reference, floor):
    """Écart relatif, plancher sur |reference| pour éviter la division par ~0."""
    denom = max(abs(reference), floor)
    return (computed - reference) / denom * 100.0


def frozen_from(series):
    """
    Détecte le palier terminal de la table de référence.

    Au-delà de la sublimation complète du char, il ne reste plus de phase
    condensée : la table de référence (TARGET) gèle alors la dernière valeur
    convergée de B'c ET de hw jusqu'à la fin de la plage. Mutation++ continue
    au contraire à résoudre l'équilibre. Ces points ne sont donc pas
    comparables et sont isolés dans les statistiques.

    `series` : {T: (B'c, hw)}. Retourne la température de début de palier,
    ou None s'il n'y en a pas.
    """
    Ts = sorted(series)
    if len(Ts) < 4:
        return None
    last_bc, last_hw = series[Ts[-1]]
    if last_bc == 0.0:              # palier bas (B'c nul), traité à part
        return None
    start = Ts[-1]
    for T in reversed(Ts[:-1]):
        bc, hw = series[T]
        if abs(bc - last_bc) > 1e-9 * max(abs(last_bc), 1e-30) or \
           abs(hw - last_hw) > 1e-9 * max(abs(last_hw), 1e-30):
            break
        start = T
    # au moins 3 points identiques pour parler de palier
    return start if sum(1 for T in Ts if T >= start) >= 3 else None


def summarize(label, errs):
    a = np.abs(np.asarray(errs))
    if a.size == 0:
        return f"{label:<28s} (aucun point)"
    return (f"{label:<28s} moy {a.mean():7.3f} %   "
            f"med {np.median(a):7.3f} %   p95 {np.percentile(a,95):8.3f} %   "
            f"max {a.max():9.3f} %   n={a.size}")


# ---------------------------------------------------------------------------
# Graphiques
# ---------------------------------------------------------------------------

def plot_overlay(results, ref, mixture, out_png):
    """Superpose calcul et référence pour quatre B'g représentatifs."""
    show_bg = [0.0, 0.1, 1.0, 10.0]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Table B' TACOT — Mutation++ ({mixture}) vs référence TACOT 3.0",
                 fontsize=14)
    colors = plt.get_cmap("plasma", len(REF_PRESSURES_ATM) + 1)

    for ax, bg in zip(axes.ravel(), show_bg):
        for i, P in enumerate(REF_PRESSURES_ATM):
            Ts = sorted(T for (p, b, T) in ref if p == P and b == bg)
            if not Ts:
                continue
            r = [ref[(P, bg, T)][0] for T in Ts]
            c = [results[(P, bg)][T][0] for T in Ts if T in results[(P, bg)]]
            Tc = [T for T in Ts if T in results[(P, bg)]]
            exp = int(round(np.log10(P)))
            lbl = rf"$10^{{{exp}}}$ atm" if exp else "1 atm"
            ax.plot(Ts, r, color=colors(i), lw=3.0, alpha=0.35,
                    label=f"réf {lbl}")
            ax.plot(Tc, c, color=colors(i), lw=1.3, ls="--",
                    label=f"M++ {lbl}")
        ax.set_yscale("log")
        ax.set_title(rf"$B'_g = {bg}$")
        ax.set_xlabel(r"$T_w$ [K]")
        ax.set_ylabel(r"$B'_c$")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=7, ncol=2, loc="upper left")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"  Figure sauvegardée : {out_png}")


def plot_error_map(results, ref, mixture, plateau_start, out_png):
    """
    Carte de l'écart relatif sur B'c, un panneau par pression.

    Les points du palier terminal de la référence sont masqués (hachures) :
    la référence y est gelée et Mutation++ y sature sur sa propre limite de
    char ajouté, les deux valeurs ne sont donc pas comparables.
    """
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.8))
    fig.suptitle(rf"Écart relatif sur $B'_c$ — Mutation++ ({mixture}) vs référence "
                 r"(zone grise : palier de la référence, non comparable)",
                 fontsize=13)
    Ts = sorted(set(T for (_, _, T) in ref))

    for ax, P in zip(axes, REF_PRESSURES_ATM):
        M = np.full((len(REF_BG), len(Ts)), np.nan)
        for i, bg in enumerate(REF_BG):
            Tp = plateau_start.get((P, bg))
            for j, T in enumerate(Ts):
                key = (P, bg, T)
                if key not in ref or T not in results.get((P, bg), {}):
                    continue
                if Tp is not None and T >= Tp:
                    continue                      # palier -> masqué
                M[i, j] = relative_error(results[(P, bg)][T][0],
                                         ref[key][0], 1e-3)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("0.85")                      # gris pour les points masqués
        im = ax.pcolormesh(Ts, REF_BG, np.ma.masked_invalid(M), cmap=cmap,
                           vmin=-5, vmax=5, shading="nearest")
        exp = int(round(np.log10(P)))
        ax.set_title(rf"$10^{{{exp}}}$ atm" if exp else "1 atm")
        ax.set_xlabel(r"$T_w$ [K]")
        ax.set_yscale("log")
        ax.set_ylabel(r"$B'_g$")
        fig.colorbar(im, ax=ax, label="écart [%]", extend="both")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"  Figure sauvegardée : {out_png}")


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mixture", default="tacot-air_25",
                    help="nom du mélange (défaut : tacot-air_25)")
    ap.add_argument("--ref", default=DEFAULT_XLS,
                    help="chemin du classeur TACOT 3.0")
    args = ap.parse_args()

    path = find_bprime()
    if path is None:
        print("Binaire 'bprime' introuvable. Compilez :\n"
              "  cmake -B build -DCMAKE_BUILD_TYPE=Release .\n"
              "  cmake --build build --target bprime")
        sys.exit(1)
    print(f"Binaire   : {path}")
    print(f"Mélange   : {args.mixture}")
    print(f"Référence : {args.ref}\n")

    ref = load_reference(args.ref)
    print(f"Table de référence : {len(ref)} points\n")

    # --- Calcul sur la grille de référence --------------------------------
    results = {}
    total = len(REF_PRESSURES_ATM) * len(REF_BG)
    n = 0
    for P in REF_PRESSURES_ATM:
        for bg in REF_BG:
            n += 1
            print(f"\r  calcul {n:3d}/{total}  (P={P:g} atm, B'g={bg:g})   ",
                  end="", flush=True)
            results[(P, bg)] = run_bprime(path, args.mixture, P * ONEATM, bg)
    print("\n")

    # --- Classement des points par régime ---------------------------------
    # ACTIF  : référence non nulle et hors palier terminal -> comparable
    # ZERO   : référence exactement nulle (char non consommé)
    # PALIER : référence gelée après sublimation complète -> non comparable
    plateau_start = {}
    for P in REF_PRESSURES_ATM:
        for bg in REF_BG:
            series = {T: ref[(P, bg, T)] for T in
                      (t for (p, b, t) in ref if p == P and b == bg)}
            plateau_start[(P, bg)] = frozen_from(series)

    out_csv = f"tacot_bprime_vs_ref_{args.mixture}.csv"
    err_bc, err_hw = [], []          # régime ACTIF uniquement
    err_hw_all = []                  # hw hors palier (y compris B'c nul)
    per_bg = {bg: [] for bg in REF_BG}
    per_p = {P: [] for P in REF_PRESSURES_ATM}
    n_zero = n_zero_ok = n_plateau = 0

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "Bg", "T_K", "regime",
                    "Bc_mpp", "Bc_ref", "Bc_err_pct",
                    "hw_mpp_MJkg", "hw_ref_MJkg", "hw_err_pct"])
        for P in REF_PRESSURES_ATM:
            for bg in REF_BG:
                Tp = plateau_start[(P, bg)]
                for T in sorted(results[(P, bg)]):
                    key = (P, bg, T)
                    if key not in ref:
                        continue
                    bc_m, hw_m = results[(P, bg)][T]
                    bc_r, hw_r = ref[key]
                    e_bc = relative_error(bc_m, bc_r, 1e-3)
                    e_hw = relative_error(hw_m, hw_r, 1e-1)

                    if Tp is not None and T >= Tp:
                        regime = "PALIER"
                        n_plateau += 1
                    elif bc_r == 0.0:
                        regime = "ZERO"
                        n_zero += 1
                        if bc_m < 1e-6:
                            n_zero_ok += 1
                        err_hw_all.append(e_hw)
                    else:
                        regime = "ACTIF"
                        err_bc.append(e_bc)
                        err_hw.append(e_hw)
                        err_hw_all.append(e_hw)
                        per_bg[bg].append(e_bc)
                        per_p[P].append(e_bc)

                    w.writerow([f"{P:g}", f"{bg:g}", f"{T:g}", regime,
                                f"{bc_m:.8e}", f"{bc_r:.8e}", f"{e_bc:.4f}",
                                f"{hw_m:.8e}", f"{hw_r:.8e}", f"{e_hw:.4f}"])

    total_pts = len(err_bc) + n_zero + n_plateau
    print(f"Comparaison point par point : {out_csv}")
    print(f"Points de la grille         : {total_pts}")
    print(f"  ACTIF  (comparables)      : {len(err_bc)}")
    print(f"  ZERO   (B'c_ref = 0)      : {n_zero}")
    print(f"  PALIER (réf. gelée)       : {n_plateau}\n")

    # --- Statistiques -----------------------------------------------------
    print("=" * 92)
    print("RÉGIME ACTIF — écart relatif Mutation++ / référence")
    print("=" * 92)
    print(summarize("B'c", err_bc))
    print(summarize("hw", err_hw))

    within = lambda e, t: 100.0 * np.mean(np.abs(np.asarray(e)) < t)
    print(f"\nB'c : {within(err_bc,1):.2f} % des points < 1 %   "
          f"{within(err_bc,5):.2f} % < 5 %   {within(err_bc,10):.2f} % < 10 %")
    print(f"hw  : {within(err_hw,1):.2f} % des points < 1 %   "
          f"{within(err_hw,5):.2f} % < 5 %   {within(err_hw,10):.2f} % < 10 %")

    print("\n" + "=" * 92)
    print("RÉGIME ZERO — la référence impose B'c = 0")
    print("=" * 92)
    frac = 100.0 * n_zero_ok / n_zero if n_zero else float("nan")
    print(f"Mutation++ donne aussi B'c < 1e-6 sur {n_zero_ok}/{n_zero} "
          f"points ({frac:.2f} %)")
    print(summarize("hw (ZERO + ACTIF)", err_hw_all))

    print("\n" + "=" * 92)
    print("ÉCART SUR B'c PAR PRESSION (régime actif)")
    print("=" * 92)
    for P in REF_PRESSURES_ATM:
        print(summarize(f"P = {P:g} atm", per_p[P]))

    print("\n" + "=" * 92)
    print("ÉCART SUR B'c PAR B'g (régime actif)")
    print("=" * 92)
    for bg in REF_BG:
        print(summarize(f"B'g = {bg:g}", per_bg[bg]))

    print("\n" + "=" * 92)
    print("DÉBUT DU PALIER DE LA RÉFÉRENCE [K] (sublimation complète)")
    print("=" * 92)
    hdr = "  B'g   " + "".join(f"{P:>12g}" for P in REF_PRESSURES_ATM)
    print(hdr.replace("B'g", "B'g \\ P atm"))
    for bg in REF_BG:
        cells = "".join(
            f"{plateau_start[(P,bg)]:>12g}" if plateau_start[(P, bg)] else f"{'-':>12s}"
            for P in REF_PRESSURES_ATM)
        print(f"{bg:>6g}   {cells}")

    # --- Graphiques -------------------------------------------------------
    print()
    plot_overlay(results, ref, args.mixture,
                 f"tacot_bprime_vs_ref_{args.mixture}.png")
    plot_error_map(results, ref, args.mixture, plateau_start,
                   f"tacot_bprime_err_{args.mixture}.png")
    print("\nTerminé.")


if __name__ == "__main__":
    main()
