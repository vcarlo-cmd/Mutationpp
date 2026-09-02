#!/usr/bin/env python
"""
Enthalpie et propriétés du gaz de pyrolyse du TACOT, à l'équilibre chimique.

Calcule h_g(T, P) — ainsi que M, Cp, gamma, rho, mu — avec l'outil MutationPP
`mppequil`, et compare à la table de référence du classeur TACOT 3.0
(feuille 'Sheet9' : 4 pressions x 152 températures).

h_g est la grandeur qui manque à la table B' pour fermer le bilan d'énergie :

    - bilan de surface   : q_conv + q_rad = q_em + q_cond
                                          + [mdot_w*h_w - mdot_g*h_g - mdot_c*h_c]
      (h_w vient de la table B', h_g de ce calcul)

    - source en profondeur : (h_g - h_s) * d(rho)/dt      <- la "chaleur de pyrolyse"

Hypothèse du classeur, reprise ici : équilibre chimique du gaz, SANS phase
condensée. Le gaz réel qui sort de la résine est métastable et n'est pas à
l'équilibre à basse température — voir la discussion en fin de sortie.

Ce script COMPARE à la référence. Pour calculer h_g sans comparaison,
voir tacot_pyrolysis_gas.py.

Usage :
    python tacot_pyrolysis_gas_validation.py [--mixture tacot-pyrogas] [--ref chemin.xls]

Prérequis :
    - binaire `mppequil` (cmake --build build --target mppequil)
    - data/mixtures/tacot-pyrogas.xml
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

ONEATM = 101325.0

# Grille de la référence (feuille Sheet9) : 0.01, 0.1 et 1 atm.
#
# ATTENTION — Sheet9 contient QUATRE blocs, mais seulement TROIS pressions.
# Le bloc étiqueté "P = 1.01 Pa" est un DOUBLON bit-identique du bloc à
# 1013.25 Pa (vérifié : écart max 0.000e+00 sur les 152 températures et les
# 6 grandeurs). Son étiquette de pression est fausse d'un facteur 1000 —
# 1.01325 y est en kPa, pas en Pa. L'en-tête de la feuille 'Pyrolysis model'
# confirme : "Pyrolysis gas properties at 0.01, 0.1, 1.0 atm".
# Ce bloc est donc ignoré ; l'inclure produirait un faux écart de 99.9 %
# sur rho (rapport exactement 1/1000).
REF_PRESSURES_PA = [1013.25, 10132.5, 101325.0]   # 0.01, 0.1, 1 atm
T_RANGE = "200:25:3975"

# Grandeurs demandées à mppequil, dans l'ordre des colonnes de sortie.
# indices : 0 Th[K], 5 Mw[kg/mol], 9 Cp_eq[J/kg-K], 10 H[J/kg],
#           18 gam_eq[-], 3 rho[kg/m3], 32 mu[Pa-s]
MPP_INDICES = "0,5,9,10,18,3,32"

# Conversion vers les unités du classeur
#   nom interne : (colonne mppequil, facteur, unité classeur, colonne Sheet9)
QUANTITIES = [
    ("M",     1, 1e3,  "kg/kmol",    2),
    ("Cp",    2, 1e-3, "kJ/kg-K",    3),
    ("h",     3, 1e-3, "kJ/kg",      5),
    ("gamma", 4, 1.0,  "-",          4),
    ("rho",   5, 1.0,  "kg/m3",      7),
    ("mu",    6, 1e4,  "millipoise", 6),
]

DEFAULT_XLS = os.path.expanduser(
    "~/.claude/uploads/214c7635-1e0d-5852-968d-a13e109a5f11/f011c908-TACOT_3.0_1.xls")


# ---------------------------------------------------------------------------
# mppequil
# ---------------------------------------------------------------------------

def find_mppequil():
    cmd = shutil.which("mppequil")
    if cmd:
        return cmd
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, "../build/src/apps/mppequil"),
              os.path.join(here, "../../build/src/apps/mppequil"),
              "build/src/apps/mppequil"):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def make_env():
    env = os.environ.copy()
    env.setdefault("MPP_DATA_DIRECTORY", os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../data")))
    return env


def run_mppequil(path, mixture, pressure_pa):
    """Retourne {T_arrondi: [Mw, Cp, H, gamma, rho, mu]} en unités SI."""
    cmd = [path, "-T", T_RANGE, "-P", str(pressure_pa),
           "-m", MPP_INDICES, "--elem-comp", "tacot_pyro", mixture]
    r = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if r.returncode != 0:
        print(f"\nERREUR mppequil à P = {pressure_pa} Pa :\n{r.stderr}")
        sys.exit(1)
    out = {}
    for line in r.stdout.strip().splitlines()[1:]:
        parts = line.split()
        try:
            vals = [float(v) for v in parts]
        except ValueError:
            continue
        if len(vals) < 7:
            continue
        out[round(vals[0])] = vals
    return out


# ---------------------------------------------------------------------------
# Référence
# ---------------------------------------------------------------------------

def load_reference(xls_path):
    """{(P_pa_arrondi, T): [M, Cp, gamma, h, mu, rho]} depuis Sheet9."""
    import xlrd
    wb = xlrd.open_workbook(xls_path, on_demand=True)
    sh = wb.sheet_by_name("Sheet9")
    ref = {}
    for r in range(1, sh.nrows):
        try:
            P = float(sh.cell_value(r, 0))
            T = float(sh.cell_value(r, 1))
            row = [float(sh.cell_value(r, c)) for c in range(2, 8)]
        except (ValueError, TypeError):
            continue
        ref[(round(P, 2), round(T))] = row
    return ref


# ---------------------------------------------------------------------------
# Comparaison
# ---------------------------------------------------------------------------

def rel_err(a, b, floor):
    return (a - b) / max(abs(b), floor) * 100.0


def summarize(label, errs):
    a = np.abs(np.asarray(errs))
    if a.size == 0:
        return f"{label:<22s} (aucun point)"
    return (f"{label:<22s} moy {a.mean():7.3f} %  med {np.median(a):7.3f} %  "
            f"p95 {np.percentile(a,95):7.3f} %  max {a.max():8.3f} %  n={a.size}")


def plot_all(results, ref, out_png):
    fig, axes = plt.subplots(2, 3, figsize=(19, 10))
    fig.suptitle("Gaz de pyrolyse du TACOT à l'équilibre — Mutation++ (traits) "
                 "vs référence TACOT 3.0 (points)", fontsize=14)
    colors = plt.get_cmap("plasma", len(REF_PRESSURES_PA) + 1)
    ylabels = {"M": "M [kg/kmol]", "Cp": r"$C_p$ [kJ/kg/K]",
               "h": r"$h_g$ [kJ/kg]", "gamma": r"$\gamma$",
               "rho": r"$\rho$ [kg/m³]", "mu": r"$\mu$ [millipoise]"}

    for ax, (name, ci, fac, unit, refcol) in zip(axes.ravel(), QUANTITIES):
        for i, P in enumerate(REF_PRESSURES_PA):
            Ts = sorted(t for (p, t) in ref if abs(p - round(P, 2)) < 1e-6)
            if not Ts:
                continue
            rv = [ref[(round(P, 2), t)][refcol - 2] for t in Ts]
            mv = [results[P][t][ci] * fac for t in Ts if t in results[P]]
            mt = [t for t in Ts if t in results[P]]
            lbl = f"{P/ONEATM:g} atm"
            ax.plot(mt, mv, color=colors(i), lw=1.8, label=f"M++ {lbl}")
            ax.plot(Ts[::6], rv[::6], "o", color=colors(i), ms=3.5,
                    mec="k", mew=0.3, label=f"réf {lbl}")
        ax.set_xlabel(r"$T$ [K]")
        ax.set_ylabel(ylabels[name])
        ax.grid(True, ls="--", alpha=0.35)
        if name in ("rho", "mu"):
            ax.set_yscale("log")
        ax.legend(fontsize=6.5, ncol=2)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"\nFigure sauvegardée : {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mixture", default="tacot-pyrogas")
    ap.add_argument("--ref", default=DEFAULT_XLS)
    args = ap.parse_args()

    path = find_mppequil()
    if path is None:
        print("Binaire 'mppequil' introuvable. Compilez :\n"
              "  cmake --build build --target mppequil")
        sys.exit(1)
    print(f"Binaire  : {path}")
    print(f"Mélange  : {args.mixture}\n")

    results = {}
    for P in REF_PRESSURES_PA:
        print(f"  mppequil à P = {P/ONEATM:g} atm ...", end=" ", flush=True)
        results[P] = run_mppequil(path, args.mixture, P)
        print(f"{len(results[P])} points")

    ref = load_reference(args.ref)
    print(f"\nRéférence Sheet9 : {len(ref)} points\n")

    # --- CSV + statistiques ----------------------------------------------
    out_csv = "tacot_pyrolysis_gas_validation.csv"
    errs = {q[0]: [] for q in QUANTITIES}
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        head = ["P_atm", "T_K"]
        for name, _, _, unit, _ in QUANTITIES:
            head += [f"{name}_mpp[{unit}]", f"{name}_ref[{unit}]", f"{name}_err_pct"]
        w.writerow(head)
        for P in REF_PRESSURES_PA:
            for T in sorted(results[P]):
                key = (round(P, 2), T)
                if key not in ref:
                    continue
                row = [f"{P/ONEATM:g}", f"{T:g}"]
                for name, ci, fac, _, refcol in QUANTITIES:
                    m = results[P][T][ci] * fac
                    r = ref[key][refcol - 2]
                    floor = {"h": 50.0}.get(name, 1e-6)
                    e = rel_err(m, r, floor)
                    errs[name].append(e)
                    row += [f"{m:.6e}", f"{r:.6e}", f"{e:.4f}"]
                w.writerow(row)
    print(f"Comparaison point par point : {out_csv}")
    print(f"Points comparés             : {len(errs['h'])}\n")

    print("=" * 84)
    print("ÉCART Mutation++ / référence TACOT 3.0")
    print("=" * 84)
    for name, _, _, unit, _ in QUANTITIES:
        print(summarize(f"{name} [{unit}]", errs[name]))
    print("\n  (h : plancher de 50 kJ/kg au dénominateur — h_g change de signe")
    print("   vers 1400 K, l'écart relatif y est intrinsèquement mal défini)")
    print("\n  gamma : la référence sort de CEA, dont 'GAMMAs' est l'exposant")
    print("   isentropique -(dlnP/dlnV)_s, alors que gam_eq de Mutation++ est")
    print("   Cp_eq/Cv_eq. Les deux ne coïncident pas pour un mélange réactif ;")
    print("   l'écart résiduel est une différence de DÉFINITION, pas d'erreur.")
    print("  mu    : modèles de transport différents (règles de mélange).")

    # --- Table lisible à 1 atm -------------------------------------------
    print("\n" + "=" * 84)
    print("ENTHALPIE DU GAZ DE PYROLYSE À 1 atm")
    print("=" * 84)
    print(f"  {'T [K]':>7} | {'h M++ [kJ/kg]':>14} {'h réf [kJ/kg]':>14} "
          f"{'écart':>9} | {'M [kg/kmol]':>12}")
    for T in (200, 300, 500, 700, 1000, 1500, 2000, 2500, 3000, 3500, 3975):
        key = (round(101325.0, 2), T)
        if T not in results[101325.0] or key not in ref:
            continue
        m = results[101325.0][T][3] * 1e-3
        r = ref[key][3]
        print(f"  {T:>7d} | {m:>14.1f} {r:>14.1f} "
              f"{rel_err(m, r, 50.0):>8.2f} % | "
              f"{results[101325.0][T][1]*1e3:>12.3f}")

    plot_all(results, ref, "tacot_pyrolysis_gas_validation.png")

    print("\n" + "=" * 84)
    print("MISE EN GARDE — équilibre vs gaz réel")
    print("=" * 84)
    print("  h_g calculé ici suppose le gaz À L'ÉQUILIBRE. Or le gaz qui sort")
    print("  réellement de la résine vers 600-900 K est métastable : sa")
    print("  composition est celle de Sykes (H2, H2O, CH4, C6H5OH...), pas la")
    print("  composition d'équilibre. À basse température l'équilibre favorise")
    print("  CH4 + H2O, plus bas en enthalpie que le mélange réel.")
    print("\n  Conséquence : h_g(équilibre) sous-estime l'enthalpie du gaz réel")
    print("  à basse T. Si votre code applique en plus sa propre cinétique de")
    print("  craquage dans le char, il compte deux fois la même énergie.")


if __name__ == "__main__":
    main()
