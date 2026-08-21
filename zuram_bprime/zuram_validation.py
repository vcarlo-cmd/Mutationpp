#!/usr/bin/env python
"""
Valide la table B' du Zuram calculée par Mutation++ contre la table de
référence AblaNTIS `Bprime_carbonPhenolInAir_AblaNTIS.txt`.

Le calcul est refait sur la grille EXACTE de la référence, puis les deux
tables sont comparées point par point sur B'c et h_w.

Référence
---------
Fichier texte à 5 colonnes, sans en-tête :

    P[Pa]   B'g   B'c   T[K]   h_w[J/kg]

Grille : 10 pressions x 24 valeurs de B'g x 65 températures = 15600 points
    P    : 100, 1000, 2500, 5000, 7500, 1e4, 2.5e4, 5e4, 7.5e4, 1e5 Pa
    B'g  : 0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4, 0.5, 0.6,
           0.7, 0.8, 1.0, 1.2, 1.5, 1.9, 2.4, 3.0, 4.0, 5.5, 7.5, 10.0
    T    : 300 -> 3500 K, pas 50 K

Mise en données associée : carbonPhenolInAir_AblaNTIS.xml, reproduite dans
data/mixtures/zuram-air.xml (20 espèces gazeuses + C(gr), air 0.21/0.79,
gaz de pyrolyse C:0.171 H:0.722 O:0.107, char C:1.0).

Précision de la comparaison
---------------------------
La référence n'est stockée qu'avec 6 chiffres significatifs, et `bprime`
imprime lui aussi 6 chiffres : la résolution de la comparaison est donc
d'environ 1e-4 % (1 ppm) en relatif. Un écart mesuré nul signifie « identique
à tous les chiffres publiés par la référence », pas « identique au bit près ».
Les écarts sont rapportés en notation scientifique pour rester lisibles dans
ce régime.

Classement des points
---------------------
ACTIF   : B'c de référence non nul et sous le plafond -> comparable
ZERO    : B'c de référence exactement nul (le char n'est pas consommé)
SATURE  : B'c = max(100*B'g, 200), le plafond de char ajouté par
          `Thermodynamics::surfaceMassBalance`. Le char y est entièrement
          sublimé : la valeur n'est pas une solution physique mais la
          signature de B'c -> l'infini. Le plafond est identique dans les
          deux codes, l'accord y est donc trivial ; ces points sont isolés
          pour ne pas flatter les statistiques.

Usage :
    python zuram_validation.py [--mixture zuram-air] [--ref chemin.txt]

Prérequis :
    - binaire `bprime` dans le PATH ou dans build/src/apps/
    - numpy, matplotlib
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
# Paramètres
# ---------------------------------------------------------------------------
ONEATM = 101325.0  # Pa

T_RANGE = "300:50:3500"

BL_COMP   = "BLedge"
PYRO_COMP = "VKIZuramPyroGas"
CHAR_COMP = "Char"
CHAR_ELEM = "C"

DEFAULT_REF = os.path.expanduser(
    "~/.claude/uploads/e9fe3e8e-3fd8-58cd-a9e2-59c127c6614f/"
    "578096ff-Bprime_carbonPhenolInAir_AblaNTIS.txt")

# Pressions retenues pour les figures (les 4 extrêmes de la grille)
PLOT_PRESSURES_PA = [100.0, 1000.0, 10000.0, 100000.0]
PLOT_BG = [0.0, 0.1, 1.0, 10.0]


def char_cap(bg):
    """Plafond numérique de B'c : max(100*B'g, 200) — cf. surfaceMassBalance."""
    return max(100.0 * bg, 200.0)


def is_saturated(bc, bg):
    return bc >= char_cap(bg) * (1.0 - 1e-9)


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
    """Exécute bprime et retourne {T_arrondi: (B'c, hw_J_par_kg)}."""
    cmd = [path, "-T", T_RANGE, "-P", repr(pressure_pa), "-b", repr(bg),
           "-m", mixture, "-bl", BL_COMP, "-py", PYRO_COMP,
           "-char", CHAR_COMP, "-char-elem", CHAR_ELEM]
    res = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if res.returncode != 0:
        print(f"\nERREUR bprime (P={pressure_pa:g} Pa, B'g={bg:g}) :\n{res.stderr}")
        sys.exit(1)
    out = {}
    for line in res.stdout.strip().splitlines()[1:]:
        parts = line.split()
        try:
            T, bc, hw = float(parts[0]), float(parts[1]), float(parts[2])
        except (ValueError, IndexError):
            continue
        out[round(T)] = (bc, hw * 1.0e6)   # bprime écrit hw en MJ/kg
    return out


# ---------------------------------------------------------------------------
# Lecture de la table de référence
# ---------------------------------------------------------------------------

def load_reference(path):
    """
    Lit la table AblaNTIS (5 colonnes : P, B'g, B'c, T, h_w).
    Retourne (ref, pressions triées, B'g triés) avec
    ref = {(P_Pa, B'g, T): (B'c, h_w_J_par_kg)}.
    """
    ref = {}
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            parts = line.split()
            if len(parts) != 5:
                continue
            try:
                P, bg, bc, T, hw = (float(v) for v in parts)
            except ValueError:
                print(f"  ligne {lineno} ignorée : {line.strip()}")
                continue
            ref[(round(P, 6), round(bg, 6), round(T))] = (bc, hw)
    pressures = sorted({k[0] for k in ref})
    bgs = sorted({k[1] for k in ref})
    return ref, pressures, bgs


# ---------------------------------------------------------------------------
# Statistiques d'écart
# ---------------------------------------------------------------------------

def relative_error(computed, reference, floor):
    """Écart relatif [%], plancher sur |reference| pour éviter la division par ~0."""
    denom = max(abs(reference), floor)
    return (computed - reference) / denom * 100.0


def summarize(label, errs):
    """Résumé d'une série d'écarts relatifs, en notation scientifique."""
    a = np.abs(np.asarray(errs, dtype=float))
    if a.size == 0:
        return f"{label:<26s} (aucun point)"
    exact = 100.0 * np.mean(a == 0.0)
    return (f"{label:<26s} moy {a.mean():9.2e} %   "
            f"p95 {np.percentile(a, 95):9.2e} %   max {a.max():9.2e} %   "
            f"identiques {exact:6.2f} %   n={a.size}")


# ---------------------------------------------------------------------------
# Graphiques
# ---------------------------------------------------------------------------

def plot_overlay(results, ref, mixture, out_png):
    """Superpose calcul et référence (B'c et h_w) pour quatre B'g représentatifs."""
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    fig.suptitle(
        f"Table B' Zuram — Mutation++ ({mixture}) vs référence AblaNTIS\n"
        "traits épais clairs : référence   —   pointillés : Mutation++",
        fontsize=14)
    colors = plt.get_cmap("plasma", len(PLOT_PRESSURES_PA) + 1)

    for col, bg in enumerate(PLOT_BG):
        ax_bc, ax_hw = axes[0, col], axes[1, col]
        for i, P in enumerate(PLOT_PRESSURES_PA):
            Ts = sorted(T for (p, b, T) in ref if p == P and b == bg)
            if not Ts:
                continue
            lbl = f"{P:g} Pa"
            ax_bc.plot(Ts, [ref[(P, bg, T)][0] for T in Ts],
                       color=colors(i), lw=3.5, alpha=0.35, label=f"réf {lbl}")
            ax_hw.plot(Ts, [ref[(P, bg, T)][1] * 1e-6 for T in Ts],
                       color=colors(i), lw=3.5, alpha=0.35, label=f"réf {lbl}")
            Tc = [T for T in Ts if T in results[(P, bg)]]
            ax_bc.plot(Tc, [results[(P, bg)][T][0] for T in Tc],
                       color=colors(i), lw=1.3, ls="--", label=f"M++ {lbl}")
            ax_hw.plot(Tc, [results[(P, bg)][T][1] * 1e-6 for T in Tc],
                       color=colors(i), lw=1.3, ls="--", label=f"M++ {lbl}")

        ax_bc.axhline(char_cap(bg), color="0.4", lw=0.8, ls=":")
        ax_bc.set_yscale("log")
        ax_bc.set_title(rf"$B'_g = {bg:g}$")
        ax_bc.set_xlabel(r"$T_w$ [K]")
        ax_bc.set_ylabel(r"$B'_c$")
        ax_bc.grid(True, which="both", ls="--", alpha=0.35)
        ax_bc.legend(fontsize=6, ncol=2, loc="upper left")

        ax_hw.set_xlabel(r"$T_w$ [K]")
        ax_hw.set_ylabel(r"$h_w$ [MJ/kg]")
        ax_hw.grid(True, ls="--", alpha=0.35)
        ax_hw.legend(fontsize=6, ncol=2, loc="upper left")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"  Figure sauvegardée : {out_png}")


def plot_error_map(results, ref, pressures, bgs, mixture, out_png):
    """Carte de l'écart relatif sur B'c, un panneau par pression tracée."""
    fig, axes = plt.subplots(1, len(PLOT_PRESSURES_PA), figsize=(19, 4.8))
    fig.suptitle(
        rf"Écart relatif sur $B'_c$ — Mutation++ ({mixture}) vs référence AblaNTIS"
        "   (gris : $B'_c$ de référence nul ou saturé sur le plafond de char)",
        fontsize=13)
    Ts = sorted({T for (_, _, T) in ref})

    for ax, P in zip(np.atleast_1d(axes), PLOT_PRESSURES_PA):
        M = np.full((len(bgs), len(Ts)), np.nan)
        for i, bg in enumerate(bgs):
            for j, T in enumerate(Ts):
                key = (P, bg, T)
                if key not in ref or T not in results.get((P, bg), {}):
                    continue
                bc_r = ref[key][0]
                if bc_r == 0.0 or is_saturated(bc_r, bg):
                    continue
                M[i, j] = relative_error(results[(P, bg)][T][0], bc_r, 1e-3)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("0.85")
        # échelle symétrique auto ; plancher à la résolution de la référence
        # (6 chiffres significatifs -> ~1e-4 %) pour ne pas amplifier du bruit
        vmax = max(np.nanmax(np.abs(M)) if np.any(np.isfinite(M)) else 0.0, 1e-4)
        im = ax.pcolormesh(Ts, bgs, np.ma.masked_invalid(M), cmap=cmap,
                           vmin=-vmax, vmax=vmax, shading="nearest")
        ax.set_title(f"{P:g} Pa")
        ax.set_xlabel(r"$T_w$ [K]")
        ax.set_yscale("symlog", linthresh=0.02)
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
    ap.add_argument("--mixture", default="zuram-air",
                    help="nom du mélange (défaut : zuram-air)")
    ap.add_argument("--ref", default=DEFAULT_REF,
                    help="table B' de référence AblaNTIS (fichier texte)")
    args = ap.parse_args()

    path = find_bprime()
    if path is None:
        print("Binaire 'bprime' introuvable. Compilez :\n"
              "  cmake -B build -DCMAKE_BUILD_TYPE=Release .\n"
              "  cmake --build build --target bprime")
        sys.exit(1)
    if not os.path.isfile(args.ref):
        print(f"Table de référence introuvable : {args.ref}")
        sys.exit(1)

    print(f"Binaire   : {path}")
    print(f"Mélange   : {args.mixture}")
    print(f"Référence : {args.ref}\n")

    ref, pressures, bgs = load_reference(args.ref)
    Ts_ref = sorted({T for (_, _, T) in ref})
    print(f"Table de référence : {len(ref)} points  "
          f"({len(pressures)} P x {len(bgs)} B'g x {len(Ts_ref)} T)")
    print(f"  P   [Pa] : {', '.join(f'{p:g}' for p in pressures)}")
    print(f"  B'g      : {', '.join(f'{b:g}' for b in bgs)}")
    print(f"  T   [K]  : {Ts_ref[0]:g} -> {Ts_ref[-1]:g}\n")

    # --- Calcul sur la grille de référence --------------------------------
    results = {}
    total = len(pressures) * len(bgs)
    n = 0
    for P in pressures:
        for bg in bgs:
            n += 1
            print(f"\r  calcul {n:3d}/{total}  (P={P:g} Pa, B'g={bg:g})        ",
                  end="", flush=True)
            results[(P, bg)] = run_bprime(path, args.mixture, P, bg)
    print("\n")

    # --- Comparaison point par point --------------------------------------
    out_csv = f"zuram_bprime_vs_ref_{args.mixture}.csv"
    err_bc, err_hw = [], []          # régime ACTIF
    err_hw_all = []                  # h_w hors saturation (ACTIF + ZERO)
    per_bg = {bg: [] for bg in bgs}
    per_p = {P: [] for P in pressures}
    n_zero = n_zero_ok = n_sat = n_sat_ok = 0

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_Pa", "Bg", "T_K", "regime",
                    "Bc_mpp", "Bc_ref", "Bc_err_pct",
                    "hw_mpp_Jkg", "hw_ref_Jkg", "hw_err_pct"])
        for P in pressures:
            for bg in bgs:
                for T in sorted(results[(P, bg)]):
                    key = (P, bg, T)
                    if key not in ref:
                        continue
                    bc_m, hw_m = results[(P, bg)][T]
                    bc_r, hw_r = ref[key]
                    e_bc = relative_error(bc_m, bc_r, 1e-3)
                    e_hw = relative_error(hw_m, hw_r, 1e5)

                    if is_saturated(bc_r, bg):
                        regime = "SATURE"
                        n_sat += 1
                        n_sat_ok += is_saturated(bc_m, bg)
                    elif bc_r == 0.0:
                        regime = "ZERO"
                        n_zero += 1
                        n_zero_ok += bc_m < 1e-6
                        err_hw_all.append(e_hw)
                    else:
                        regime = "ACTIF"
                        err_bc.append(e_bc)
                        err_hw.append(e_hw)
                        err_hw_all.append(e_hw)
                        per_bg[bg].append(e_bc)
                        per_p[P].append(e_bc)

                    w.writerow([f"{P:g}", f"{bg:g}", f"{T:g}", regime,
                                f"{bc_m:.8e}", f"{bc_r:.8e}", f"{e_bc:.6f}",
                                f"{hw_m:.8e}", f"{hw_r:.8e}", f"{e_hw:.6f}"])

    total_pts = len(err_bc) + n_zero + n_sat
    print(f"Comparaison point par point : {out_csv}")
    print(f"Points comparés             : {total_pts} / {len(ref)}")
    print(f"  ACTIF  (comparables)      : {len(err_bc)}")
    print(f"  ZERO   (B'c_ref = 0)      : {n_zero}")
    print(f"  SATURE (plafond de char)  : {n_sat}\n")

    # --- Statistiques -----------------------------------------------------
    print("=" * 96)
    print("RÉGIME ACTIF — écart relatif Mutation++ / référence AblaNTIS")
    print("=" * 96)
    print(summarize("B'c", err_bc))
    print(summarize("hw", err_hw))

    def within(e, t):
        return 100.0 * np.mean(np.abs(np.asarray(e, dtype=float)) < t)

    if err_bc:
        print(f"\nB'c : {within(err_bc, 0.01):.2f} % des points < 0.01 %   "
              f"{within(err_bc, 0.1):.2f} % < 0.1 %   "
              f"{within(err_bc, 1):.2f} % < 1 %")
        print(f"hw  : {within(err_hw, 0.01):.2f} % des points < 0.01 %   "
              f"{within(err_hw, 0.1):.2f} % < 0.1 %   "
              f"{within(err_hw, 1):.2f} % < 1 %")

    print("\n" + "=" * 96)
    print("RÉGIME ZERO — la référence impose B'c = 0")
    print("=" * 96)
    frac = 100.0 * n_zero_ok / n_zero if n_zero else float("nan")
    print(f"Mutation++ donne aussi B'c < 1e-6 sur {n_zero_ok}/{n_zero} "
          f"points ({frac:.2f} %)")
    print(summarize("hw (ZERO + ACTIF)", err_hw_all))

    print("\n" + "=" * 96)
    print("RÉGIME SATURÉ — plafond max(100*B'g, 200) de surfaceMassBalance")
    print("=" * 96)
    frac = 100.0 * n_sat_ok / n_sat if n_sat else float("nan")
    print(f"Mutation++ sature sur le même plafond sur {n_sat_ok}/{n_sat} "
          f"points ({frac:.2f} %)")

    print("\n" + "=" * 96)
    print("ÉCART SUR B'c PAR PRESSION (régime actif)")
    print("=" * 96)
    for P in pressures:
        print(summarize(f"P = {P:g} Pa", per_p[P]))

    print("\n" + "=" * 96)
    print("ÉCART SUR B'c PAR B'g (régime actif)")
    print("=" * 96)
    for bg in bgs:
        print(summarize(f"B'g = {bg:g}", per_bg[bg]))

    # --- Verdict ----------------------------------------------------------
    worst = max((abs(e) for e in err_bc + err_hw), default=0.0)
    print("\n" + "=" * 96)
    print("VERDICT")
    print("=" * 96)
    print(f"Écart maximal sur B'c et h_w, tous régimes comparables confondus : "
          f"{worst:.3e} %")
    if worst == 0.0:
        print(f"Les deux tables sont identiques à tous les chiffres publiés par la\n"
              f"référence (6 chiffres significatifs) sur les {len(ref)} points de la\n"
              f"grille : B'c, h_w, les points à B'c nul et le plafond de saturation.")
    elif worst < 1e-3:
        print("Accord au niveau du ppm : l'écart est sous la résolution de la\n"
              "référence (6 chiffres significatifs).")
    else:
        print("Écart au-dessus de la résolution de la référence : vérifier la\n"
              "liste d'espèces et les compositions élémentaires du mélange.")

    # --- Graphiques -------------------------------------------------------
    print()
    plot_overlay(results, ref, args.mixture,
                 f"zuram_bprime_vs_ref_{args.mixture}.png")
    plot_error_map(results, ref, pressures, bgs, args.mixture,
                   f"zuram_bprime_err_{args.mixture}.png")
    print("\nTerminé.")


if __name__ == "__main__":
    main()
