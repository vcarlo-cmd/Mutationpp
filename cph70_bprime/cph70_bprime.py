#!/usr/bin/env python
"""
Génère et trace la table B' du CPh70 dans l'air.

CPh70 : composite carbone/phénolique dense
    - 70 % fibres de carbone / 30 % résine phénolique (volume solide)
    - porosité 0.01
    - fibres et résine identiques à celles du TACOT
      (fibres ex-cellulose 1600 kg/m3, novolac/formaldéhyde 1200 kg/m3)

    rho_vierge = 1465.2 kg/m3      rho_char = 1287.0 kg/m3
    couplage stationnaire k = B'g/B'c = (rho_v - rho_c)/rho_c = 0.1385

Composition du gaz de pyrolyse (résine identique au TACOT) :
    C:0.206, H:0.679, O:0.115 (fractions molaires élémentaires)
Char : carbone pur (fibres + résine carbonisée), C:1.0

Le mélange cph70-air_25 reproduit la liste de 25 espèces de la table B' de
référence TACOT 3.0 (+ C(gr) pour la phase condensée).

Sorties :
    - une table CSV par valeur de B'g (25 isobares x plage de température)
    - un graphique isobares par valeur de B'g
    - la comparaison des B'g à 1 atm
    - B'c en fonction de B'g à température fixée, avec le point de
      fonctionnement stationnaire du CPh70 et celui du TACOT
    - la table du point de fonctionnement stationnaire

Usage :
    python cph70_bprime.py

Prérequis :
    - binaire `bprime` dans le PATH ou dans build/src/apps/
    - data/mixtures/cph70-air_25.xml
    - numpy, matplotlib
"""

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
# Paramètres du calcul
# ---------------------------------------------------------------------------
BPRIME_CMD = "bprime"
T_RANGE    = "300:25:5000"
MIXTURE    = "cph70-air_25"
BL_COMP    = "air"
PYRO_COMP  = "cph70_pyro"
CHAR_COMP  = "cph70_char"
CHAR_ELEM  = "C"

BG_VALUES = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0]

ONEATM = 101325.0

PRESSURES_ATM      = np.logspace(-3, 3, 25)
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)

# Couplage matériau : B'g = k * B'c en ablation stationnaire
K_CPH70  = 0.1385      # (1465.2 - 1287.0) / 1287.0
K_TACOT = 0.2727      # (  280.0 -  220.0) /  220.0

# Balayage fin en B'g pour la courbe B'c(B'g)
BG_SWEEP = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4,
            0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4,
            3.0, 4.0, 5.5, 7.5, 10.0]
SWEEP_T_RANGE   = "300:25:4000"
SWEEP_PRESSURES = [0.01, 1.0]          # atm
SWEEP_TEMPS     = [2000, 3000, 3400, 3700]   # K


# ---------------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------------

def find_bprime():
    """Localise le binaire bprime (PATH ou répertoires build courants)."""
    cmd = shutil.which(BPRIME_CMD)
    if cmd:
        return cmd
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, "../build/src/apps/bprime"),
              os.path.join(here, "../../build/src/apps/bprime"),
              "build/src/apps/bprime",
              "../build/src/apps/bprime"):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def make_env():
    """Environnement d'exécution avec MPP_DATA_DIRECTORY défini."""
    env = os.environ.copy()
    if "MPP_DATA_DIRECTORY" not in env:
        here = os.path.dirname(os.path.abspath(__file__))
        env["MPP_DATA_DIRECTORY"] = os.path.abspath(
            os.path.join(here, "../data"))
    return env


def run_bprime(bprime_path, pressure_pa, bg, t_range=T_RANGE):
    """Exécute bprime pour une pression (Pa) et un B'g donnés."""
    cmd = [
        bprime_path,
        "-T", t_range,
        "-P", str(pressure_pa),
        "-b", str(bg),
        "-m", MIXTURE,
        "-bl", BL_COMP,
        "-py", PYRO_COMP,
        "-char", CHAR_COMP,
        "-char-elem", CHAR_ELEM,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if result.returncode != 0:
        print(f"\nERREUR bprime à P={pressure_pa:.2f} Pa, B'g={bg} :")
        print(result.stderr)
        sys.exit(1)
    return result.stdout


def parse_output(output):
    """Parse la sortie de bprime. Retourne (header, data ndarray)."""
    lines = [l.strip() for l in output.strip().splitlines() if l.strip()]
    header = [h.strip('"') for h in lines[0].split()]
    data = []
    for line in lines[1:]:
        try:
            data.append([float(v) for v in line.split()])
        except ValueError:
            continue
    return header, np.array(data)


def bg_label(bg):
    return rf"$B'_g = {bg}$"


def bg_filename(bg):
    return str(bg).replace(".", "p")


# ---------------------------------------------------------------------------
# Table isobares (un graphique par B'g)
# ---------------------------------------------------------------------------

def plot_isobar_table(all_data, pressures_atm, bg, out_png):
    data_map = {P: d for P, d in zip(pressures_atm, all_data)}
    colors = plt.get_cmap("plasma", len(PLOT_PRESSURES_ATM) + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        rf"Table B' — CPh70 dans l'air  ({bg_label(bg)},  "
        r"$P \in [10^{-3},\,10^3]$ atm)",
        fontsize=13,
    )

    for idx, P_atm in enumerate(PLOT_PRESSURES_ATM):
        closest = min(pressures_atm,
                      key=lambda p: abs(np.log10(p) - np.log10(P_atm)))
        _, data = data_map[closest]
        Tw, Bc, hw = data[:, 0], data[:, 1], data[:, 2]
        exp = int(round(np.log10(P_atm)))
        lbl = rf"$10^{{{exp}}}$ atm" if exp != 0 else "1 atm"
        ax1.plot(Tw, Bc, color=colors(idx), lw=2, label=lbl)
        ax2.plot(Tw, hw, color=colors(idx), lw=2, label=lbl)

    ax1.set_yscale("log")
    ax1.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(rf"Taux d'ablation char $B'_c$  ({bg_label(bg)})")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    ax2.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(rf"Enthalpie de paroi $h_w$  ({bg_label(bg)})")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"  Figure isobares sauvegardée : {out_png}")
    plt.close()


def plot_bg_comparison(results_by_bg):
    colors = plt.get_cmap("tab10", len(BG_VALUES))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        r"Influence de $B'_g$ sur la table B' — CPh70 dans l'air (P = 1 atm)",
        fontsize=13,
    )
    for idx, bg in enumerate(BG_VALUES):
        pressures_atm, all_data = results_by_bg[bg]
        _, data = min(zip(pressures_atm, all_data),
                      key=lambda x: abs(np.log10(x[0])))
        _, data = data
        Tw, Bc, hw = data[:, 0], data[:, 1], data[:, 2]
        ax1.plot(Tw, Bc, color=colors(idx), lw=2, label=bg_label(bg))
        ax2.plot(Tw, hw, color=colors(idx), lw=2, label=bg_label(bg))

    ax1.set_yscale("log")
    ax1.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(r"Taux d'ablation $B'_c$ vs $B'_g$")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=10, loc="upper left")

    ax2.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"Enthalpie de paroi $h_w$ vs $B'_g$")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=10, loc="upper left")

    plt.tight_layout()
    out_png = "cph70_bprime_bg_comparison.png"
    plt.savefig(out_png, dpi=150)
    print(f"\nFigure de comparaison B'g sauvegardée : {out_png}")
    plt.close()


# ---------------------------------------------------------------------------
# B'c en fonction de B'g, à température fixée
# ---------------------------------------------------------------------------

def solve_operating_point(bgs, bc_of_bg, k, iters=300):
    """Point fixe B'c = f(B'g = k*B'c), amorti (la table est très raide)."""
    x = bc_of_bg[0]
    for _ in range(iters):
        xn = np.interp(min(k * x, bgs[-1]), bgs, bc_of_bg)
        if not np.isfinite(xn):
            break
        if abs(xn - x) < 1e-12 * max(abs(x), 1e-12):
            x = xn
            break
        x = 0.5 * (x + xn)
    return x


def bc_vs_bg(bprime_path):
    """
    Calcule B'c(B'g) à températures fixées, pour deux pressions.
    Retourne {(P, T): (bgs, bc, hw)}.
    """
    out = {}
    raw = {}
    for P in SWEEP_PRESSURES:
        for bg in BG_SWEEP:
            print(f"  balayage  P = {P:g} atm, B'g = {bg:g} ...",
                  end=" ", flush=True)
            o = run_bprime(bprime_path, P * ONEATM, bg, SWEEP_T_RANGE)
            _, d = parse_output(o)
            raw[(P, bg)] = d
            print(f"{len(d)} points")
    for P in SWEEP_PRESSURES:
        Ts = raw[(P, BG_SWEEP[0])][:, 0]
        for T in SWEEP_TEMPS:
            j = int(np.argmin(np.abs(Ts - T)))
            bcs = np.array([raw[(P, bg)][j, 1] for bg in BG_SWEEP])
            hws = np.array([raw[(P, bg)][j, 2] for bg in BG_SWEEP])
            out[(P, T)] = (np.array(BG_SWEEP), bcs, hws)
    return out


def plot_bc_vs_bg(sweep, out_png):
    """
    B'c en fonction de B'g à T fixée, avec les points de fonctionnement.

    La courbe B'c(B'g) est la MEME pour CPh70 et TACOT (même table) ; seuls
    les points de fonctionnement diffèrent, car ils sont l'intersection de
    cette courbe avec la droite B'g = k*B'c propre à chaque matériau.
    """
    fig, axes = plt.subplots(1, len(SWEEP_PRESSURES),
                             figsize=(8 * len(SWEEP_PRESSURES), 6))
    axes = np.atleast_1d(axes)
    fig.suptitle(
        r"$B'_c$ en fonction de $B'_g$ — la courbe est commune CPh70 / TACOT, "
        r"seuls les points de fonctionnement $B'_g = k\,B'_c$ diffèrent",
        fontsize=13)
    colors = plt.get_cmap("viridis", len(SWEEP_TEMPS) + 1)

    for ax, P in zip(axes, SWEEP_PRESSURES):
        for i, T in enumerate(SWEEP_TEMPS):
            bgs, bcs, _ = sweep[(P, T)]
            c = colors(i)
            ax.plot(bgs, bcs, color=c, lw=2, marker="o", ms=3,
                    label=rf"$T_w = {T}$ K")
            for k, mark, lab in ((K_CPH70, "*", "CPh70"),
                                 (K_TACOT, "s", "TACOT")):
                x = solve_operating_point(bgs, bcs, k)
                ax.plot(k * x, x, mark, color=c, ms=13 if mark == "*" else 7,
                        mec="k", mew=0.8, zorder=5)
        # droites de fonctionnement B'g = k*B'c  ->  B'c = B'g/k
        bgline = np.array([1e-3, 10.0])
        ax.plot(bgline, bgline / K_CPH70, "k--", lw=1.2, alpha=0.7,
                label=rf"CPh70 : $B'_g = {K_CPH70}\,B'_c$  ($\star$)")
        ax.plot(bgline, bgline / K_TACOT, "k:", lw=1.2, alpha=0.7,
                label=rf"TACOT : $B'_g = {K_TACOT}\,B'_c$  ($\blacksquare$)")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1.5e-2, 11)
        ax.set_ylim(5e-2, 3e2)
        exp = int(round(np.log10(P)))
        ax.set_title(rf"$10^{{{exp}}}$ atm" if exp else "1 atm")
        ax.set_xlabel(r"$B'_g$")
        ax.set_ylabel(r"$B'_c$")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=8, loc="upper left")

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Figure B'c(B'g) sauvegardée : {out_png}")
    plt.close()


# ---------------------------------------------------------------------------
# Point de fonctionnement stationnaire
# ---------------------------------------------------------------------------

def steady_state_table(bprime_path, out_csv):
    """
    Table du point de fonctionnement stationnaire du CPh70 :
    pour chaque (P, T), résout B'c = table(T, P, B'g = k*B'c).
    """
    raw = {}
    for P in SWEEP_PRESSURES:
        for bg in BG_SWEEP:
            raw[(P, bg)] = parse_output(
                run_bprime(bprime_path, P * ONEATM, bg, SWEEP_T_RANGE))[1]

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "T_K", "Bc_ss", "Bg_ss", "hw_ss_MJkg",
                    "Bc_Bg0", "recession_over_mdote_m3_per_kg"])
        rho_c = 1287.0
        for P in SWEEP_PRESSURES:
            Ts = raw[(P, BG_SWEEP[0])][:, 0]
            for j, T in enumerate(Ts):
                bcs = np.array([raw[(P, bg)][j, 1] for bg in BG_SWEEP])
                hws = np.array([raw[(P, bg)][j, 2] for bg in BG_SWEEP])
                bc = solve_operating_point(np.array(BG_SWEEP), bcs, K_CPH70)
                bg = min(K_CPH70 * bc, BG_SWEEP[-1])
                hw = np.interp(bg, BG_SWEEP, hws)
                w.writerow([f"{P:g}", f"{T:g}", f"{bc:.6e}", f"{bg:.6e}",
                            f"{hw:.6e}", f"{bcs[0]:.6e}", f"{bc/rho_c:.6e}"])
    print(f"Table du point de fonctionnement : {out_csv}")


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

if __name__ == "__main__":
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

    # --- Tables principales ------------------------------------------------
    results_by_bg = {}
    for bg in BG_VALUES:
        print(f"=== B'g = {bg} ===")
        all_data = []
        for P_atm in PRESSURES_ATM:
            P_pa = P_atm * ONEATM
            print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
                  end=" ", flush=True)
            header, data = parse_output(run_bprime(bprime_path, P_pa, bg))
            all_data.append((header, data))
            print(f"{len(data)} points")

        out_csv = f"cph70_bprime_Bg{bg_filename(bg)}.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Bg", "P_atm"] + header)
            for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
                for row in data:
                    writer.writerow(
                        [f"{bg}", f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row])
        print(f"  Table CSV sauvegardée : {out_csv}")

        plot_isobar_table(all_data, PRESSURES_ATM, bg,
                          f"cph70_bprime_Bg{bg_filename(bg)}.png")
        results_by_bg[bg] = (PRESSURES_ATM, all_data)
        print()

    plot_bg_comparison(results_by_bg)

    # --- B'c en fonction de B'g -------------------------------------------
    print(f"\n=== Balayage B'c(B'g) ===")
    sweep = bc_vs_bg(bprime_path)

    out_csv = "cph70_bprime_Bc_vs_Bg.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "T_K", "Bg", "Bc", "hw_MJkg"])
        for (P, T), (bgs, bcs, hws) in sorted(sweep.items()):
            for bg, bc, hw in zip(bgs, bcs, hws):
                w.writerow([f"{P:g}", f"{T:g}", f"{bg:g}",
                            f"{bc:.6e}", f"{hw:.6e}"])
    print(f"\nTable B'c(B'g) sauvegardée : {out_csv}")

    plot_bc_vs_bg(sweep, "cph70_bprime_Bc_vs_Bg.png")

    # --- Point de fonctionnement ------------------------------------------
    print()
    steady_state_table(bprime_path, "cph70_bprime_steady_state.csv")

    # --- Récapitulatif -----------------------------------------------------
    print("\n" + "=" * 78)
    print("POINTS DE FONCTIONNEMENT — CPh70 (k = 0.1385) vs TACOT (k = 0.2727)")
    print("=" * 78)
    print(f"{'P atm':>7} {'T K':>6} | {'CPh70 Bc':>10} {'CPh70 Bg':>9} | "
          f"{'TACOT Bc':>10} {'TACOT Bg':>9} | {'ecart Bc':>9}")
    for (P, T), (bgs, bcs, _) in sorted(sweep.items()):
        a = solve_operating_point(bgs, bcs, K_CPH70)
        b = solve_operating_point(bgs, bcs, K_TACOT)
        rel = (a - b) / max(b, 1e-12) * 100.0
        print(f"{P:>7g} {T:>6g} | {a:>10.5f} {K_CPH70*a:>9.4f} | "
              f"{b:>10.5f} {K_TACOT*b:>9.4f} | {rel:>8.2f} %")

    print("\nTerminé.")
