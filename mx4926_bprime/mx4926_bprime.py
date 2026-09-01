#!/usr/bin/env python
"""
Génère et trace la table B' du MX-4926 dans l'air.

MX-4926 : composite carbone/phénolique chargé, matériau de tuyère SRM
    - renfort  52.2 % masse : fibres de carbone ex-rayon, satin de 5 ou 8
    - matrice  34.2 % masse : résol phénol-formaldéhyde SC-1008
    - charge   13.6 % masse : noir de carbone
      (centroïde exact du domaine admissible — cf. composition_mx4926.md)

    rho_vierge = 1451 kg/m3      rho_char = 1228 kg/m3   (porosité 0.02)
    couplage stationnaire k = B'g/B'c = (1 - Y_comp)/Y_comp = 0.1817

Composition du gaz de pyrolyse : HÉRITÉE du SC-1008, la matrice étant la même
    C:0.2526, H:0.6407, O:0.1068 (fractions molaires élémentaires)
Char : carbone pur — fibres ex-rayon, noir de carbone et résol carbonisé
    C:1.0

La table B'c(T, P, B'g) est donc rigoureusement celle du SC-1008 (identité
vérifiée bit à bit par composition_mx4926_verification.py §5). Ce qui est
propre au MX-4926 est le POINT DE FONCTIONNEMENT : l'intersection de la
courbe commune B'c(B'g) avec la droite B'g = k·B'c, k = 0.1817.

Sorties :
    - une table CSV par valeur de B'g (25 isobares x plage de température)
    - un graphique isobares par valeur de B'g
    - la comparaison des B'g à 1 atm
    - B'c en fonction de B'g à température fixée, avec les points de
      fonctionnement du MX-4926, du PICA et du TACOT
    - la table du point de fonctionnement stationnaire

Usage :
    python mx4926_bprime.py

Prérequis :
    - binaire `bprime` dans le PATH ou dans build/src/apps/
    - data/mixtures/mx4926-air.xml
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
MIXTURE    = "mx4926-air"
BL_COMP    = "air"
PYRO_COMP  = "mx4926_pyro"
CHAR_COMP  = "mx4926_char"
CHAR_ELEM  = "C"

BG_VALUES = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0]

ONEATM = 101325.0

PRESSURES_ATM      = np.logspace(-3, 3, 25)
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)

# Couplage matériau : B'g = k * B'c en ablation stationnaire.
#   k = (rho_v - rho_c)/rho_c = (1 - Y_comp)/Y_comp — un rapport de MASSES,
#   indépendant des densités et de la porosité.
K_MX4926 = 0.1817     # Y_comp = 1 - 0.45*0.341667 = 0.84625
K_PICA   = 0.2346     # PICA / SC-1008 : 0.19/0.81  (resine_sc1008.md §5.1)
K_TACOT  = 0.2727     # (280.0 - 220.0)/220.0

# Enveloppe de k sur le domaine de composition admissible
K_MX4926_MIN = 0.1621   # w_resine = 0.31
K_MX4926_MAX = 0.1998   # w_resine = 0.37

RHO_CHAR = 1228.0     # kg/m3, porosité 0.02

# Balayage fin en B'g pour la courbe B'c(B'g)
BG_SWEEP = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4,
            0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4,
            3.0, 4.0, 5.5, 7.5, 10.0]
SWEEP_T_RANGE   = "300:25:4000"
SWEEP_PRESSURES = [0.01, 1.0]                # atm
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
        rf"Table B' — MX-4926 dans l'air  ({bg_label(bg)},  "
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
        r"Influence de $B'_g$ sur la table B' — MX-4926 dans l'air (P = 1 atm)",
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
    out_png = "mx4926_bprime_bg_comparison.png"
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

    La courbe B'c(B'g) est COMMUNE au MX-4926, au PICA et au TACOT — c'est la
    même résine phénolique et le même char de carbone pur. Seuls les points de
    fonctionnement diffèrent : ils sont l'intersection de cette courbe avec la
    droite B'g = k·B'c propre à chaque matériau, et k n'est fixé que par le
    rendement en char du COMPOSITE.
    """
    fig, axes = plt.subplots(1, len(SWEEP_PRESSURES),
                             figsize=(8 * len(SWEEP_PRESSURES), 6))
    axes = np.atleast_1d(axes)
    fig.suptitle(
        r"$B'_c$ en fonction de $B'_g$ — courbe commune MX-4926 / PICA / TACOT, "
        r"seuls les points de fonctionnement $B'_g = k\,B'_c$ diffèrent",
        fontsize=13)
    colors = plt.get_cmap("viridis", len(SWEEP_TEMPS) + 1)

    for ax, P in zip(axes, SWEEP_PRESSURES):
        for i, T in enumerate(SWEEP_TEMPS):
            bgs, bcs, _ = sweep[(P, T)]
            c = colors(i)
            ax.plot(bgs, bcs, color=c, lw=2, marker="o", ms=3,
                    label=rf"$T_w = {T}$ K")
            for k, mark, ms in ((K_MX4926, "*", 13), (K_PICA, "D", 6),
                                (K_TACOT, "s", 7)):
                x = solve_operating_point(bgs, bcs, k)
                ax.plot(k * x, x, mark, color=c, ms=ms, mec="k", mew=0.8,
                        zorder=5)
        # droites de fonctionnement B'g = k*B'c  ->  B'c = B'g/k
        bgline = np.array([1e-3, 10.0])
        ax.plot(bgline, bgline / K_MX4926, "k--", lw=1.2, alpha=0.7,
                label=rf"MX-4926 : $B'_g = {K_MX4926}\,B'_c$  ($\star$)")
        ax.plot(bgline, bgline / K_PICA, "k-.", lw=1.0, alpha=0.6,
                label=rf"PICA : $B'_g = {K_PICA}\,B'_c$  ($\diamond$)")
        ax.plot(bgline, bgline / K_TACOT, "k:", lw=1.2, alpha=0.7,
                label=rf"TACOT : $B'_g = {K_TACOT}\,B'_c$  ($\blacksquare$)")
        # enveloppe du domaine de composition du MX-4926
        ax.fill_betweenx([5e-2, 3e2],
                         K_MX4926_MIN * np.array([5e-2, 3e2]),
                         K_MX4926_MAX * np.array([5e-2, 3e2]),
                         color="0.5", alpha=0.15, zorder=0,
                         label=r"MX-4926 : enveloppe $k \in [0.162,\,0.200]$")

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
    Table du point de fonctionnement stationnaire du MX-4926 :
    pour chaque (P, T), résout B'c = table(T, P, B'g = k*B'c).

    Les colonnes Bc_kmin / Bc_kmax donnent le même point de fonctionnement aux
    deux bornes du domaine de composition admissible : c'est l'incertitude
    RÉELLE induite par les plages de composition annoncées.
    """
    raw = {}
    for P in SWEEP_PRESSURES:
        for bg in BG_SWEEP:
            raw[(P, bg)] = parse_output(
                run_bprime(bprime_path, P * ONEATM, bg, SWEEP_T_RANGE))[1]

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "T_K", "Bc_ss", "Bg_ss", "hw_ss_MJkg",
                    "Bc_Bg0", "Bc_kmin", "Bc_kmax",
                    "recession_over_mdote_m3_per_kg"])
        for P in SWEEP_PRESSURES:
            Ts = raw[(P, BG_SWEEP[0])][:, 0]
            for j, T in enumerate(Ts):
                bcs = np.array([raw[(P, bg)][j, 1] for bg in BG_SWEEP])
                hws = np.array([raw[(P, bg)][j, 2] for bg in BG_SWEEP])
                bc = solve_operating_point(np.array(BG_SWEEP), bcs, K_MX4926)
                lo = solve_operating_point(np.array(BG_SWEEP), bcs, K_MX4926_MIN)
                hi = solve_operating_point(np.array(BG_SWEEP), bcs, K_MX4926_MAX)
                bg = min(K_MX4926 * bc, BG_SWEEP[-1])
                hw = np.interp(bg, BG_SWEEP, hws)
                w.writerow([f"{P:g}", f"{T:g}", f"{bc:.6e}", f"{bg:.6e}",
                            f"{hw:.6e}", f"{bcs[0]:.6e}",
                            f"{lo:.6e}", f"{hi:.6e}",
                            f"{bc/RHO_CHAR:.6e}"])
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

        out_csv = f"mx4926_bprime_Bg{bg_filename(bg)}.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Bg", "P_atm"] + header)
            for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
                for row in data:
                    writer.writerow(
                        [f"{bg}", f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row])
        print(f"  Table CSV sauvegardée : {out_csv}")

        plot_isobar_table(all_data, PRESSURES_ATM, bg,
                          f"mx4926_bprime_Bg{bg_filename(bg)}.png")
        results_by_bg[bg] = (PRESSURES_ATM, all_data)
        print()

    plot_bg_comparison(results_by_bg)

    # Table B'c au format long, unites SI : Tw_K, P_bar, Bg, Bc
    # (nTw x nP x nBg lignes)
    out_csv = "mx4926_bprime_bc_table.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Tw_K", "P_bar", "Bg", "Bc"])
        for bg in BG_VALUES:
            pressures_atm, all_data_bg = results_by_bg[bg]
            for P_atm, (_, data) in zip(pressures_atm, all_data_bg):
                P_bar = P_atm * ONEATM / 1.0e5
                for row in data:
                    w.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}", f"{bg:g}",
                                f"{row[1]:.6e}"])
    print(f"Table B'c (nTw x nP x nBg) sauvegardée : {out_csv}")

    # Table h_w au format long, unites SI, a B'g = 0 : Tw_K, P_bar, hw_Jkg
    # (nTw x nP lignes)
    out_csv = "mx4926_bprime_hw_table.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Tw_K", "P_bar", "hw_Jkg"])
        pressures_atm, all_data_bg0 = results_by_bg[0.0]
        for P_atm, (_, data) in zip(pressures_atm, all_data_bg0):
            P_bar = P_atm * ONEATM / 1.0e5
            for row in data:
                w.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}",
                            f"{row[2] * 1.0e6:.6e}"])
    print(f"Table h_w (nTw x nP, B'g=0) sauvegardée : {out_csv}")

    # --- B'c en fonction de B'g -------------------------------------------
    print("\n=== Balayage B'c(B'g) ===")
    sweep = bc_vs_bg(bprime_path)

    out_csv = "mx4926_bprime_Bc_vs_Bg.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "T_K", "Bg", "Bc", "hw_MJkg"])
        for (P, T), (bgs, bcs, hws) in sorted(sweep.items()):
            for bg, bc, hw in zip(bgs, bcs, hws):
                w.writerow([f"{P:g}", f"{T:g}", f"{bg:g}",
                            f"{bc:.6e}", f"{hw:.6e}"])
    print(f"\nTable B'c(B'g) sauvegardée : {out_csv}")

    plot_bc_vs_bg(sweep, "mx4926_bprime_Bc_vs_Bg.png")

    # --- Point de fonctionnement ------------------------------------------
    print()
    steady_state_table(bprime_path, "mx4926_bprime_steady_state.csv")

    # --- Récapitulatif -----------------------------------------------------
    print("\n" + "=" * 78)
    print(f"POINTS DE FONCTIONNEMENT — MX-4926 (k = {K_MX4926}) vs "
          f"PICA (k = {K_PICA}) vs TACOT (k = {K_TACOT})")
    print("=" * 78)
    print(f"{'P atm':>7} {'T K':>6} | {'MX Bc':>9} {'MX Bg':>8} | "
          f"{'PICA Bc':>9} | {'TACOT Bc':>9} | {'MX/TACOT':>9}")
    for (P, T), (bgs, bcs, _) in sorted(sweep.items()):
        a = solve_operating_point(bgs, bcs, K_MX4926)
        p = solve_operating_point(bgs, bcs, K_PICA)
        b = solve_operating_point(bgs, bcs, K_TACOT)
        rel = (a - b) / max(b, 1e-12) * 100.0
        print(f"{P:>7g} {T:>6g} | {a:>9.5f} {K_MX4926*a:>8.4f} | "
              f"{p:>9.5f} | {b:>9.5f} | {rel:>8.2f} %")

    print("\n" + "=" * 78)
    print("ENVELOPPE DE COMPOSITION — B'c stationnaire aux bornes du domaine")
    print("=" * 78)
    print(f"{'P atm':>7} {'T K':>6} | {'k=0.1621':>10} {'k=0.1817':>10} "
          f"{'k=0.1998':>10} | {'étendue':>9}")
    for (P, T), (bgs, bcs, _) in sorted(sweep.items()):
        lo = solve_operating_point(bgs, bcs, K_MX4926_MIN)
        mid = solve_operating_point(bgs, bcs, K_MX4926)
        hi = solve_operating_point(bgs, bcs, K_MX4926_MAX)
        sp = (max(lo, hi) - min(lo, hi)) / max(mid, 1e-12) * 100.0
        print(f"{P:>7g} {T:>6g} | {lo:>10.5f} {mid:>10.5f} {hi:>10.5f} | "
              f"{sp:>8.2f} %")

    print("\nTerminé.")
