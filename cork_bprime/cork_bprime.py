#!/usr/bin/env python
"""
Génère et trace la table B' d'un liège/phénolique (cork phenolic) dans l'air.

Matériau :
    - 80 % liège / 20 % résine phénolique (fractions MASSIQUES du solide)
    - rendement en char du COMPOSITE : 20 % (TGA du cork P50, argon,
      10 K/min — Sakraker et al., CEAS Space Journal 14:377-393, 2022)
    - résine phénolique (novolac C7H6O) : rendement en char 50 %
    - liège : 12.5 %, déduit des deux précédents

    Le liège PYROLYSE, contrairement aux fibres de carbone d'un TACOT ou
    d'un CPh70 : les deux constituants produisent du gaz, et le rapport
    liège/résine entre donc dans la composition du gaz de pyrolyse.

    Sur 100 g de vierge : 20 g de char + 80 g de gaz
    couplage stationnaire k = B'g/B'c = m_gaz/m_char = (1-y)/y = 4.0
    (contre 0.2727 pour le TACOT : quinze fois plus de soufflage pyrolytique)

    k se calcule sur les MASSES : l'identité (rho_v - rho_c)/rho_c suppose
    un volume constant, que ce matériau ne respecte pas (rho_c/rho_v = 0.62
    mesuré pour une masse résiduelle de 0.20).

Composition du gaz de pyrolyse (fermeture élémentaire liège + résine) :
    C:0.300, H:0.594, O:0.107 (fractions molaires élémentaires)
Char : carbone pur (liège carbonisé + résine carbonisée), C:1.0

Voir cork_pyrolysis_data.py pour le calcul de ces compositions et
mise_en_donnees_cork.md pour la mise en données complète.

Le mélange cork-air reprend la liste de 25 espèces de la table B' de
référence TACOT 3.0 (+ C(gr) pour la phase condensée).

Sorties :
    - une table CSV par valeur de B'g (25 isobares x plage de température)
    - un graphique isobares par valeur de B'g
    - la comparaison des B'g à 1 atm
    - B'c en fonction de B'g à température fixée, avec le point de
      fonctionnement stationnaire du liège/phénolique et celui du TACOT
    - la table du point de fonctionnement stationnaire

Usage :
    python cork_bprime.py

Prérequis :
    - binaire `bprime` dans le PATH ou dans build/src/apps/
    - data/mixtures/cork-air.xml
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
MIXTURE    = "cork-air"
BL_COMP    = "air"
PYRO_COMP  = "cork_pyro"
CHAR_COMP  = "cork_char"
CHAR_ELEM  = "C"

# Le liège/phénolique fonctionne à fort B'g : on pousse le balayage à 5.
BG_VALUES = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

ONEATM = 101325.0

PRESSURES_ATM      = np.logspace(-3, 3, 25)
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)

# Couplage matériau : B'g = k * B'c en ablation stationnaire
K_CORK  = 4.000       # m_gaz/m_char = 80 g gaz / 20 g char
K_TACOT = 0.2727      # (  280.0 -  220.0) /  220.0

# Masses volumiques (kg/m3) : n'entrent PAS dans le XML, seulement dans la
# récession s_dot = B'c * mdot_e / rho_c.
# Mesures P50 (Sakraker et al. 2022) : le char se rétracte, rho_c mesurée
# n'est pas 0.20*rho_v.
RHO_VIRGIN = 465.6
RHO_CHAR   = 289.1

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
        rf"Table B' — liège/phénolique dans l'air  ({bg_label(bg)},  "
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
        r"Influence de $B'_g$ sur la table B' — liège/phénolique dans l'air (P = 1 atm)",
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
    out_png = "cork_bprime_bg_comparison.png"
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

    La courbe tracée est celle du liège/phénolique. Les deux jeux de
    symboles montrent où la MEME courbe est lue selon le couplage matériau :
    k = 4.0 (liège/phénolique, étoiles) ou k = 0.2727 (couplage typique
    d'un carbone/phénolique, carrés). Le point de fonctionnement est
    l'intersection de la courbe avec la droite B'g = k*B'c.
    """
    fig, axes = plt.subplots(1, len(SWEEP_PRESSURES),
                             figsize=(8 * len(SWEEP_PRESSURES), 6))
    axes = np.atleast_1d(axes)
    fig.suptitle(
        r"$B'_c$ en fonction de $B'_g$ — liège/phénolique ; "
        r"le point de fonctionnement est l'intersection avec $B'_g = k\,B'_c$",
        fontsize=13)
    colors = plt.get_cmap("viridis", len(SWEEP_TEMPS) + 1)

    for ax, P in zip(axes, SWEEP_PRESSURES):
        for i, T in enumerate(SWEEP_TEMPS):
            bgs, bcs, _ = sweep[(P, T)]
            c = colors(i)
            ax.plot(bgs, bcs, color=c, lw=2, marker="o", ms=3,
                    label=rf"$T_w = {T}$ K")
            for k, mark, lab in ((K_CORK, "*", "liège/phénolique"),
                                 (K_TACOT, "s", "réf. carbone/phéno.")):
                x = solve_operating_point(bgs, bcs, k)
                ax.plot(k * x, x, mark, color=c, ms=13 if mark == "*" else 7,
                        mec="k", mew=0.8, zorder=5)
        # droites de fonctionnement B'g = k*B'c  ->  B'c = B'g/k
        bgline = np.array([1e-3, 10.0])
        ax.plot(bgline, bgline / K_CORK, "k--", lw=1.2, alpha=0.7,
                label=rf"liège/phéno. : $B'_g = {K_CORK:g}\,B'_c$  ($\star$)")
        ax.plot(bgline, bgline / K_TACOT, "k:", lw=1.2, alpha=0.7,
                label=rf"réf. carbone/phéno. : $B'_g = {K_TACOT}\,B'_c$"
                      + r"  ($\blacksquare$)")

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
    Table du point de fonctionnement stationnaire du liège/phénolique :
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
        rho_c = RHO_CHAR
        for P in SWEEP_PRESSURES:
            Ts = raw[(P, BG_SWEEP[0])][:, 0]
            for j, T in enumerate(Ts):
                bcs = np.array([raw[(P, bg)][j, 1] for bg in BG_SWEEP])
                hws = np.array([raw[(P, bg)][j, 2] for bg in BG_SWEEP])
                bc = solve_operating_point(np.array(BG_SWEEP), bcs, K_CORK)
                bg = min(K_CORK * bc, BG_SWEEP[-1])
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

        out_csv = f"cork_bprime_Bg{bg_filename(bg)}.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Bg", "P_atm"] + header)
            for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
                for row in data:
                    writer.writerow(
                        [f"{bg}", f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row])
        print(f"  Table CSV sauvegardée : {out_csv}")

        plot_isobar_table(all_data, PRESSURES_ATM, bg,
                          f"cork_bprime_Bg{bg_filename(bg)}.png")
        results_by_bg[bg] = (PRESSURES_ATM, all_data)
        print()

    plot_bg_comparison(results_by_bg)

    # --- B'c en fonction de B'g -------------------------------------------
    print(f"\n=== Balayage B'c(B'g) ===")
    sweep = bc_vs_bg(bprime_path)

    out_csv = "cork_bprime_Bc_vs_Bg.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "T_K", "Bg", "Bc", "hw_MJkg"])
        for (P, T), (bgs, bcs, hws) in sorted(sweep.items()):
            for bg, bc, hw in zip(bgs, bcs, hws):
                w.writerow([f"{P:g}", f"{T:g}", f"{bg:g}",
                            f"{bc:.6e}", f"{hw:.6e}"])
    print(f"\nTable B'c(B'g) sauvegardée : {out_csv}")

    plot_bc_vs_bg(sweep, "cork_bprime_Bc_vs_Bg.png")

    # --- Point de fonctionnement ------------------------------------------
    print()
    steady_state_table(bprime_path, "cork_bprime_steady_state.csv")

    # --- Récapitulatif -----------------------------------------------------
    print("\n" + "=" * 78)
    print("POINTS DE FONCTIONNEMENT — liège/phéno. (k = 4.0) vs TACOT (k = 0.2727)")
    print("=" * 78)
    print(f"{'P atm':>7} {'T K':>6} | {'cork Bc':>10} {'cork Bg':>9} | "
          f"{'k=0.273 Bc':>11} {'Bg':>9} | {'ecart Bc':>9}")
    for (P, T), (bgs, bcs, _) in sorted(sweep.items()):
        a = solve_operating_point(bgs, bcs, K_CORK)
        b = solve_operating_point(bgs, bcs, K_TACOT)
        rel = (a - b) / max(b, 1e-12) * 100.0
        print(f"{P:>7g} {T:>6g} | {a:>10.5f} {K_CORK*a:>9.4f} | "
              f"{b:>11.5f} {K_TACOT*b:>9.4f} | {rel:>8.2f} %")

    print("\nTerminé.")
