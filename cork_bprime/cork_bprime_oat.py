#!/usr/bin/env python
"""
Table B' du liège/phénolique (cork P50) sous flamme oxyacétylénique —
sortie de torche OAT (type ASTM E285) — et comparaison à la table air.

Seul le bord de couche limite change par rapport à cork_bprime.py : même
gaz de pyrolyse (cork_pyro), même char (C pur), même liste d'espèces
(data/mixtures/cork-oat.xml reprend celle de cork-air.xml).

Bord de couche limite = produits de C2H2 + r O2, sans entraînement d'air :
éléments C:H:O = 2:2:2r. Trois réglages de torche :

    r = 1.0  flamme neutre         -> C/O = 1 : pas d'oxygène « libre »
    r = 1.3  réglage OAT usuel     -> 0.6 mol O en excès pour 2 C
    r = 2.5  stoechiométrique      -> 3 mol O en excès pour 2 C

Plateau d'oxydation attendu (tout l'O en excès sur C forme du CO) :

    B'c = M_C * (2r - 2) / (M_C2H2 + r M_O2)
        r = 1.0 : 0         r = 1.3 : 0.106     r = 2.5 : 0.340
        (air    : 0.175)

Sorties :
    cork_oat_bprime_bc_table.csv     Tw_K, P_bar, O2_C2H2, Bg, Bc
    cork_oat_bprime_hw_table.csv     Tw_K, P_bar, O2_C2H2, Bg, hw_Jkg
    cork_oat_bprime_steady_state.csv point de fonctionnement k = 4 (+ air)
    cork_oat_bprime_vs_air.png       B'c et h_w à 1 atm, B'g = 0 et 5
    cork_oat_bprime_steady_state.png B'c stationnaire, air vs OAT

Usage :
    python cork_bprime_oat.py
"""

import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cork_bprime import (find_bprime, make_env, parse_output,
                         solve_operating_point, K_CORK, RHO_CHAR, ONEATM)
import subprocess

MIXTURE   = "cork-oat"
PYRO_COMP = "cork_pyro"
CHAR_COMP = "cork_char"
CHAR_ELEM = "C"

# (étiquette composition XML, rapport volumique O2/C2H2, libellé)
EDGES = [
    ("air",      None, "air"),
    ("oat_r1p0", 1.0,  "OAT O$_2$/C$_2$H$_2$ = 1.0 (neutre)"),
    ("oat_r1p3", 1.3,  "OAT O$_2$/C$_2$H$_2$ = 1.3"),
    ("oat_r2p5", 2.5,  "OAT O$_2$/C$_2$H$_2$ = 2.5 (stoech.)"),
]

T_RANGE       = "300:25:5000"
PRESSURES_ATM = np.logspace(-3, 3, 25)
BG_VALUES     = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

BG_SWEEP = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.25, 0.32, 0.4,
            0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.9, 2.4,
            3.0, 4.0, 5.5, 7.5, 10.0]
SWEEP_T_RANGE   = "300:25:4000"
SWEEP_PRESSURES = [0.01, 1.0]

M_C, M_H, M_O = 12.011, 1.008, 15.999


def plateau(r):
    """B'c de plateau d'oxydation (O en excès -> CO) pour C2H2 + r O2."""
    return M_C * max(2 * r - 2, 0) / (2 * M_C + 2 * M_H + 2 * r * M_O)


def run(bprime_path, bl, P_pa, bg, t_range=T_RANGE):
    cmd = [bprime_path, "-T", t_range, "-P", str(P_pa), "-b", str(bg),
           "-m", MIXTURE, "-bl", bl, "-py", PYRO_COMP,
           "-char", CHAR_COMP, "-char-elem", CHAR_ELEM]
    res = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if res.returncode != 0:
        print(f"ERREUR bprime ({bl}, P={P_pa} Pa, B'g={bg}) :\n{res.stderr}")
        sys.exit(1)
    return parse_output(res.stdout)[1]


def main():
    bprime_path = find_bprime()
    if bprime_path is None:
        sys.exit("Binaire 'bprime' introuvable (cf. cork_bprime.py).")

    # --- Tables B'c / h_w, format long ------------------------------------
    tables = {}   # (bl, P_atm, bg) -> data
    for bl, r, _ in EDGES:
        for bg in BG_VALUES:
            print(f"  {bl:9s} B'g = {bg:<4g}", flush=True)
            for P in PRESSURES_ATM:
                tables[(bl, P, bg)] = run(bprime_path, bl, P * ONEATM, bg)

    with open("cork_oat_bprime_bc_table.csv", "w", newline="") as fb, \
         open("cork_oat_bprime_hw_table.csv", "w", newline="") as fh:
        wb, wh = csv.writer(fb), csv.writer(fh)
        wb.writerow(["Tw_K", "P_bar", "O2_C2H2", "Bg", "Bc"])
        wh.writerow(["Tw_K", "P_bar", "O2_C2H2", "Bg", "hw_Jkg"])
        for bl, r, _ in EDGES:
            if r is None:
                continue   # la table air est déjà dans cork_bprime_*.csv
            for bg in BG_VALUES:
                for P in PRESSURES_ATM:
                    P_bar = P * ONEATM / 1.0e5
                    for row in tables[(bl, P, bg)]:
                        wb.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}", f"{r:g}",
                                     f"{bg:g}", f"{row[1]:.6e}"])
                        wh.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}", f"{r:g}",
                                     f"{bg:g}", f"{row[2] * 1e6:.6e}"])
    print("Tables : cork_oat_bprime_bc_table.csv, cork_oat_bprime_hw_table.csv")

    # --- Figure B'c / h_w à 1 atm -----------------------------------------
    P1 = PRESSURES_ATM[np.argmin(np.abs(np.log10(PRESSURES_ATM)))]
    colors = {"air": "k", "oat_r1p0": "tab:blue", "oat_r1p3": "tab:red",
              "oat_r2p5": "tab:orange"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Table B' — liège/phénolique P50 : air vs sortie de torche "
                 "oxyacétylénique (1 atm)", fontsize=13)
    for bl, r, lab in EDGES:
        for bg, ls in ((0.0, "-"), (5.0, "--")):
            d = tables[(bl, P1, bg)]
            ax1.plot(d[:, 0], np.maximum(d[:, 1], 1e-5), ls, color=colors[bl],
                     lw=2, label=f"{lab}, B'g = {bg:g}")
            ax2.plot(d[:, 0], d[:, 2], ls, color=colors[bl], lw=2)
        if r is not None and r > 1:
            ax1.axhline(plateau(r), color=colors[bl], lw=0.8, ls=":")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-4, 3e2)
    ax1.set_xlabel(r"$T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$")
    ax1.set_title(r"$B'_c$ (pointillés : plateau O en excès $\to$ CO)")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=7.5, loc="upper left")
    ax2.set_xlabel(r"$T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"$h_w$ (trait plein B'g = 0, tirets B'g = 5)")
    ax2.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig("cork_oat_bprime_vs_air.png", dpi=150)
    plt.close()

    # --- Point de fonctionnement stationnaire (k = 4) ---------------------
    ss = {}
    for bl, r, _ in EDGES:
        for P in SWEEP_PRESSURES:
            raw = {bg: run(bprime_path, bl, P * ONEATM, bg, SWEEP_T_RANGE)
                   for bg in BG_SWEEP}
            Ts = raw[0.0][:, 0]
            rows = []
            for j, T in enumerate(Ts):
                bcs = np.array([raw[bg][j, 1] for bg in BG_SWEEP])
                hws = np.array([raw[bg][j, 2] for bg in BG_SWEEP])
                bc = solve_operating_point(np.array(BG_SWEEP), bcs, K_CORK)
                bg = min(K_CORK * bc, BG_SWEEP[-1])
                rows.append((T, bc, bg, np.interp(bg, BG_SWEEP, hws), bcs[0]))
            ss[(bl, P)] = np.array(rows)
            print(f"  point de fonctionnement {bl:9s} P = {P:g} atm")

    with open("cork_oat_bprime_steady_state.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["edge", "O2_C2H2", "P_atm", "T_K", "Bc_ss", "Bg_ss",
                    "hw_ss_MJkg", "Bc_Bg0", "recession_over_mdote_m3_per_kg"])
        for bl, r, _ in EDGES:
            for P in SWEEP_PRESSURES:
                for T, bc, bg, hw, bc0 in ss[(bl, P)]:
                    w.writerow([bl, "" if r is None else f"{r:g}", f"{P:g}",
                                f"{T:g}", f"{bc:.6e}", f"{bg:.6e}",
                                f"{hw:.6e}", f"{bc0:.6e}",
                                f"{bc / RHO_CHAR:.6e}"])
    print("Point de fonctionnement : cork_oat_bprime_steady_state.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(r"Liège/phénolique P50 — $B'_c$ stationnaire "
                 r"($B'_g = 4\,B'_c$) : air vs torche oxyacétylénique",
                 fontsize=13)
    for ax, P in zip(axes, SWEEP_PRESSURES):
        for bl, r, lab in EDGES:
            d = ss[(bl, P)]
            ok = d[:, 2] < BG_SWEEP[-1]
            ax.plot(d[ok, 0], np.maximum(d[ok, 1], 1e-5), color=colors[bl],
                    lw=2, label=lab)
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 3)
        ax.set_xlabel(r"$T_w$ [K]")
        ax.set_ylabel(r"$B'_c$ stationnaire")
        ax.set_title("1 atm" if P == 1.0 else f"{P:g} atm")
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=9, loc="upper left")
    plt.tight_layout()
    plt.savefig("cork_oat_bprime_steady_state.png", dpi=150)
    plt.close()

    # --- Récapitulatif ----------------------------------------------------
    print("\nPlateaux théoriques : " + ", ".join(
        f"r = {r:g} -> {plateau(r):.4f}" for _, r, _ in EDGES if r))
    print(f"\n1 atm   {'':9s}" + "".join(f"{T:>10d}" for T in
                                         (1000, 1500, 2000, 2500, 3000, 3400)))
    for kind, bg in (("B'c B'g=0", 0.0), ("B'c B'g=5", 5.0)):
        for bl, _, _ in EDGES:
            d = tables[(bl, P1, bg)]
            vals = [np.interp(T, d[:, 0], d[:, 1])
                    for T in (1000, 1500, 2000, 2500, 3000, 3400)]
            print(f"{kind} {bl:9s}" + "".join(f"{v:>10.4f}" for v in vals))
    for bl, _, _ in EDGES:
        d = ss[(bl, 1.0)]
        vals = [np.interp(T, d[:, 0], d[:, 1])
                for T in (1000, 1500, 2000, 2500, 3000, 3400)]
        hws = [np.interp(T, d[:, 0], d[:, 3])
               for T in (1000, 1500, 2000, 2500, 3000, 3400)]
        print(f"B'c ss    {bl:9s}" + "".join(f"{v:>10.4f}" for v in vals))
        print(f"hw ss     {bl:9s}" + "".join(f"{v:>10.2f}" for v in hws))


if __name__ == "__main__":
    main()
