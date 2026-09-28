#!/usr/bin/env python
"""
Table B' du liège SEUL (sans résine), rendement en char imposé à 20 %.

Cas de comparaison, pas un matériau du dépôt : même analyse élémentaire que
le liège du P50 (C 62.4 / H 8.5 / O 28.4), même char (C pur), même air,
mais sans résine et avec le rendement en char du composite P50 et du chêne
(20 %), pour isoler l'effet de la chimie du gaz de pyrolyse.

    gaz (`corkpure_pyro`) : C:0.2576, H:0.6133, O:0.1291   (O/C = 0.50)
    k = B'g/B'c = (1 - y)/y = 4.0

La composition est recalculée ici par cork_pyrolysis_data.composite_balance
et vérifiée contre data/mixtures/corkpure-air.xml.

Reprend les fonctions de cork_bprime.py (appel de bprime, lecture, point
fixe stationnaire) avec le mélange corkpure-air.

Sorties (format identique aux fichiers du P50) :
    cork_seul_bprime_bc_table.csv     Tw_K, P_bar, Bg, Bc
    cork_seul_bprime_steady_state.csv P_atm, T_K, Bc_ss, Bg_ss, hw_ss_MJkg,
                                      Bc_Bg0, mass_loss_over_mdote
"""

import csv
import os
import re
import sys

import numpy as np

import cork_bprime as cb
import cork_pyrolysis_data as cpd

HERE = os.path.dirname(os.path.abspath(__file__))

Y_CHAR = 0.20
K = (1.0 - Y_CHAR) / Y_CHAR

cb.MIXTURE = "corkpure-air"
cb.PYRO_COMP = "corkpure_pyro"
cb.CHAR_COMP = "corkpure_char"


def check_composition():
    b = cpd.composite_balance(w_cork=1.0, w_resin=0.0, cork_char_yield=Y_CHAR)
    gas = cpd.normalize(b["gas"])
    xml = open(os.path.join(HERE, "../data/mixtures/corkpure-air.xml")).read()
    m = re.search(r'name="corkpure_pyro">([^<]+)<', xml)
    ref = {e: float(v) for e, v in
           (kv.split(":") for kv in m.group(1).replace(" ", "").split(","))}
    print(f"gaz liège seul, y = {Y_CHAR:.0%} : {cpd.fmt(gas, 4)}  "
          f"O/C = {gas['O'] / gas['C']:.3f}  k = {b['k']:.3f}")
    for e, v in ref.items():
        if abs(gas[e] - v) > 1e-4:
            sys.exit(f"corkpure-air.xml incohérent sur {e} : {v} vs {gas[e]:.4f}")


def main():
    check_composition()
    bprime = cb.find_bprime()
    if bprime is None:
        sys.exit("binaire bprime introuvable (cmake --build build --target bprime)")

    # table B'c(T, P, B'g), mêmes isobares et B'g que le P50
    out = os.path.join(HERE, "cork_seul_bprime_bc_table.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Tw_K", "P_bar", "Bg", "Bc"])
        for bg in cb.BG_VALUES:
            for P in cb.PRESSURES_ATM:
                _, d = cb.parse_output(cb.run_bprime(bprime, P * cb.ONEATM, bg))
                for row in d:
                    w.writerow([f"{row[0]:.6g}", f"{P * cb.ONEATM / 1e5:.6g}",
                                f"{bg:g}", f"{row[1]:.6e}"])
    print(f"-> {os.path.relpath(out, HERE)}")

    # point de fonctionnement stationnaire B'c = table(T, P, B'g = k B'c)
    bgs = np.array(cb.BG_SWEEP)
    out = os.path.join(HERE, "cork_seul_bprime_steady_state.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["P_atm", "T_K", "Bc_ss", "Bg_ss", "hw_ss_MJkg",
                    "Bc_Bg0", "mass_loss_over_mdote"])
        for P in cb.SWEEP_PRESSURES:
            raw = [cb.parse_output(cb.run_bprime(
                bprime, P * cb.ONEATM, bg, cb.SWEEP_T_RANGE))[1]
                for bg in cb.BG_SWEEP]
            for j, T in enumerate(raw[0][:, 0]):
                bcs = np.array([d[j, 1] for d in raw])
                hws = np.array([d[j, 2] for d in raw])
                bc = cb.solve_operating_point(bgs, bcs, K)
                bg = min(K * bc, bgs[-1])
                hw = np.interp(bg, bgs, hws)
                w.writerow([f"{P:g}", f"{T:g}", f"{bc:.6e}", f"{bg:.6e}",
                            f"{hw:.6e}", f"{bcs[0]:.6e}",
                            f"{(1 + K) * bc:.6e}"])
    print(f"-> {os.path.relpath(out, HERE)}")


if __name__ == "__main__":
    main()
