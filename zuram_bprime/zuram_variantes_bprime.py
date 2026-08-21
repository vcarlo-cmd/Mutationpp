#!/usr/bin/env python
"""
Variantes ZURAM XX/YY : vérification numérique par `bprime`.

Complète `zuram_variantes.py`, qui fait le bilan de phase solide. Ici on
appelle réellement le solveur pour établir deux choses :

  1. LA TABLE B' EST LA MÊME pour toutes les variantes. On ne se contente pas
     de réutiliser le même XML : pour chaque variante on RECALCULE la
     composition élémentaire du gaz par fermeture (résine - char) à partir de
     son propre bilan de masse, on écrit un mélange temporaire, et on compare
     les tables point par point. L'écart doit être rigoureusement nul.

  2. LE POINT DE FONCTIONNEMENT, LUI, CHANGE. En ablation stationnaire le
     matériau impose B'g = k·B'c avec k = (rho_v - rho_c)/rho_c. On résout ce
     point fixe pour chaque variante et on compare B'c, h_w et la vitesse de
     recul.

Usage :
    python zuram_variantes_bprime.py
    python zuram_variantes_bprime.py 14/40 18/50 18/80
"""

import os
import shutil
import subprocess
import sys

from zuram_variantes import (CHAR_YIELD, RHO_FIBER, RHO_RESIN, parse, variant)

ONEATM = 101325.0
M = {"C": 12.0107, "H": 1.00794, "O": 15.999, "N": 14.007}

# Résine ZURAM mesurée (Torres-Herrador fig. 7), renormalisée sur C+H+N+O.
# Identique dans toutes les variantes : c'est le même polymère.
RESIN = {"C": 77.91, "H": 5.96, "N": 1.49, "O": 14.64}

SPECIES = ("C(gr) C H O N CH CH4 CO CO2 CN C2H C2H2,acetylene C3 C3H C4 "
           "C4H2,butadiyne C5 HCN H2 H2O N2")


def pyro_gas_from_balance(v):
    """
    Composition élémentaire du gaz, recalculée depuis le bilan de masse de la
    variante : sur une base de 1 m³, la résine perd m_r - m_rc de masse, dont
    tout le carbone manquant part avec le gaz.

    Le résultat ne doit PAS dépendre de XX ni de YY — c'est le point à montrer.
    """
    basis = v["m_r"]                       # kg de résine vierge par m³
    mass = {e: basis * RESIN[e] / 100.0 for e in RESIN}
    m_char = v["m_rc"]                     # kg de char, supposé carbone pur
    gas = dict(mass)
    gas["C"] -= m_char
    if gas["C"] < 0:
        raise ValueError("rendement en char trop élevé pour cette résine")
    n = {e: gas[e] / M[e] for e in ("C", "H", "O")}
    tot = sum(n.values())
    return {e: n[e] / tot for e in n}


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


def run_bprime(path, mixdir, name, bg, pressure, trange):
    env = os.environ.copy()
    env["MPP_DATA_DIRECTORY"] = mixdir
    cmd = [path, "-T", trange, "-P", str(pressure), "-b", f"{bg:.10f}",
           "-m", name, "-bl", "air", "-py", "pyro",
           "-char", "char", "-char-elem", "C"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"ERREUR bprime :\n{r.stderr}")
        sys.exit(1)
    out = {}
    for line in r.stdout.strip().splitlines()[1:]:
        p = line.split()
        try:
            out[round(float(p[0]))] = (float(p[1]), float(p[2]))
        except (ValueError, IndexError):
            continue
    return out


def write_mixture(mixdir, name, gas):
    comp = ", ".join(f"{e}:{v:.12f}" for e, v in gas.items())
    xml = f"""<mixture thermo_db="NASA-9">
    <species> {SPECIES} </species>
    <element_compositions default="air">
        <composition name="air">O:0.21, N:0.79</composition>
        <composition name="pyro">{comp}</composition>
        <composition name="char">C:1.0</composition>
    </element_compositions>
</mixture>
"""
    with open(os.path.join(mixdir, "mixtures", f"{name}.xml"), "w") as f:
        f.write(xml)


def main():
    args = sys.argv[1:] or ["14/40", "18/50", "18/80"]
    variants = [variant(*parse(a)) for a in args]
    labels = [f"{int(v['xx'])}/{int(v['yy'])}" for v in variants]

    path = find_bprime()
    if path is None:
        print("binaire `bprime` introuvable — compiler d'abord "
              "(make -C build bprime)")
        return 1
    here = os.path.dirname(os.path.abspath(__file__))
    mixdir = os.path.abspath(os.path.join(here, "../data"))

    print()
    print("=" * 79)
    print("  1. LE GAZ DE PYROLYSE NE DÉPEND PAS DE XX/YY")
    print("=" * 79)
    print("  Recalculé pour chaque variante par fermeture résine - char,")
    print("  depuis son propre bilan de masse (kg/m³) :\n")
    print(f"  {'variante':>10} {'résine':>10} {'char':>10} {'gaz':>10} | "
          f"{'C':>9} {'H':>9} {'O':>9}")
    gases = []
    for lab, v in zip(labels, variants):
        g = pyro_gas_from_balance(v)
        gases.append(g)
        print(f"  {lab:>10} {v['m_r']:10.1f} {v['m_rc']:10.1f} "
              f"{v['m_r']-v['m_rc']:10.1f} | "
              f"{g['C']:9.6f} {g['H']:9.6f} {g['O']:9.6f}")
    spread = max(max(abs(g[e] - gases[0][e]) for e in g) for g in gases)
    print(f"\n  dispersion maximale entre variantes : {spread:.1e}")
    print(f"  -> {'IDENTIQUE' if spread < 1e-12 else 'DIVERGENCE'} : les masses "
          f"changent, les proportions non.")

    # --- 2. tables B' ------------------------------------------------------
    print("\n" + "=" * 79)
    print("  2. TABLES B' COMPARÉES POINT PAR POINT")
    print("=" * 79)
    trange = "500:50:4000"
    bgs = [0.0, 0.1, 0.5, 1.0]
    tmp = os.path.join(mixdir, "mixtures", "tmp_zuram_variant.xml")
    tables = []
    try:
        for g in gases:
            write_mixture(mixdir, "tmp_zuram_variant", g)
            tables.append({bg: run_bprime(path, mixdir, "tmp_zuram_variant",
                                          bg, ONEATM, trange) for bg in bgs})
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    print(f"  1 atm, T de 500 à 4000 K par 50 K, B'g = {bgs}")
    print(f"  {len(bgs) * len(tables[0][bgs[0]])} points par variante\n")
    print(f"  {'variante':>10} | {'ecart max Bc':>16} {'ecart max h_w':>16}")
    ref = tables[0]
    worst = 0.0
    for lab, tab in zip(labels, tables):
        dc = dh = 0.0
        for bg in bgs:
            for T, (bc, hw) in tab[bg].items():
                rbc, rhw = ref[bg][T]
                dc = max(dc, abs(bc - rbc))
                dh = max(dh, abs(hw - rhw))
        worst = max(worst, dc)
        print(f"  {lab:>10} | {dc:16.3e} {dh:16.3e}")
    print(f"\n  -> {'TABLE RIGOUREUSEMENT IDENTIQUE' if worst == 0.0 else 'DIVERGENCE'}"
          f" (référence : {labels[0]})")

    # --- 3. points de fonctionnement ---------------------------------------
    print("\n" + "=" * 79)
    print("  3. LE POINT DE FONCTIONNEMENT, LUI, CHANGE")
    print("=" * 79)
    print("  Ablation stationnaire : le matériau impose B'g = k·B'c.")
    print("  Point fixe résolu par itérations successives sur B'g.\n")

    write_mixture(mixdir, "tmp_zuram_variant", gases[0])
    temps = (1500, 2000, 2500, 3000, 3300, 3600)
    try:
        def fixed_point(v, T):
            """Résout B'g = k·B'c(B'g) à T et 1 atm. Retourne (B'c, B'g, h_w)."""
            bg, bc, hw = 0.1, 0.0, 0.0
            for _ in range(60):
                tab = run_bprime(path, mixdir, "tmp_zuram_variant", bg,
                                 ONEATM, "%d:50:%d" % (T, T + 50))
                bc, hw = tab.get(T, (0.0, 0.0))
                new = v["k"] * bc
                if abs(new - bg) < 1e-9:
                    bg = new
                    break
                bg = 0.5 * (bg + new)
            return bc, bg, hw

        head = "".join(("%s  k=%.3f" % (lab, v["k"])).rjust(24)
                       for lab, v in zip(labels, variants))
        print("  %7s |%s" % ("T [K]", head))
        sub = "".join("%8s%8s%8s" % ("B'c", "B'g", "recul") for _ in variants)
        print("  %7s |%s" % ("", sub))
        print("  " + "-" * (9 + 24 * len(variants)))
        for T in temps:
            sol = [fixed_point(v, T) for v in variants]
            # recul relatif : s_dot = B'c/rho_c, rapporté à la 1re variante
            ref = sol[0][0] / variants[0]["rho_c"]
            cells = "".join(
                "%8.4f%8.4f%8.2f" % (bc, bg, (bc / v["rho_c"]) / ref)
                for (bc, bg, _), v in zip(sol, variants))
            print("  %7d |%s" % (T, cells))
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    print("\n  s_dot rel : vitesse de recul rapportée à celle de la 1re variante")
    print("  à flux de masse d'arête identique (s_dot = B'c·m_e/rho_c).")
    print("\n" + "=" * 79)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
