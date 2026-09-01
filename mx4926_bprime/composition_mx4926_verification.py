#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rejoue toutes les vérifications de `composition_mx4926.md`.

Pendant de `sc1008_bprime/resine_sc1008_verification.py`, mais un cran plus
haut : la résine y était l'objet, ici c'est le COMPOSITE. La composition
élémentaire du gaz n'est donc pas recalculée — elle est héritée du SC-1008 —
et ce script contrôle ce qui est propre au MX-4926 :

  1. fermeture du domaine de composition  (les trois plages ne sont pas
     indépendantes : bornes basses 83 %, bornes hautes 109 %)
  2. centroïde EXACT du domaine admissible (centroïde du polygone, pas
     une moyenne de grille)
  3. rendement en char du composite et couplage k = B'g/B'c
  4. densités vierge et char, et porosité qui restitue la densité publiée
  5. identité BIT À BIT de la table B' du MX-4926 et de celle du SC-1008
  6. contrôle physique : B'c(300 K, B'g = 0, air) = 0.0874

Les contrôles 5 et 6 nécessitent le binaire `bprime` ; ils sont sautés (et
signalés) s'il est introuvable.

Usage :
    python composition_mx4926_verification.py
"""

import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# Données d'entrée — les seules valeurs en dur de ce script
# ---------------------------------------------------------------------------

# Plages de fractions MASSIQUES du composite (données utilisateur)
RANGES = {
    "fibre":  (0.41, 0.56),   # renfort carbone ex-rayon, satin de 5 ou 8
    "resine": (0.31, 0.37),   # résol phénol-formaldéhyde SC-1008
    "charge": (0.11, 0.16),   # noir de carbone
}

CHAR_YIELD_RESIN = 0.55       # sc1008_bprime/resine_sc1008.md §5.2

RHO = {                       # masses volumiques intrinsèques [kg/m3]
    "fibre":  1600.0,         # fibres ex-cellulose — Description!A11, TACOT 3.0
    "resine": 1250.0,         # résol cuit
    "charge": 1800.0,         # noir de carbone
}

RHO_MX4926_PUBLIE = 1450.0    # densité publiée du MX-4926, ~1.45 g/cm3

# Gaz de pyrolyse hérité du SC-1008 (fractions molaires élémentaires)
PYRO = {"C": 0.2526, "H": 0.6407, "O": 0.1068}

ONEATM = 101325.0
TOL = 5e-4


# ---------------------------------------------------------------------------
# Géométrie du domaine admissible
# ---------------------------------------------------------------------------

def clip(poly, keep):
    """Sutherland-Hodgman : ne garde du polygone que le demi-plan keep(p) >= 0."""
    out = []
    n = len(poly)
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        fa, fb = keep(a), keep(b)
        if fa >= -1e-15:
            out.append(a)
        if fa * fb < -1e-15:
            t = fa / (fa - fb)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out


def admissible_polygon():
    """
    Domaine admissible dans le plan (w_fibre, w_resine).

    La boîte des trois plages, intersectée avec la contrainte de fermeture
    Sum w = 1 réécrite sur la charge :

        w_charge_min <= 1 - w_fibre - w_resine <= w_charge_max
    """
    (fl, fh), (rl, rh), (cl, ch) = (RANGES["fibre"], RANGES["resine"],
                                    RANGES["charge"])
    poly = [(fl, rl), (fh, rl), (fh, rh), (fl, rh)]
    poly = clip(poly, lambda p: (1 - p[0] - p[1]) - cl)   # charge >= min
    poly = clip(poly, lambda p: ch - (1 - p[0] - p[1]))   # charge <= max
    return poly


def polygon_centroid(poly):
    """Centroïde d'aire exact (formule du polygone), et son aire."""
    a = cx = cy = 0.0
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        cr = x0 * y1 - x1 * y0
        a += cr
        cx += (x0 + x1) * cr
        cy += (y0 + y1) * cr
    a /= 2.0
    return (cx / (6 * a), cy / (6 * a)), a


# ---------------------------------------------------------------------------
# Réponse matériau
# ---------------------------------------------------------------------------

def char_yield_composite(w_resine, y=CHAR_YIELD_RESIN):
    """
    Rendement en char du COMPOSITE.

        Y_comp = w_fibre + w_charge + Y * w_resine
               = (1 - w_resine) + Y * w_resine
               = 1 - (1 - Y) * w_resine

    Le renfort ex-rayon est déjà carbonisé et le noir de carbone est inerte :
    ni l'un ni l'autre ne perd de masse. La fermeture Sum w = 1 fait alors
    disparaître w_fibre et w_charge — Y_comp ne dépend QUE de la fraction de
    résine, pas du partage fibres/charge.
    """
    return 1.0 - (1.0 - y) * w_resine


def coupling(y_comp):
    """k = B'g/B'c = (rho_v - rho_c)/rho_c = (1 - Y_comp)/Y_comp."""
    return (1.0 - y_comp) / y_comp


def densities(w, porosity):
    """rho vierge et rho char [kg/m3] pour une composition massique w."""
    v = sum(w[k] / RHO[k] for k in w)          # volume de solide par kg
    rho_v = (1.0 - porosity) / v
    return rho_v, rho_v * char_yield_composite(w["resine"])


def porosity_for(w, rho_target):
    """Porosité qui restitue une densité vierge donnée."""
    return 1.0 - rho_target * sum(w[k] / RHO[k] for k in w)


# ---------------------------------------------------------------------------
# bprime
# ---------------------------------------------------------------------------

def find_bprime():
    c = shutil.which("bprime")
    if c:
        return c
    here = os.path.dirname(os.path.abspath(__file__))
    for p in (os.path.join(here, "../build/src/apps/bprime"),
              os.path.join(here, "../../build/src/apps/bprime"),
              "build/src/apps/bprime"):
        if os.path.isfile(p):
            return os.path.abspath(p)
    return None


def run_bprime(path, mixture, pyro, char, bg, t_range="300:25:4000",
               pressure=ONEATM):
    here = os.path.dirname(os.path.abspath(__file__))
    env = os.environ.copy()
    env.setdefault("MPP_DATA_DIRECTORY",
                   os.path.abspath(os.path.join(here, "../data")))
    cmd = [path, "-T", t_range, "-P", str(pressure), "-b", str(bg),
           "-m", mixture, "-bl", "air", "-py", pyro,
           "-char", char, "-char-elem", "C"]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        print("ERREUR bprime :\n" + res.stderr)
        sys.exit(1)
    return res.stdout


# ---------------------------------------------------------------------------
# Programme
# ---------------------------------------------------------------------------

def main():
    ko = 0
    line = "=" * 78

    # --- 1. fermeture du domaine ------------------------------------------
    print(line)
    print("  1. FERMETURE DU DOMAINE DE COMPOSITION")
    print(line)
    lo = sum(v[0] for v in RANGES.values())
    hi = sum(v[1] for v in RANGES.values())
    print("  plages annoncées (fractions massiques) :")
    for k, (a, b) in RANGES.items():
        print(f"     {k:7s} : {100*a:5.1f} – {100*b:5.1f} %")
    print(f"\n  somme des bornes basses : {100*lo:5.1f} %")
    print(f"  somme des bornes hautes : {100*hi:5.1f} %")
    ok = lo < 1.0 < hi
    print(f"  [{'OK ' if ok else 'KO '}] les trois plages ne sont donc PAS "
          "indépendantes : toute composition")
    print("       réelle doit sommer à 100 %, ce qui restreint le domaine.")
    ko += 0 if ok else 1

    poly = admissible_polygon()
    print(f"\n  domaine admissible : polygone à {len(poly)} sommets")
    for x, y in poly:
        print(f"     fibre {100*x:5.2f} %   résine {100*y:5.2f} %   "
              f"charge {100*(1-x-y):5.2f} %")

    eff = {
        "fibre":  (min(p[0] for p in poly), max(p[0] for p in poly)),
        "resine": (min(p[1] for p in poly), max(p[1] for p in poly)),
        "charge": (min(1 - p[0] - p[1] for p in poly),
                   max(1 - p[0] - p[1] for p in poly)),
    }
    print("\n  plages EFFECTIVES après fermeture :")
    for k in RANGES:
        a, b = eff[k]
        tag = "" if abs(a - RANGES[k][0]) < 1e-9 and abs(b - RANGES[k][1]) < 1e-9 \
              else "  <-- tronquée"
        print(f"     {k:7s} : {100*a:5.2f} – {100*b:5.2f} %{tag}")
    trunc = abs(eff["fibre"][0] - RANGES["fibre"][0]) > 1e-9
    print(f"\n  [{'OK ' if trunc else 'KO '}] seule la borne basse du renfort "
          f"est déplacée : {100*RANGES['fibre'][0]:.0f} % -> "
          f"{100*eff['fibre'][0]:.0f} %.")
    print("       Les 41 % annoncés ne sont atteignables avec AUCUN couple")
    print("       (résine, charge) admissible.")
    ko += 0 if trunc else 1

    # --- 2. centroïde exact ------------------------------------------------
    print("\n" + line)
    print("  2. CENTROÏDE EXACT DU DOMAINE ADMISSIBLE")
    print(line)
    (wf, wr), area = polygon_centroid(poly)
    wc = 1.0 - wf - wr
    w = {"fibre": wf, "resine": wr, "charge": wc}
    print(f"  aire du polygone : {area:.6f}")
    print(f"  centroïde        : fibre {100*wf:.4f} %   résine {100*wr:.4f} % "
          f"  charge {100*wc:.4f} %")
    s = wf + wr + wc
    print(f"  somme            : {s:.12f}")
    inside = all(eff[k][0] - 1e-12 <= w[k] <= eff[k][1] + 1e-12 for k in w)
    print(f"  [{'OK ' if inside and abs(s-1) < 1e-12 else 'KO '}] "
          "il somme à 1 et tombe dans les plages effectives")
    ko += 0 if inside and abs(s - 1) < 1e-12 else 1
    print("\n  Valeurs arrondies portées dans l'en-tête du XML : "
          f"{100*wf:.1f} / {100*wr:.1f} / {100*wc:.1f} %")

    # --- 3. rendement en char et couplage ---------------------------------
    print("\n" + line)
    print("  3. RENDEMENT EN CHAR DU COMPOSITE ET COUPLAGE k")
    print(line)
    print("  Y_comp = w_fibre + w_charge + Y*w_resine")
    print("         = 1 - (1 - Y)*w_resine          (par la fermeture Sum w = 1)")
    print(f"\n  [!] Y_comp ne dépend QUE de w_resine : le partage fibres/charge")
    print("      n'a aucun effet, les deux étant du carbone à perte de masse nulle.")

    # contrôle : la forme longue et la forme courte coïncident
    longue = w["fibre"] + w["charge"] + CHAR_YIELD_RESIN * w["resine"]
    courte = char_yield_composite(w["resine"])
    okf = abs(longue - courte) < 1e-15
    print(f"\n  forme longue  : {longue:.10f}")
    print(f"  forme courte  : {courte:.10f}")
    print(f"  [{'OK ' if okf else 'KO '}] les deux écritures coïncident")
    ko += 0 if okf else 1

    yc = courte
    k = coupling(yc)
    print(f"\n  au centroïde  : Y_comp = {yc:.4f}    k = {k:.4f}")
    print("\n  enveloppe sur le domaine (pilotée par w_resine seule) :")
    print(f"  {'w_resine':>10} {'Y_comp':>10} {'k':>10}")
    for wri in (RANGES["resine"][0], wr, RANGES["resine"][1]):
        yi = char_yield_composite(wri)
        tag = "   <- centroïde" if abs(wri - wr) < 1e-12 else ""
        print(f"  {wri:>10.4f} {yi:>10.4f} {coupling(yi):>10.4f}{tag}")

    print("\n  comparaison aux autres matériaux du dépôt :")
    for lab, kk, note in [("TACOT   (10/10/80 vol.)", 0.2727, "tacot_bprime/material_response.py"),
                          ("PICA    (SC-1008)",       0.2346, "sc1008 §7 : 0.19/0.81"),
                          ("MX-4926 (52/34/14 wt)",   k,      "ce script"),
                          ("CPh70   (69/30/1 vol.)",  0.1385, "cph70_bprime/README.md")]:
        print(f"     {lab:26s} k = {kk:.4f}   ({note})")

    # --- 4. densités -------------------------------------------------------
    print("\n" + line)
    print("  4. DENSITÉS ET POROSITÉ")
    print(line)
    print("  masses volumiques intrinsèques retenues [kg/m3] :")
    for kk, v in RHO.items():
        print(f"     {kk:7s} : {v:7.1f}")
    phi = porosity_for(w, RHO_MX4926_PUBLIE)
    print(f"\n  porosité qui restitue la densité publiée "
          f"({RHO_MX4926_PUBLIE:.0f} kg/m3) : {phi:.4f}")
    okphi = 0.0 <= phi <= 0.10
    print(f"  [{'OK ' if okphi else 'KO '}] elle est faible et positive — "
          "cohérent avec un composite dense pressé")
    ko += 0 if okphi else 1

    print(f"\n  {'porosité':>10} {'rho_v':>10} {'rho_c':>10} {'gaz':>10} "
          f"{'% masse':>10}")
    for p in (0.00, 0.02, 0.04, 0.08):
        rv, rc = densities(w, p)
        print(f"  {p:>10.2f} {rv:>10.1f} {rc:>10.1f} {rv-rc:>10.1f} "
              f"{100*(1-yc):>9.1f} %")
    rv, rc = densities(w, 0.02)
    print(f"\n  Valeurs portées dans l'en-tête du XML (porosité 0.02) :")
    print(f"     rho_v = {rv:.0f} kg/m3   rho_c = {rc:.0f} kg/m3   "
          f"gaz = {rv-rc:.0f} kg/m3")
    print("\n  [!] k ne dépend PAS de ces densités : k = (1-Y_comp)/Y_comp est")
    print("      un rapport de MASSES. Les densités ne servent qu'à la")
    print("      récession s_dot = B'c*mdot_e/rho_c.")

    # --- 5 & 6. table B' ---------------------------------------------------
    print("\n" + line)
    print("  5. IDENTITÉ DE LA TABLE B' AVEC CELLE DU SC-1008")
    print(line)
    path = find_bprime()
    if path is None:
        print("  binaire `bprime` introuvable — contrôles 5 et 6 sautés.")
        print("    cmake -B build -DCMAKE_BUILD_TYPE=Release .")
        print("    cmake --build build --target bprime")
    else:
        print(f"  binaire : {path}")
        print("\n  Même liste d'espèces, même gaz, même char => la table doit")
        print("  être identique CARACTÈRE POUR CARACTÈRE.\n")
        for bg in (0.0, 0.1, 0.5, 2.0):
            a = run_bprime(path, "mx4926-air", "mx4926_pyro", "mx4926_char", bg)
            b = run_bprime(path, "sc1008-air", "sc1008_pyro", "sc1008_char", bg)
            same = a == b
            ko += 0 if same else 1
            print(f"  [{'OK ' if same else 'KO '}] B'g = {bg:<4} : "
                  f"{'tables identiques' if same else 'TABLES DIFFERENTES'} "
                  f"({len(a.splitlines())-1} lignes)")

        print("\n" + line)
        print("  6. CONTRÔLE PHYSIQUE : LIMITE D'OXYDATION À 300 K")
        print(line)
        out = run_bprime(path, "mx4926-air", "mx4926_pyro", "mx4926_char", 0.0,
                         t_range="300:25:400")
        bc = float(out.splitlines()[1].split()[1])
        okb = abs(bc - 0.0874) < 1e-3
        print(f"  B'c(300 K, B'g = 0, 1 atm) = {bc:.6f}")
        print(f"  [{'OK ' if okb else 'KO '}] attendu 0.0874 — limite "
              "C + O2 -> CO2, indépendante de la pression")
        print("       (une valeur ~200 signalerait l'absence de C(gr))")
        ko += 0 if okb else 1

        # la composition du gaz est bien celle du SC-1008
        print("\n  Rappel — composition héritée, non recalculée ici :")
        print("     mx4926_pyro = sc1008_pyro = "
              + ", ".join(f"{e}:{v}" for e, v in PYRO.items()))
        print(f"     H/O = {PYRO['H']/PYRO['O']:.2f}  (invariant phénolique, "
              "cf. resine_sc1008.md §3)")

    print("\n" + line)
    print("  [OK] toutes les vérifications passent" if ko == 0
          else f"  [KO] {ko} vérification(s) en échec")
    print(line)
    return 1 if ko else 0


if __name__ == "__main__":
    sys.exit(main())
