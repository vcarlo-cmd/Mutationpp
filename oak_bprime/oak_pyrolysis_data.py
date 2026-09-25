#!/usr/bin/env python
"""
Mise en donnees d'un bois de chene : du materiau aux trois compositions
elementaires que consomme `bprime`.

Materiau considere
------------------
    bois de chene massif (Quercus robur / Q. petraea), SEC, sans resine
    rendement en char : 20 % (masse) -- donnee de l'enonce

C'est le cas le plus simple des ablateurs pyrolysants du depot : UN SEUL
constituant, qui pyrolyse entierement vers char + gaz. Pas de melange de gaz
comme pour le liege/phenolique (cf. cork_bprime/cork_pyrolysis_data.py), pas
de renfort inerte comme pour le TACOT.

Principe : fermeture elementaire (aucun atome n'est cree)

        bois  =  char  +  gaz
        n_E(gaz)  =  n_E(bois) - n_E(char)

Usage :
    python oak_pyrolysis_data.py
"""

import csv
import os
import re

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}

# ---------------------------------------------------------------------------
# Chene : composition BIOCHIMIQUE et unites de repetition
# ---------------------------------------------------------------------------
#   (formule de l'unite, % masse)
#
# Choix des unites -- chacun est un representant, pas une verite :
#   cellulose      : unite anhydroglucose C6H10O5
#   hemicelluloses : glucuronoxylane des feuillus -> unite anhydroxylose
#                    C5H8O4 (acetyles et acide 4-O-methylglucuronique negliges)
#   lignine        : lignine de feuillu syringyle / guaiacyle, S/G ~ 1.5
#                    molaire. Deux unites :
#                      syringyle  = alcool sinapylique   C11H14O4 (210.2)
#                      guaiacyle  = alcool coniferylique C10H12O3 (180.2)
#                    16 % S + 9 % G en masse -> S/G = 1.52 molaire
#   tanins         : ellagitanins du bois de coeur (vescalagine / castalagine,
#                    isomeres) -> C41H26O26
OAK_BIOCHEM = {
    "cellulose":      ({"C":  6, "H": 10, "O":  5}, 43.0),
    "hemicelluloses": ({"C":  5, "H":  8, "O":  4}, 22.0),
    "lignine S":      ({"C": 11, "H": 14, "O":  4}, 16.0),
    "lignine G":      ({"C": 10, "H": 12, "O":  3},  9.0),
    "tanins":         ({"C": 41, "H": 26, "O": 26},  8.0),
}
# Les parts somment a 98 % : le complement (cendres ~0.4 %, extractibles
# mineurs, groupes acetyles) n'est pas de la matiere C/H/O identifiee. On
# renormalise sur les 98 % declares.

# Analyse elementaire du chene, valeurs de litterature (sec, sans cendres).
# C'est la donnee RETENUE ; la reconstruction depuis les unites de repetition
# (oak_elemental) sert de controle.
#
# ATTENTION -- PROVENANCE A CONFIRMER. Ces trois nombres sont des valeurs
# usuelles citees pour le bois de chene (C ~50, H ~6, O ~44 % masse) ; ils
# n'ont PAS ete verifies sur une source primaire accessible depuis cet
# environnement. L'azote (0.1-0.3 %) est neglige.
OAK_MASS_PCT_LITERATURE = {"C": 50.0, "H": 6.1, "O": 43.9}

# Rendements en char des constituants (ordre de grandeur, ATG lente sous
# inerte) : sert de CONTROLE au rendement de 20 % de l'enonce.
OAK_CONSTITUENT_CHAR_YIELD = {
    "cellulose":      0.10,   # depolymerisation en levoglucosan, peu de char
    "hemicelluloses": 0.20,
    "lignine S":      0.35,   # la lignine S charbonne moins que la G
    "lignine G":      0.45,
    "tanins":         0.40,   # polyphenols
}

# Rendement en char du bois : DONNEE DE L'ENONCE.
OAK_CHAR_YIELD = 0.20

# Char suppose purement carbone (cas de base).
CHAR_COMP = {"C": 1.0}

# Masse volumique (kg/m3) pour la reponse materiau -- n'entre PAS dans le XML,
# seulement dans la recession s_dot = B'c mdot_e / rho_c.
# rho_vierge : chene sec, ordre de grandeur 650-750 kg/m3.
# rho_char   : NON MESUREE ici. Le bois se retracte en carbonisant ; on fait
#              l'hypothese V_char/V_vierge = 0.55, a remplacer par une mesure.
RHO_VIRGIN = 700.0
VOLUME_RATIO_CHAR = 0.55
RHO_CHAR = RHO_VIRGIN * OAK_CHAR_YIELD / VOLUME_RATIO_CHAR


# ---------------------------------------------------------------------------
# Fermeture elementaire
# ---------------------------------------------------------------------------

def molar_mass(comp):
    return sum(n * M[e] for e, n in comp.items())


def normalize(comp):
    tot = sum(comp.values())
    return {e: n / tot for e, n in comp.items()}


def mass_fractions(comp):
    m = {e: n * M[e] for e, n in comp.items()}
    tot = sum(m.values())
    return {e: v / tot for e, v in m.items()}


def oak_elemental(biochem=None):
    """Analyse elementaire (% masse) du chene depuis ses unites de repetition."""
    biochem = biochem or OAK_BIOCHEM
    total = sum(w for _, w in biochem.values())
    out = {"C": 0.0, "H": 0.0, "O": 0.0}
    for unit, w in biochem.values():
        y = mass_fractions(unit)
        for e in out:
            out[e] += (w / total) * y.get(e, 0.0)
    return {e: 100.0 * v for e, v in out.items()}


def oak_char_yield_from_constituents(biochem=None, yields=None):
    """Rendement en char du bois par additivite des constituants."""
    biochem = biochem or OAK_BIOCHEM
    yields = yields or OAK_CONSTITUENT_CHAR_YIELD
    total = sum(w for _, w in biochem.values())
    return sum((w / total) * yields[name] for name, (_, w) in biochem.items())


def repeat_unit(mass_pct, n_carbon=6):
    """Unite de repetition globale C_n H_x O_y equivalente a une analyse."""
    n = {e: mass_pct[e] / M[e] for e in mass_pct}
    return {e: n_carbon * n[e] / n["C"] for e in n}


# Analyse elementaire retenue : celle de la LITTERATURE.
# Mettre OAK_MASS_PCT = oak_elemental() pour passer a la reconstruction.
OAK_MASS_PCT = dict(OAK_MASS_PCT_LITERATURE)


def moles_from_mass_pct(mass_pct, mass):
    """Analyse elementaire (% masse) + masse [g] -> moles d'atomes."""
    y = normalize(mass_pct)
    return {e: mass * y[e] / M[e] for e in y}


def split_char_gas(constituent_moles, char_yield, char_comp=None):
    """Retourne (moles_char, moles_gaz) pour un rendement en char massique."""
    char_comp = normalize(char_comp or CHAR_COMP)
    m_virgin = sum(n * M[e] for e, n in constituent_moles.items())
    m_char = char_yield * m_virgin
    m_per_atom = sum(char_comp[e] * M[e] for e in char_comp)
    n_atoms = m_char / m_per_atom
    char = {e: n_atoms * char_comp[e] for e in char_comp}

    gas = {}
    for e in set(constituent_moles) | set(char):
        v = constituent_moles.get(e, 0.0) - char.get(e, 0.0)
        if v < -1e-9:
            raise ValueError(
                f"fermeture impossible : le char demande plus de {e} que le "
                f"bois n'en contient (rendement trop eleve ?)")
        gas[e] = max(v, 0.0)
    return char, gas


def fmt(comp, digits=3):
    return ", ".join(f"{e}:{v:.{digits}f}"
                     for e, v in sorted(comp.items()) if v > 1e-12)


def wood_balance(char_yield=OAK_CHAR_YIELD, char_comp=None, mass_pct=None,
                 moisture=0.0, basis=100.0):
    """Bilan sur `basis` grammes de bois vierge.

    moisture : fraction massique d'eau du bois tel que charge (base humide).
               L'eau part integralement dans le gaz ; le rendement en char
               reste rapporte au bois SEC.
    """
    wood = moles_from_mass_pct(mass_pct or OAK_MASS_PCT, basis)
    char, gas = split_char_gas(wood, char_yield, char_comp)
    m_char = sum(n * M[e] for e, n in char.items())
    m_gas = sum(n * M[e] for e, n in gas.items())
    if moisture:
        m_water = basis * moisture / (1.0 - moisture)
        n_w = m_water / (2 * M["H"] + M["O"])
        gas["H"] += 2 * n_w
        gas["O"] += n_w
        m_gas += m_water
        basis += m_water
    return dict(wood=wood, gas=gas, char=char, m_char=m_char, m_gas=m_gas,
                char_yield=m_char / basis, k=m_gas / m_char)


# ---------------------------------------------------------------------------
# Programme
# ---------------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    line = "=" * 76

    print(line)
    print("0. LE CHENE RECONSTRUIT DEPUIS SES UNITES DE REPETITION")
    print(line)
    tot = sum(w for _, w in OAK_BIOCHEM.values())
    print(f"  {'constituant':15s} {'% masse':>8}  {'unite':>12}  "
          f"{'M':>8}  {'C':>6} {'H':>6} {'O':>6}  (% masse de l'unite)")
    for name, (unit, w) in OAK_BIOCHEM.items():
        y = mass_fractions(unit)
        f = "".join(f"{e}{int(n)}" for e, n in unit.items())
        print(f"  {name:15s} {w:8.1f}  {f:>12}  {molar_mass(unit):8.2f}  "
              f"{100*y['C']:6.1f} {100*y['H']:6.1f} {100*y['O']:6.1f}")
    print(f"  {'':15s} {tot:8.1f}   <- complement a 100 % : cendres, "
          "acetyles, extractibles mineurs")

    rec = oak_elemental()
    lit = OAK_MASS_PCT_LITERATURE
    print(f"\n  chene reconstruit (parts renormalisees) : "
          f"C {rec['C']:.2f}  H {rec['H']:.2f}  O {rec['O']:.2f}")
    print(f"  analyse de litterature  <- RETENUE      : "
          f"C {lit['C']:.2f}  H {lit['H']:.2f}  O {lit['O']:.2f}")
    print(f"  ecart                                   : "
          f"C {rec['C']-lit['C']:+.2f}  H {rec['H']-lit['H']:+.2f}  "
          f"O {rec['O']-lit['O']:+.2f}")
    for lab, mp in (("litterature", lit), ("reconstruction", rec)):
        u = repeat_unit(mp)
        print(f"  unite de repetition globale ({lab:14s}) : "
              f"C6 H{u['H']:.2f} O{u['O']:.2f}  "
              f"(M = {molar_mass(u):.1f} g/mol ; CH{u['H']/6:.3f}O{u['O']/6:.3f})")
    b = wood_balance()
    br = wood_balance(mass_pct=rec)
    print(f"  Effet sur le gaz de pyrolyse (a char inchange) :")
    print(f"     avec la litterature (retenue) : {fmt(normalize(b['gas']))}")
    print(f"     avec la reconstruction        : {fmt(normalize(br['gas']))}")

    print("\n" + line)
    print("0'. CONTROLE DU RENDEMENT EN CHAR PAR ADDITIVITE DES CONSTITUANTS")
    print(line)
    add = oak_char_yield_from_constituents()
    for name, (_, w) in OAK_BIOCHEM.items():
        print(f"  {name:15s} {w:6.1f} %   char {100*OAK_CONSTITUENT_CHAR_YIELD[name]:3.0f} %")
    print(f"  -> chene : {100*add:.1f} %   (contre {100*OAK_CHAR_YIELD:.1f} %"
          " de l'enonce) : les deux voies sont coherentes.")

    print("\n" + line)
    print("CHENE SEC -- fermeture elementaire, rendement en char "
          f"{100*OAK_CHAR_YIELD:.0f} %")
    print(line)
    print(f"\n1. Bois vierge (base 100 g)")
    print(f"   moles : {fmt(b['wood'], 4)}   soit x = {fmt(normalize(b['wood']))}")
    print(f"\n2. Char : {fmt(b['char'], 4)} mol   ({b['m_char']:.2f} g, C pur)")
    print(f"\n3. Gaz de pyrolyse = bois - char")
    print(f"   moles      : {fmt(b['gas'], 4)}   ({b['m_gas']:.2f} g)")
    print(f"   x molaires : {fmt(normalize(b['gas']))}")
    print(f"   y massiques: {fmt(mass_fractions(b['gas']))}")
    x = normalize(b["gas"])
    print(f"   H/O = {x['H']/x['O']:.3f}   M moyen par atome = "
          f"{sum(x[e]*M[e] for e in x):.3f} g/mol")
    print(f"\n4. Couplage stationnaire")
    print(f"   k = B'g/B'c = m_gaz/m_char = {b['k']:.3f}")
    print(f"   rho_v = {RHO_VIRGIN:.0f} kg/m3 ; rho_c = {RHO_CHAR:.0f} kg/m3 "
          f"(hypothese V_char/V_v = {VOLUME_RATIO_CHAR})")
    print(f"   (rho_v - rho_c)/rho_c = {(RHO_VIRGIN-RHO_CHAR)/RHO_CHAR:.2f} :"
          " FAUX des que le char se retracte -> k sur les MASSES.")
    print(f"\n5. A recopier dans data/mixtures/oak-air.xml")
    print(f'   <composition name="oak_pyro">{fmt(normalize(b["gas"]))}'
          "</composition>")
    print(f'   <composition name="oak_char">{fmt(normalize(b["char"]), 1)}'
          "</composition>")

    print("\n" + line)
    print("SENSIBILITE AU RENDEMENT EN CHAR")
    print(line)
    print(f"  {'char':>6} | {'C':>7} {'H':>7} {'O':>7} | {'k':>6}")
    rows = []
    for cy in (0.10, 0.15, 0.20, 0.25, 0.30):
        bb = wood_balance(char_yield=cy)
        xx = normalize(bb["gas"])
        print(f"  {100*cy:>5.0f}% | {xx['C']:>7.3f} {xx['H']:>7.3f} "
              f"{xx['O']:>7.3f} | {bb['k']:>6.3f}")
        rows.append(dict(char_yield=cy, x_C=xx["C"], x_H=xx["H"], x_O=xx["O"],
                         k=bb["k"]))

    print("\n" + line)
    print("VARIANTE : BOIS HUMIDE (l'eau part dans le gaz)")
    print(line)
    for mo in (0.0, 0.08, 0.12):
        bb = wood_balance(moisture=mo)
        xx = normalize(bb["gas"])
        print(f"  humidite {100*mo:4.0f} % | gaz {fmt(xx)} | "
              f"char {100*bb['char_yield']:.1f} % du bois humide | k = {bb['k']:.3f}")

    print("\n" + line)
    print("VARIANTE : CHAR RETENANT DE L'OXYGENE ET DE L'HYDROGENE")
    print(line)
    cc = {e: v / M[e] for e, v in (("C", 0.90), ("H", 0.02), ("O", 0.08))}
    bv = wood_balance(char_comp=cc)
    print("  char 90 % C / 2 % H / 8 % O en masse")
    print(f"  char : {fmt(normalize(bv['char']))}")
    print(f"  gaz  : {fmt(normalize(bv['gas']))}   (contre "
          f"{fmt(normalize(b['gas']))} avec un char C pur)")

    out = os.path.join(here, "oak_pyrolysis_data.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n  -> {out}")


if __name__ == "__main__":
    main()
