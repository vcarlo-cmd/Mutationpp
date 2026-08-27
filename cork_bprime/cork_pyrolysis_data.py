#!/usr/bin/env python
"""
Mise en donnees d'un liege/phenolique (cork phenolic) : du materiau aux trois
compositions elementaires que consomme `bprime`.

Materiau considere
------------------
    80 % liege / 20 % resine phenolique  (fractions MASSIQUES du solide)
    resine phenolique : rendement en char 50 % (masse)
    liege             : rendement en char 25 % (masse) -- le liege pyrolyse
                        AUSSI, contrairement aux fibres de carbone d'un
                        TACOT / PICA / CPh70.

C'est toute la difference avec un carbone/phenolique : les deux constituants
produisent du gaz, et leurs gaz n'ont pas la meme composition elementaire.
Le rapport renfort/resine entre donc DANS la table B', ce qui n'etait pas le
cas pour le CPh70 (cf. tacot_bprime/cph70_vs_tacot.md).

Principe : fermeture elementaire, constituant par constituant.

        constituant  =  char  +  gaz            (aucun atome n'est cree)
        n_E(gaz, i)  =  n_E(i) - n_E(char, i)

puis melange des deux gaz au prorata des masses reellement degazees :

        n_E(gaz) = n_E(gaz, liege) + n_E(gaz, resine)

Usage :
    python cork_pyrolysis_data.py
"""

import csv
import os
import re

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}

# ---------------------------------------------------------------------------
# Donnees materiau
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Liege : composition BIOCHIMIQUE et unites de repetition
# ---------------------------------------------------------------------------
# Plutot que de prendre une analyse elementaire toute faite, on reconstruit le
# liege depuis ses constituants et l'unite de repetition de chacun.
#
#   (formule de l'unite, % masse, commentaire)
#
# Choix des unites -- chacun est un representant, pas une verite :
#   suberine  : polyester glycerol / acides gras en C18. Unite = glycerol
#               esterifie par 3 acides 9,10-epoxy-18-hydroxyoctadecanoiques
#               (le monomere dominant du liege de Quercus suber), moins 3 H2O
#               de condensation :  C3H8O3 + 3 C18H34O4 - 3 H2O = C57H104O12
#   lignine   : lignine de liege essentiellement guaiacyle -> alcool
#               coniferylique C10H12O3
#   polysacch.: unite anhydroglucose C6H10O5 (cellulose et hemicelluloses)
#   tanins    : tanins condenses = proanthocyanidines -> unite flavan-3-ol
#               (catechine) C15H14O6
#   ceroides  : cires du liege, dominees par la friedeline C30H50O
CORK_BIOCHEM = {
    "suberine":     ({"C": 57, "H": 104, "O": 12}, 45.0),
    "lignine":      ({"C": 10, "H":  12, "O":  3}, 27.0),
    "polysacch.":   ({"C":  6, "H":  10, "O":  5}, 12.0),
    "tanins":       ({"C": 15, "H":  14, "O":  6},  6.0),
    "ceroides":     ({"C": 30, "H":  50, "O":  1},  6.0),
}
# Les parts somment a 96 % : le complement (cendres, humidite, extractibles
# mineurs) n'est pas de la matiere C/H/O identifiee. On renormalise sur les
# 96 % declares -- voir la discussion du programme.

# Variante d'unite pour la suberine : acide 9,10-epoxy-18-hydroxyoctadecanoique
# esterifie sans glycerol (- H2O). Plus reduit, donc plus riche en C.
SUBERIN_ALT = {"C": 18, "H": 32, "O": 3}

# Analyse elementaire du liege, valeurs de litterature (Quercus suber, sec
# sans cendres). C'est la donnee RETENUE ; la reconstruction depuis les unites
# de repetition (cork_elemental) sert de controle.
#
# ATTENTION -- PROVENANCE A CONFIRMER. Ces trois nombres sont des valeurs
# usuelles citees pour le liege ; ils n'ont PAS ete verifies sur une source
# primaire accessible depuis cet environnement. A confronter a une analyse
# elementaire publiee (ou mesuree) avant tout usage engageant. La
# reconstruction biochimique, elle, est tracable ligne a ligne.
CORK_MASS_PCT_LITERATURE = {"C": 62.4, "H": 8.5, "O": 28.4}

# Rendements en char des constituants du liege (ordre de grandeur, ATG lente
# sous inerte) : sert de CONTROLE au rendement deduit de la TGA du composite.
CORK_CONSTITUENT_CHAR_YIELD = {
    "suberine":   0.15,   # polyester aliphatique, tres volatil
    "lignine":    0.45,   # aromatique, le plus charbonnant
    "polysacch.": 0.15,   # cellulose ~10-20 %
    "tanins":     0.40,   # polyphenols
    "ceroides":   0.02,   # cires : partent entierement
}

# Resine phenolique : motif novolac C7H6O (phenol + formaldehyde),
# meme resine que celle retenue pour le TACOT.
RESIN_FORMULA = "C7H6O"

# Rendement en char du liege : DEDUIT du rendement COMPOSITE mesure par TGA
# sur le P50 (20 % masse) et du rendement de la resine (50 %) :
#       0.80 * y_liege + 0.20 * 0.50 = 0.20   =>   y_liege = 12.5 %
CORK_CHAR_YIELD = 0.125
RESIN_CHAR_YIELD = 0.50   # donnee de l'enonce (valeur classique novolac)

# Fractions massiques du solide vierge
W_CORK = 0.80
W_RESIN = 0.20

# Char suppose purement carbone pour les deux constituants (cas de base).
CHAR_COMP = {"C": 1.0}

# Masses volumiques (kg/m3) pour la reponse materiau -- n'entrent PAS dans le
# XML, seulement dans la recession s_dot = B'c mdot_e / rho_c.
# Mesures publiees sur le cork P50 (Sakraker et al., CEAS Space J. 14:377-393,
# 2022) : rho_vierge 464.5 et 466.7, rho_char 279.9 et 298.4.
RHO_VIRGIN = 465.6
RHO_CHAR_MEASURED = 289.1

# Rendement en char COMPOSITE mesure par TGA sur le P50 (argon, 10 K/min) :
# 24.5 % a 780 K puis plateau a 20 % jusqu'a 1650 K.
P50_CHAR_YIELD = 0.20


# ---------------------------------------------------------------------------
# Fermeture elementaire
# ---------------------------------------------------------------------------

def parse_formula(formula):
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula):
        if not el:
            continue
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out


def molar_mass(comp):
    return sum(n * M[e] for e, n in comp.items())


def normalize(comp):
    tot = sum(comp.values())
    return {e: n / tot for e, n in comp.items()}


def mass_fractions(comp):
    m = {e: n * M[e] for e, n in comp.items()}
    tot = sum(m.values())
    return {e: v / tot for e, v in m.items()}


def cork_elemental(biochem=None):
    """Analyse elementaire (% masse) du liege depuis ses unites de repetition.

    Chaque constituant apporte la composition elementaire de son unite de
    repetition, ponderee par sa part massique ; les parts sont renormalisees.
    """
    biochem = biochem or CORK_BIOCHEM
    total = sum(w for _, w in biochem.values())
    out = {"C": 0.0, "H": 0.0, "O": 0.0}
    for unit, w in biochem.values():
        y = mass_fractions(unit)
        for e in out:
            out[e] += (w / total) * y.get(e, 0.0)
    return {e: 100.0 * v for e, v in out.items()}


def cork_char_yield_from_constituents(biochem=None, yields=None):
    """Rendement en char du liege par additivite des constituants."""
    biochem = biochem or CORK_BIOCHEM
    yields = yields or CORK_CONSTITUENT_CHAR_YIELD
    total = sum(w for _, w in biochem.values())
    return sum((w / total) * yields[name] for name, (_, w) in biochem.items())


# Analyse elementaire retenue : celle de la LITTERATURE (choix utilisateur).
# Mettre CORK_MASS_PCT = cork_elemental() pour repasser a la reconstruction.
CORK_MASS_PCT = dict(CORK_MASS_PCT_LITERATURE)


def moles_from_mass_pct(mass_pct, mass):
    """Analyse elementaire (% masse) + masse [g] -> moles d'atomes."""
    y = normalize(mass_pct)
    return {e: mass * y[e] / M[e] for e in y}


def moles_from_formula(formula, mass):
    """Motif chimique + masse [g] -> moles d'atomes."""
    unit = parse_formula(formula)
    n_unit = mass / molar_mass(unit)
    return {e: n_unit * v for e, v in unit.items()}


def split_char_gas(constituent_moles, char_yield, char_comp=None):
    """
    constituent_moles : moles d'atomes du constituant vierge
    char_yield        : fraction MASSIQUE restant en char
    char_comp         : composition elementaire du char (defaut C pur)

    Retourne (moles_char, moles_gaz).
    """
    char_comp = normalize(char_comp or CHAR_COMP)
    m_virgin = sum(n * M[e] for e, n in constituent_moles.items())
    m_char = char_yield * m_virgin

    # moles d'atomes du char a partir de sa composition elementaire
    m_per_atom = sum(char_comp[e] * M[e] for e in char_comp)
    n_atoms = m_char / m_per_atom
    char = {e: n_atoms * char_comp[e] for e in char_comp}

    gas = {}
    for e in set(constituent_moles) | set(char):
        v = constituent_moles.get(e, 0.0) - char.get(e, 0.0)
        if v < -1e-9:
            raise ValueError(
                f"fermeture impossible : le char demande plus de {e} que le "
                f"constituant n'en contient (rendement trop eleve ?)")
        gas[e] = max(v, 0.0)
    return char, gas


def fmt(comp, digits=3):
    return ", ".join(f"{e}:{v:.{digits}f}"
                     for e, v in sorted(comp.items()) if v > 1e-12)


# ---------------------------------------------------------------------------
# Bilan complet du composite
# ---------------------------------------------------------------------------

def composite_balance(w_cork=W_CORK, w_resin=W_RESIN,
                      cork_char_yield=CORK_CHAR_YIELD,
                      resin_char_yield=RESIN_CHAR_YIELD,
                      cork_char_comp=None, resin_char_comp=None,
                      basis=100.0):
    """Bilan sur `basis` grammes de composite vierge."""
    m_cork = basis * w_cork / (w_cork + w_resin)
    m_resin = basis * w_resin / (w_cork + w_resin)

    cork = moles_from_mass_pct(CORK_MASS_PCT, m_cork)
    resin = moles_from_formula(RESIN_FORMULA, m_resin)

    cork_char, cork_gas = split_char_gas(cork, cork_char_yield,
                                         cork_char_comp)
    resin_char, resin_gas = split_char_gas(resin, resin_char_yield,
                                           resin_char_comp)

    gas = {e: cork_gas.get(e, 0.0) + resin_gas.get(e, 0.0)
           for e in set(cork_gas) | set(resin_gas)}
    char = {e: cork_char.get(e, 0.0) + resin_char.get(e, 0.0)
            for e in set(cork_char) | set(resin_char)}

    m_char = sum(n * M[e] for e, n in char.items())
    m_gas = sum(n * M[e] for e, n in gas.items())

    return dict(
        m_cork=m_cork, m_resin=m_resin,
        cork_gas=cork_gas, resin_gas=resin_gas,
        gas=gas, char=char,
        m_char=m_char, m_gas=m_gas,
        char_yield=m_char / basis,
        # couplage stationnaire : B'g/B'c = m_gaz/m_char = (rho_v-rho_c)/rho_c
        k=m_gas / m_char,
    )


# ---------------------------------------------------------------------------
# Programme
# ---------------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    line = "=" * 76

    print(line)
    print("0. LE LIEGE RECONSTRUIT DEPUIS SES UNITES DE REPETITION")
    print(line)
    tot = sum(w for _, w in CORK_BIOCHEM.values())
    print(f"  {'constituant':12s} {'% masse':>8}  {'unite':>12}  "
          f"{'M':>8}  {'C':>6} {'H':>6} {'O':>6}  (% masse de l'unite)")
    for name, (unit, w) in CORK_BIOCHEM.items():
        y = mass_fractions(unit)
        f = "".join(f"{e}{int(n)}" for e, n in
                    sorted(unit.items(), key=lambda kv: "CHO".index(kv[0])))
        print(f"  {name:12s} {w:8.1f}  {f:>12}  {molar_mass(unit):8.2f}  "
              f"{100*y['C']:6.1f} {100*y['H']:6.1f} {100*y['O']:6.1f}")
    print(f"  {'':12s} {tot:8.1f}   <- le complement a 100 % (cendres,"
          " humidite, extractibles) n'est pas de la matiere C/H/O identifiee")

    rec = cork_elemental()
    lit = CORK_MASS_PCT_LITERATURE
    print(f"\n  liege reconstruit (parts renormalisees)  : "
          f"C {rec['C']:.2f}  H {rec['H']:.2f}  O {rec['O']:.2f}")
    print(f"  analyse de litterature  <- RETENUE       : "
          f"C {lit['C']:.2f}  H {lit['H']:.2f}  O {lit['O']:.2f}")
    print(f"  ecart                                    : "
          f"C {rec['C']-lit['C']:+.2f}  H {rec['H']-lit['H']:+.2f}  "
          f"O {rec['O']-lit['O']:+.2f}")
    alt = cork_elemental({**CORK_BIOCHEM, "suberine": (SUBERIN_ALT, 45.0)})
    print(f"  variante suberine sans glycerol (C18H32O3) : "
          f"C {alt['C']:.2f}  H {alt['H']:.2f}  O {alt['O']:.2f}")
    print("  -> H est retrouve a 0.2 point pres ; C est surestime de ~4 points")
    print("     (et O sous-estime d'autant). Si les 4 % manquants sont de la")
    print("     matiere oxygenee, la reconstruction rejoint la mesure.")
    bb = composite_balance()
    br = composite_balance()
    import copy
    saved = dict(CORK_MASS_PCT)
    CORK_MASS_PCT.update(rec)
    br = composite_balance()
    CORK_MASS_PCT.clear(); CORK_MASS_PCT.update(saved)
    print(f"  Effet sur le gaz de pyrolyse (a char inchange) :")
    print(f"     avec la litterature (retenue) : {fmt(normalize(bb['gas']))}")
    print(f"     avec la reconstruction        : {fmt(normalize(br['gas']))}")

    print("\n" + line)
    print("1'. CONTROLE DU RENDEMENT EN CHAR PAR ADDITIVITE DES CONSTITUANTS")
    print(line)
    add = cork_char_yield_from_constituents()
    print(f"  {'constituant':12s} {'% masse':>8} {'char':>7}")
    for name, (_, w) in CORK_BIOCHEM.items():
        print(f"  {name:12s} {w:8.1f} {100*CORK_CONSTITUENT_CHAR_YIELD[name]:6.0f}%")
    print(f"  -> liege : {100*add:.1f} %   (contre {100*CORK_CHAR_YIELD:.1f} %"
          f" deduit de la TGA du composite)")
    comp_add = W_CORK * add + W_RESIN * RESIN_CHAR_YIELD
    print(f"  -> composite : {100*comp_add:.1f} %  contre "
          f"{100*P50_CHAR_YIELD:.1f} % mesure sur le P50.")
    print("  Les deux estimations ne se rejoignent pas : le P50 contient un")
    print("  plastifiant glycol qui part entierement en gaz, son rapport")
    print("  liege/resine n'est pas force 80/20, et les rendements par")
    print("  constituant ci-dessus sont des ordres de grandeur. On garde la")
    print("  MESURE (composite 20 %) : c'est elle qui fixe k.")
    bx = composite_balance(cork_char_yield=add)
    print(f"  Pour information, avec un liege a {100*add:.1f} % : "
          f"gaz {fmt(normalize(bx['gas']))}, k = {bx['k']:.3f}")

    print("\n" + line)
    print("LIEGE / PHENOLIQUE  80 / 20 (masse) -- fermeture elementaire")
    print(line)

    b = composite_balance()

    print(f"\n1. Constituants (base 100 g de composite vierge)")
    print(f"   liege  : {b['m_cork']:6.2f} g   "
          f"{fmt(normalize(moles_from_mass_pct(CORK_MASS_PCT, 1.0)))} "
          f"(x molaires)   rendement char {100*CORK_CHAR_YIELD:.1f} %")
    print(f"   resine : {b['m_resin']:6.2f} g   {RESIN_FORMULA} "
          f"(M = {molar_mass(parse_formula(RESIN_FORMULA)):.2f} g/mol)"
          f"        rendement char {100*RESIN_CHAR_YIELD:.0f} %")

    print(f"\n2. Gaz produit par chaque constituant (moles d'atomes)")
    print(f"   liege  -> {fmt(b['cork_gas'], 4)}"
          f"   soit x = {fmt(normalize(b['cork_gas']))}")
    print(f"   resine -> {fmt(b['resin_gas'], 4)}"
          f"   soit x = {fmt(normalize(b['resin_gas']))}")

    print(f"\n3. Gaz de pyrolyse du composite")
    print(f"   moles      : {fmt(b['gas'], 4)}   ({b['m_gas']:.2f} g)")
    print(f"   x molaires : {fmt(normalize(b['gas']))}")
    print(f"   y massiques: {fmt(mass_fractions(b['gas']))}")

    print(f"\n4. Char")
    print(f"   {fmt(normalize(b['char']))}   ({b['m_char']:.2f} g)")

    print(f"\n5. Bilan et couplage stationnaire")
    print(f"   rendement en char du composite : {100*b['char_yield']:.1f} %")
    print(f"   gaz degaze                     : {100*(1-b['char_yield']):.1f} %"
          " de la masse vierge")
    print(f"   k = B'g/B'c = m_gaz/m_char     : {b['k']:.3f}")
    print(f"   rho_v = {RHO_VIRGIN:.0f} kg/m3  ->  rho_c = "
          f"{RHO_VIRGIN*b['char_yield']:.0f} kg/m3 SI le volume est conserve,"
          f"\n   mais la mesure donne {RHO_CHAR_MEASURED:.0f} kg/m3 : le char se"
          f" retracte (V_char/V_v = {RHO_VIRGIN*b['char_yield']/RHO_CHAR_MEASURED:.2f})."
          f"\n   k se calcule donc sur les MASSES, pas sur les densites.")

    print(f"\n6. A recopier dans data/mixtures/cork-air.xml")
    print(f'   <composition name="cork_pyro">{fmt(normalize(b["gas"]))}'
          "</composition>")
    print(f'   <composition name="cork_char">{fmt(normalize(b["char"]), 1)}'
          "</composition>")

    # -----------------------------------------------------------------------
    # Sensibilite au rendement en char du liege (le parametre choisi)
    # -----------------------------------------------------------------------
    print("\n" + line)
    print("SENSIBILITE AU RENDEMENT EN CHAR DU LIEGE")
    print(line)
    print(f"  {'char liege':>10} | {'C':>7} {'H':>7} {'O':>7} | "
          f"{'char comp.':>10} | {'k':>6}")
    rows = []
    for cy in (0.05, 0.10, 0.125, 0.15, 0.20, 0.25):
        bb = composite_balance(cork_char_yield=cy)
        x = normalize(bb["gas"])
        print(f"  {100*cy:>9.1f}% | {x['C']:>7.3f} {x['H']:>7.3f} "
              f"{x['O']:>7.3f} | {100*bb['char_yield']:>9.1f}% | "
              f"{bb['k']:>6.3f}")
        rows.append(dict(cork_char_yield=cy, x_C=x["C"], x_H=x["H"],
                         x_O=x["O"], composite_char_yield=bb["char_yield"],
                         k=bb["k"]))

    # -----------------------------------------------------------------------
    # Sensibilite au rapport liege / resine
    # -----------------------------------------------------------------------
    print("\n" + line)
    print("SENSIBILITE AU RAPPORT LIEGE / RESINE  (a rendements fixes)")
    print(line)
    print(f"  {'liege':>6} | {'C':>7} {'H':>7} {'O':>7} | {'k':>6}"
          "   <- ici le rapport CHANGE la table B'")
    for wc in (0.0, 0.2, 0.5, 0.8, 1.0):
        bb = composite_balance(w_cork=wc, w_resin=1.0 - wc)
        x = normalize(bb["gas"])
        print(f"  {100*wc:>5.0f}% | {x['C']:>7.3f} {x['H']:>7.3f} "
              f"{x['O']:>7.3f} | {bb['k']:>6.3f}")

    # -----------------------------------------------------------------------
    # Variante : char de liege non purement carbone
    # -----------------------------------------------------------------------
    print("\n" + line)
    print("VARIANTE : CHAR DE LIEGE RETENANT DE L'OXYGENE ET DE L'HYDROGENE")
    print(line)
    print("  char liege 90 % C / 2 % H / 8 % O en masse, char resine C pur")
    cork_char_comp = {e: v / M[e] for e, v in
                      (("C", 0.90), ("H", 0.02), ("O", 0.08))}
    bv = composite_balance(cork_char_comp=cork_char_comp)
    print(f"  char     : {fmt(normalize(bv['char']))}"
          f"   -> -char-elem C reste valable (C majoritaire)")
    print(f"  gaz      : {fmt(normalize(bv['gas']))}"
          f"   (contre {fmt(normalize(b['gas']))} avec un char C pur)")
    print(f"  k        : {bv['k']:.3f}")
    print("  -> le char devient MULTI-ELEMENT : `bprime` le gere via")
    print("     -char cork_char -char-elem C, et la table B' change.")

    # -----------------------------------------------------------------------
    # Variante : 80 % en VOLUME et non en masse
    # -----------------------------------------------------------------------
    print("\n" + line)
    print("VARIANTE : 80 % EN VOLUME (liege 120 kg/m3, resine 1200 kg/m3)")
    print(line)
    m_c, m_r = 0.80 * 120.0, 0.20 * 1200.0
    w = m_c / (m_c + m_r)
    bvol = composite_balance(w_cork=w, w_resin=1.0 - w)
    print(f"  -> {100*w:.1f} % de liege EN MASSE seulement")
    print(f"  gaz : {fmt(normalize(bvol['gas']))}   "
          f"char {100*bvol['char_yield']:.1f} %   k = {bvol['k']:.3f}")
    print("  -> l'enonce '80 % de liege' doit imperativement preciser"
          " masse ou volume.")

    # -----------------------------------------------------------------------
    # Calage sur la TGA publiee du cork P50
    # -----------------------------------------------------------------------
    print("\n" + line)
    print("ORIGINE DE L'HYPOTHESE : TGA DU CORK P50 (Sakraker et al. 2022)")
    print(line)
    print("  TGA sous argon, 10 K/min : masse residuelle 24.5 % a 780 K, puis")
    print("  plateau a 20 % jusqu'a 1650 K  ->  rendement en char COMPOSITE 20 %")
    print("  rho_vierge 464.5 / 466.7 kg/m3, rho_char 279.9 / 298.4 kg/m3")
    print("  C'est le cas de base ci-dessus. Variante : le plateau a 780 K.")
    for target in (P50_CHAR_YIELD, 0.245):
        y_cork = (target - W_RESIN * RESIN_CHAR_YIELD) / W_CORK
        bp = composite_balance(cork_char_yield=y_cork)
        print(f"\n  composite {100*target:.1f} %  =>  liege {100*y_cork:.1f} % "
              f"(a 80/20 masse et resine 50 %)")
        print(f"    gaz x molaires : {fmt(normalize(bp['gas']))}")
        print(f"    k = m_gaz/m_char : {bp['k']:.3f}")

    rho_v, rho_c = 465.6, 289.1      # moyennes des deux echantillons
    print(f"\n  ATTENTION : k = (rho_v - rho_c)/rho_c = "
          f"{(rho_v - rho_c)/rho_c:.2f} donnerait un tout autre chiffre.")
    print(f"  Cette identite suppose un volume constant, ce que le P50 ne")
    print(f"  respecte pas : rho_c/rho_v = {rho_c/rho_v:.2f} alors que la masse")
    print(f"  residuelle est 0.20, soit V_char/V_vierge = "
          f"{P50_CHAR_YIELD*rho_v/rho_c:.2f}.")
    print(f"  Le couplage stationnaire se calcule sur les MASSES : k = "
          f"(1-y)/y = {(1-P50_CHAR_YIELD)/P50_CHAR_YIELD:.1f}.")

    out = os.path.join(here, "cork_pyrolysis_data.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n  -> {out}")


if __name__ == "__main__":
    main()
