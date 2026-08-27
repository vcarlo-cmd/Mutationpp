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

# Liege : analyse elementaire (% masse, sec, sans cendres). Valeurs usuelles
# de la litterature sur le liege de Quercus suber (subrine + lignine +
# polysaccharides) : C 62-63 %, H 8.5-9 %, O 28-29 %, N ~0.6 % (neglige ici,
# il serait sinon a ajouter comme element du gaz de pyrolyse).
CORK_MASS_PCT = {"C": 62.4, "H": 8.5, "O": 28.4}

# Resine phenolique : motif novolac C7H6O (phenol + formaldehyde),
# meme resine que celle retenue pour le TACOT.
RESIN_FORMULA = "C7H6O"

CORK_CHAR_YIELD = 0.25    # choix : ATG liege ~ 20-30 % a 1000 K
                          # (le calage sur la TGA du P50 donne plutot 12.5 %,
                          #  cf. la section P50 en fin de programme)
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
    print("LIEGE / PHENOLIQUE  80 / 20 (masse) -- fermeture elementaire")
    print(line)

    b = composite_balance()

    print(f"\n1. Constituants (base 100 g de composite vierge)")
    print(f"   liege  : {b['m_cork']:6.2f} g   "
          f"{fmt(normalize(moles_from_mass_pct(CORK_MASS_PCT, 1.0)))} "
          f"(x molaires)   rendement char {100*CORK_CHAR_YIELD:.0f} %")
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
          f"{RHO_VIRGIN*b['char_yield']:.0f} kg/m3 (sans retrait)")

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
    for cy in (0.15, 0.20, 0.25, 0.30, 0.35, 0.40):
        bb = composite_balance(cork_char_yield=cy)
        x = normalize(bb["gas"])
        print(f"  {100*cy:>9.0f}% | {x['C']:>7.3f} {x['H']:>7.3f} "
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
    print("CALAGE SUR LA TGA DU CORK P50 (Sakraker et al., CEAS Space J. 2022)")
    print(line)
    print("  TGA sous argon, 10 K/min : masse residuelle 24.5 % a 780 K, puis")
    print("  plateau a 20 % jusqu'a 1650 K  ->  rendement en char COMPOSITE 20 %")
    print("  rho_vierge 464.5 / 466.7 kg/m3, rho_char 279.9 / 298.4 kg/m3")
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
