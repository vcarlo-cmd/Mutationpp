#!/usr/bin/env python
"""
Composition élémentaire du gaz de pyrolyse à partir de celle de la résine.

Principe — fermeture élémentaire. La pyrolyse ne crée ni ne détruit d'atomes :

        résine  =  char  +  gaz de pyrolyse

donc, élément par élément :

        n_E(gaz)  =  n_E(résine)  -  n_E(char)

C'est exactement ce dont la table B' a besoin : le solveur d'équilibre ne
consomme que des fractions élémentaires, jamais la spéciation moléculaire.

Il faut donc connaître trois choses, et trois seulement :
    1. la composition élémentaire de la résine vierge   (formule brute)
    2. le rendement en char de la résine                (fraction massique)
    3. la composition élémentaire du char               (souvent C pur)

Le rendement en char est le paramètre critique : il se mesure par ATG et
dépend de l'histoire thermique. Ce script quantifie cette sensibilité.

Ce qu'on ne peut PAS obtenir par cette voie : la composition MOLÉCULAIRE du
gaz (quelles espèces, en quelles proportions). Elle demande de l'expérience
(ATG-SM, Py-GC/MS) ou un modèle cinétique détaillé — mais elle est inutile
pour la table B'.

Usage :
    python pyrolysis_gas_from_resin.py
"""

import re

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "Si": 28.086}


# ---------------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------------

def parse_formula(formula):
    """'C6H6O' ou 'C7H6O' -> {'C':6, 'H':6, 'O':1}."""
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula):
        if not el:
            continue
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    if not out:
        raise ValueError(f"formule illisible : {formula!r}")
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


def fmt(comp, digits=4):
    return ", ".join(f"{e}:{v:.{digits}f}"
                     for e, v in sorted(comp.items()) if v > 1e-12)


# ---------------------------------------------------------------------------
# Fermeture directe : résine + char -> gaz
# ---------------------------------------------------------------------------

def pyrolysis_gas(resin, char_yield, char=None):
    """
    resin      : dict {élément: nb d'atomes} de la résine vierge
    char_yield : fraction MASSIQUE de résine restant en char
    char       : dict de la composition élémentaire du char
                 (défaut : carbone pur)

    Retourne (gaz_moles, gaz_frac_molaires, diagnostic).
    Lève ValueError si la fermeture est impossible (char trop riche en un
    élément que la résine ne contient pas en quantité suffisante).
    """
    if char is None:
        char = {"C": 1.0}
    if not 0.0 <= char_yield < 1.0:
        raise ValueError("char_yield doit être dans [0, 1[")

    m_resin = molar_mass(resin)              # g par mole de motif
    m_char = char_yield * m_resin            # g de char par mole de motif

    # moles d'atomes dans le char, à partir de sa composition élémentaire
    char_n = normalize(char)
    m_char_per_atom = sum(char_n[e] * M[e] for e in char_n)
    n_char_atoms = m_char / m_char_per_atom
    char_moles = {e: n_char_atoms * char_n[e] for e in char_n}

    # gaz = résine - char
    gas = {}
    for e in set(resin) | set(char_moles):
        v = resin.get(e, 0.0) - char_moles.get(e, 0.0)
        if v < -1e-9:
            raise ValueError(
                f"fermeture impossible : le char demande plus de {e} que la "
                f"résine n'en contient (rendement en char trop élevé, ou "
                f"composition de char incompatible)")
        gas[e] = max(v, 0.0)

    if sum(gas.values()) <= 0:
        raise ValueError("gaz nul : rendement en char de 100 % ?")

    diag = dict(m_resin=m_resin, m_char=m_char, m_gas=m_resin - m_char,
                char_moles=char_moles)
    return gas, normalize(gas), diag


# ---------------------------------------------------------------------------
# Fermeture inverse : gaz + char -> résine  (validation)
# ---------------------------------------------------------------------------

def resin_from_gas(gas_frac, char_yield, char=None, basis=100.0):
    """Reconstruit la résine à partir du gaz — sert à valider une donnée."""
    if char is None:
        char = {"C": 1.0}
    m_char = char_yield * basis
    m_gas = basis - m_char

    g = normalize(gas_frac)
    n_gas = m_gas / sum(g[e] * M[e] for e in g)
    n = {e: n_gas * g[e] for e in g}

    c = normalize(char)
    n_char = m_char / sum(c[e] * M[e] for e in c)
    for e in c:
        n[e] = n.get(e, 0.0) + n_char * c[e]
    return n


# ---------------------------------------------------------------------------
# Programme
# ---------------------------------------------------------------------------

def show(label, resin_formula, char_yield, char=None):
    resin = parse_formula(resin_formula)
    gas, gasf, d = pyrolysis_gas(resin, char_yield, char)
    print(f"\n  {label}")
    print(f"    résine {resin_formula}  (M = {d['m_resin']:.2f} g/mol), "
          f"rendement char {100*char_yield:.0f} %")
    print(f"    -> gaz  moles      : {fmt(gas)}")
    print(f"    -> gaz  x molaires : {fmt(gasf)}")
    print(f"    -> gaz  y massique : {fmt(mass_fractions(gas))}")
    print(f"    XML : <composition name=\"..._pyro\">{fmt(gasf, 3)}"
          "</composition>")
    return gasf


if __name__ == "__main__":
    print("=" * 76)
    print("1. VALIDATION — reconstruire la résine du TACOT depuis son gaz")
    print("=" * 76)
    tacot_gas = {"C": 0.206, "H": 0.679, "O": 0.115}
    n = resin_from_gas(tacot_gas, 0.50)
    yf = mass_fractions(n)
    print(f"  gaz TACOT (donnée)  : {fmt(tacot_gas, 3)}, rendement char 50 %")
    print(f"  -> résine reconstruite, ratios sur O : "
          f"C{n['C']/n['O']:.2f} H{n['H']/n['O']:.2f} O1")
    print(f"  -> % masse : C {100*yf['C']:.1f}  H {100*yf['H']:.1f}  "
          f"O {100*yf['O']:.1f}")
    print()
    for nom, f in (("phénol   C6H6O", "C6H6O"),
                   ("novolac  C7H6O", "C7H6O")):
        y = mass_fractions(parse_formula(f))
        print(f"     comparaison {nom} : C {100*y['C']:.1f}  "
              f"H {100*y['H']:.1f}  O {100*y['O']:.1f}")
    print("\n  -> la résine du TACOT est très proche du phénol C6H6O.")
    print("     La fermeture est donc cohérente avec la donnée expérimentale.")

    print("\n" + "=" * 76)
    print("2. SENS DIRECT — du motif de résine au gaz de pyrolyse")
    print("=" * 76)
    show("Phénol, char 50 %",  "C6H6O", 0.50)
    show("Novolac, char 50 %", "C7H6O", 0.50)
    show("Phénol, char 60 %",  "C6H6O", 0.60)

    print("\n" + "=" * 76)
    print("3. SENSIBILITÉ AU RENDEMENT EN CHAR  (résine phénol C6H6O)")
    print("=" * 76)
    print(f"  {'char [%]':>9} | {'C':>7} {'H':>7} {'O':>7} | "
          f"{'M_gaz':>7} | commentaire")
    _, ref, _ = pyrolysis_gas(parse_formula("C6H6O"), 0.50)
    for cy in (0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70):
        try:
            gas, gasf, d = pyrolysis_gas(parse_formula("C6H6O"), cy)
        except ValueError as e:
            print(f"  {100*cy:>9.0f} | {'--':>7} {'--':>7} {'--':>7} | "
                  f"{'--':>7} | {e}")
            continue
        Mg = d["m_gas"] / sum(gas.values())
        if cy == 0.50:
            note = "référence"
        else:
            dev = max(abs(gasf[e] - ref[e]) / ref[e] for e in ref) * 100
            note = f"écart max / 50 % : {dev:5.1f} %"
        print(f"  {100*cy:>9.0f} | {gasf['C']:>7.4f} {gasf['H']:>7.4f} "
              f"{gasf['O']:>7.4f} | {Mg:>7.3f} | {note}")

    print("\n  Le rendement en char pilote fortement la composition du gaz :")
    print("  c'est LA donnée à mesurer (ATG), bien avant la spéciation.")
    print("\n  Note : à rendement croissant, le char capte de plus en plus de")
    print("  carbone, donc le gaz s'enrichit relativement en H et en O.")
