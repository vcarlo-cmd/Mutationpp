#!/usr/bin/env python
"""
Traçabilité de la résine SC-1008 : des articles publiés au XML de Mutation++.

Rejoue, à partir des seules données publiées, la chaîne qui mène des mesures
de pyrolyse aux trois nombres inscrits dans `data/mixtures/sc1008-*.xml` :

    <composition name="sc1008_pyro">C:0.2526, H:0.6407, O:0.1068</composition>

Références primaires :

  [TS] K. A. Trick, T. E. Saliba, "Mechanisms of the pyrolysis of phenolic
       resin in a carbon/phenolic composite", Carbon 33 (1995) 1509-1515.
       - table 2 : % molaire DANS chaque région réactionnelle
       - table 4 : distribution de chaque espèce SUR les trois régions
       ATTENTION : les données d'évolution gazeuse ne sont PAS une mesure des
       auteurs. Elles proviennent de la littérature (réfs [2,3] de l'article),
       sur une résine phénolique NUE, et sont déconvoluées par eux. Leur
       travail expérimental propre est la FTIR sur préimprégné Amoco T300.

  [Wo] H.-W. Wong, J. Peck, J. Assif, F. Panerai, J. Lachaud, N. N. Mansour,
       "Detailed analysis of species production from the pyrolysis of the
       Phenolic Impregnated Carbon Ablator", J. Anal. Appl. Pyrolysis 122
       (2016) 258-267.
       - table 1 : rendements molaires, mmol/100 mg PICA, 323 -> 1252 K
       - PICA = FiberForm + SC-1008 (Hexion) : c'est LA source SC-1008

Les cinq vérifications :

  1. cohérence interne de [TS]      (tables 2 et 4 se recoupent-elles ?)
  2. [TS] -> composition élémentaire, comparée à la fermeture élémentaire
  3. [Wo] -> composition élémentaire, et diagnostic des écarts
  4. invariant H/O                  (H/O du gaz = H/O de la résine)
  5. convergence du rendement en char sur trois sources indépendantes

Voir `resine_sc1008.md` pour le commentaire détaillé.

Usage :
    python sc1008_bprime/resine_sc1008_verification.py
"""

import os
import re

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007}

# ---------------------------------------------------------------------------
# Données publiées — saisies telles quelles
# ---------------------------------------------------------------------------

# [TS] table 2 : pourcentage MOLAIRE du total de la région
TS_TABLE_2 = {
    1: {"H2O": 49.8, "LMS": 50.1, "CO2": 0.1},
    2: {"H2": 59.4, "CH4": 14.9, "CO": 12.7, "H2O": 12.7, "CO2": 0.2,
        "C2H6": 0.1},
    3: {"H2": 85.7, "CO": 9.5, "H2O": 4.7, "CO2": 0.1},
}

# [TS] table 4 : distribution de chaque espèce sur les régions (% du total)
TS_TABLE_4 = {
    "H2":   {2: 57.1, 3: 42.9},
    "H2O":  {1: 66.3, 2: 28.2, 3: 5.5},
    "LMS":  {1: 100.0},
    "CO":   {2: 71.9, 3: 28.1},
    "CH4":  {2: 100.0},
    "CO2":  {1: 20.5, 2: 66.8, 3: 12.7},
    "C2H6": {2: 100.0},
}

# [Wo] table 1 : totaux intégrés, mmol / 100 mg PICA
WONG_TOTALS = {
    "H2O":    7.97e-1,
    "H2":     2.36e-1,
    "CH4":    6.60e-2,
    "CO":     1.28e-1,
    "CO2":    1.77e-2,
    "C6H6O":  6.65e-3,   # phénol
    "C7H8O":  1.37e-2,   # crésol
    "C8H10O": 8.20e-3,   # diméthylphénol
    "C9H12O": 1.82e-4,   # triméthylphénol
    "C6H6":   3.62e-3,   # benzène
    "C7H8":   4.85e-3,   # toluène
    "C8H10":  1.78e-3,   # xylène
    "C2H4":   8.51e-4,
    "C2H6":   2.38e-4,
    "C3H6":   3.69e-3,
    "C3H8":   3.72e-4,
}

# [Wo] table 1, colonne H2O palier par palier (mmol/100 mg PICA)
WONG_T = [323, 370, 422, 468, 521, 565, 612, 664, 717, 769,
          814, 861, 909, 948, 1009, 1049, 1110, 1148, 1200, 1252]
WONG_H2O = [1.79e-2, 1.37e-2, 4.92e-3, 2.26e-2, 2.53e-2, 2.78e-2, 3.84e-2,
            3.23e-2, 3.20e-2, 5.03e-2, 7.67e-2, 4.08e-2, 3.59e-2, 6.66e-2,
            6.22e-2, 3.92e-2, 5.37e-2, 8.38e-2, 1.95e-2, 5.37e-2]

# [Wo] § 2.2 et 3.2
WONG_PICA_MASS = 100.0      # mg chargés dans le réacteur
WONG_MASS_LOSS = 19.0       # % de perte de masse à 1250 K (balance)
WONG_RESIN_FRACTION = 0.50  # « phenolic resin consists of approximately half
                            #   of the PICA mass »

# Formules brutes (nC, nH, nO)
FORMULA = {
    "H2O": (0, 2, 1), "H2": (0, 2, 0), "CH4": (1, 4, 0), "CO": (1, 0, 1),
    "CO2": (1, 0, 2), "C2H4": (2, 4, 0), "C2H6": (2, 6, 0), "C3H6": (3, 6, 0),
    "C3H8": (3, 8, 0), "C6H6": (6, 6, 0), "C7H8": (7, 8, 0), "C8H10": (8, 10, 0),
    "C6H6O": (6, 6, 1), "C7H8O": (7, 8, 1), "C8H10O": (8, 10, 1),
    "C9H12O": (9, 12, 1),
    # LMS = « low molecular weight substances » = phénol + crésol, non résolus
    # par [TS]. Hypothèse centrale : mélange 50/50.
    "LMS": (6.5, 7.0, 1.0),
}

# Motif de la résine résol cuite retenu (cf. resine_sc1008.md § 2)
MOTIF = (7.5, 6.0, 1.0)   # C7.5 H6 O

# Rendement en char retenu (cf. § 5 : recoupement [TS] / ATG)
Y_RETAINED = 0.55


# ---------------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------------

def elemental(mol_fractions, formula=FORMULA):
    """{espèce: quantité} -> moles atomiques {C, H, O}."""
    n = {"C": 0.0, "H": 0.0, "O": 0.0}
    for sp, q in mol_fractions.items():
        c, h, o = formula[sp]
        n["C"] += c * q
        n["H"] += h * q
        n["O"] += o * q
    return n


def normalize(n):
    t = sum(n.values())
    return {e: v / t for e, v in n.items()}


def mass_of(n):
    return sum(v * M[e] for e, v in n.items())


def gas_from_resin(motif, char_yield):
    """Fermeture élémentaire : gaz = résine - char, char de carbone pur.

    Retourne les moles atomiques du gaz par mole de motif.
    """
    c, h, o = motif
    m_resin = c * M["C"] + h * M["H"] + o * M["O"]
    return {"C": c - char_yield * m_resin / M["C"], "H": h, "O": o}


def fmt(x, ho=None):
    s = f"C:{x['C']:.4f}  H:{x['H']:.4f}  O:{x['O']:.4f}"
    return s + (f"   H/O = {ho:5.2f}" if ho is not None else "")


def rel(a, b):
    return abs(a - b) / b * 100.0


# ---------------------------------------------------------------------------
# 1. Cohérence interne de Trick & Saliba
# ---------------------------------------------------------------------------

def check_ts_consistency():
    """Les tables 2 et 4 se recoupent : elles déterminent N1:N2:N3.

    Pour une espèce j présente dans les régions i et k :
        f_ij * N_i / (f_kj * N_k) = d_ij / d_kj
    donne N_i/N_k. Chaque espèce fournit une estimation indépendante.
    """
    print("=" * 78)
    print("  1. COHÉRENCE INTERNE DE TRICK & SALIBA (tables 2 x 4)")
    print("=" * 78)
    print("""
  Les deux tables ne donnent pas la composition globale, mais elles se
  recoupent : chaque espèce présente dans deux régions fixe le rapport des
  tailles de ces régions.

      f_ij * N_i / (f_kj * N_k) = d_ij / d_kj     ->     N_i / N_k
""")

    def ratio(sp, i, k):
        return (TS_TABLE_4[sp][i] / TS_TABLE_4[sp][k]) * \
               (TS_TABLE_2[k][sp] / TS_TABLE_2[i][sp])

    r23 = {sp: ratio(sp, 2, 3) for sp in ("H2O", "H2", "CO", "CO2")}
    print("  N2/N3, une estimation par espèce commune aux régions 2 et 3 :")
    for sp, v in r23.items():
        flag = "   <- écarté (voir ci-dessous)" if sp == "CO2" else ""
        print(f"      via {sp:4s} : {v:.4f}{flag}")

    kept = [v for sp, v in r23.items() if sp != "CO2"]
    spread = (max(kept) - min(kept)) / (sum(kept) / len(kept)) * 100
    print(f"""
  CO2 pèse 0.1-0.2 % dans la table 2 : l'arrondi à la première décimale y
  domine le signal, son estimation ({r23['CO2']:.2f}) n'a pas de valeur.
  Les trois autres espèces s'accordent à {spread:.1f} % — c'est le test :
  deux tableaux construits indépendamment redonnent le même découpage.
""")

    n23 = sum(kept) / len(kept)
    n12 = ratio("H2O", 1, 2)
    print(f"  N1/N2 (seule H2O est commune aux régions 1 et 2) : {n12:.4f}")
    print(f"\n  RETENU :  N1 : N2 : N3  =  {n12 * n23:.3f} : {n23:.3f} : 1.000")

    ok = spread < 2.0
    print(f"\n  [{'OK' if ok else 'ÉCHEC'}] dispersion {spread:.1f} % < 2 %")
    return {3: 1.0, 2: n23, 1: n12 * n23}, ok


# ---------------------------------------------------------------------------
# 2. Trick & Saliba -> élémentaire, comparé à la fermeture
# ---------------------------------------------------------------------------

def check_ts_vs_closure(regions):
    print()
    print("=" * 78)
    print("  2. TRICK & SALIBA  vs  FERMETURE ÉLÉMENTAIRE")
    print("=" * 78)

    total = {}
    for i, comp in TS_TABLE_2.items():
        for sp, pc in comp.items():
            total[sp] = total.get(sp, 0.0) + pc / 100.0 * regions[i]
    s = sum(total.values())
    print("\n  Composition moléculaire globale reconstruite :")
    for sp, v in sorted(total.items(), key=lambda kv: -kv[1]):
        print(f"      {sp:5s} {v / s:.4f}")

    x_ts = normalize(elemental(total))
    ho_ts = x_ts["H"] / x_ts["O"]
    print(f"\n  -> élémentaire (LMS = C6.5H7O, phénol/crésol 50/50) :")
    print(f"     {fmt(x_ts, ho_ts)}")

    print("\n  Sensibilité à l'hypothèse LMS — c'est la seule vraiment libre :")
    for label, lms in (("100 % phénol C6H6O", (6, 6, 1)),
                       ("50/50        C6.5H7O", (6.5, 7, 1)),
                       ("100 % crésol C7H8O", (7, 8, 1))):
        f2 = dict(FORMULA, LMS=lms)
        x = normalize(elemental(total, f2))
        print(f"      {label:22s} {fmt(x, x['H'] / x['O'])}")

    print("""
  Le motif de résine et le rendement en char n'entrent PAS dans ce calcul :
  il ne sort que des deux tableaux publiés. La fermeture ci-dessous est donc
  une prédiction indépendante.
""")
    best_y, best_d = None, 1e9
    for i in range(300, 750):
        y = i / 1000.0
        x = normalize(gas_from_resin(MOTIF, y))
        d = rel(x["C"], x_ts["C"])
        if d < best_d:
            best_d, best_y = d, y

    print("  Fermeture élémentaire, motif C7.5H6O, char de carbone pur :")
    for y in (0.50, 0.55, 0.575, 0.60):
        x = normalize(gas_from_resin(MOTIF, y))
        mark = ""
        print(f"      Y = {y:.3f}   {fmt(x, x['H'] / x['O'])}"
              f"   écart/[TS] sur C : {rel(x['C'], x_ts['C']):5.1f} %{mark}")

    x_ref = normalize(gas_from_resin(MOTIF, 0.55))
    print(f"\n  Comparaison au point retenu Y = 0.55 :")
    for e in ("C", "H", "O"):
        print(f"      {e} : [TS] {x_ts[e]:.4f}   fermeture {x_ref[e]:.4f}"
              f"   écart {rel(x_ref[e], x_ts[e]):4.1f} %")
    print(f"      H/O : [TS] {ho_ts:.3f}   fermeture "
          f"{x_ref['H'] / x_ref['O']:.3f}   écart "
          f"{rel(x_ref['H'] / x_ref['O'], ho_ts):.1f} %")

    print(f"\n  Rendement en char qui reproduit au mieux [TS] : Y = {best_y:.3f}")

    ok = rel(x_ref["C"], x_ts["C"]) < 5.0
    print(f"\n  [{'OK' if ok else 'ÉCHEC'}] accord fermeture / [TS] sous 5 % sur C")
    return x_ts, best_y, ok


# ---------------------------------------------------------------------------
# 3. Wong et al. -> élémentaire, et diagnostic
# ---------------------------------------------------------------------------

def check_wong():
    print()
    print("=" * 78)
    print("  3. WONG et al. — PICA (SC-1008 / FiberForm), mesure GC")
    print("=" * 78)

    n_raw = elemental(WONG_TOTALS)
    m_raw = mass_of(n_raw)
    x_raw = normalize(n_raw)
    n_mol = sum(WONG_TOTALS.values())

    print(f"""
  Table 1 intégrée sur 323-1252 K : {n_mol:.4f} mmol de gaz / 100 mg de PICA,
  soit {m_raw:.2f} mg — masse molaire moyenne {m_raw / n_mol:.2f} g/mol.
""")
    print(f"  BRUT                        {fmt(x_raw, x_raw['H'] / x_raw['O'])}")

    # -- correction A : eau antérieure à toute décomposition
    d_a = sum(w for t, w in zip(WONG_T, WONG_H2O) if t <= 664)
    n_a = dict(n_raw); n_a["H"] -= 2 * d_a; n_a["O"] -= d_a
    x_a = normalize(n_a)
    print(f"  A  - eau T<=664 K ({d_a:.3f} mmol)  "
          f"{fmt(x_a, x_a['H'] / x_a['O'])}   {mass_of(n_a):5.2f} mg")

    # -- correction B : + eau >850 K, que [Wo] attribue à la désorption
    d_b = d_a + sum(w for t, w in zip(WONG_T, WONG_H2O) if t >= 861)
    n_b = dict(n_raw); n_b["H"] -= 2 * d_b; n_b["O"] -= d_b
    x_b = normalize(n_b)
    print(f"  B  + eau >850 K   ({d_b:.3f} mmol)  "
          f"{fmt(x_b, x_b['H'] / x_b['O'])}   {mass_of(n_b):5.2f} mg")

    # -- correction C : calage sur la balance
    d_c = (m_raw - WONG_MASS_LOSS) / (2 * M["H"] + M["O"])
    n_c = dict(n_raw); n_c["H"] -= 2 * d_c; n_c["O"] -= d_c
    x_c = normalize(n_c)
    print(f"  C  calage balance ({d_c:.3f} mmol)  "
          f"{fmt(x_c, x_c['H'] / x_c['O'])}   {mass_of(n_c):5.2f} mg")

    print(f"""
  Trois causes d'écart, toutes documentées par les auteurs :

  (a) HUMIDITÉ. « The highly porous architecture of the nano-dispersed
      phenolic is responsible for adsorption of atmospheric moisture » ;
      « unavoidable atmospheric moisture leaking into the system at room
      temperature between runs ». L'eau au-delà de 850 K est explicitement
      attribuée à de la désorption. Le GC totalise {m_raw:.1f} mg contre
      {WONG_MASS_LOSS:.0f} mg à la balance : l'excédent de {m_raw - WONG_MASS_LOSS:.1f} mg
      vaut {d_c:.3f} mmol d'eau.
""")

    # -- bilan carbone
    m_resin = WONG_PICA_MASS * WONG_RESIN_FRACTION
    m_char = m_resin - WONG_MASS_LOSS
    y_pica = m_char / m_resin
    c, h, o = MOTIF
    n_motif = m_resin / (c * M["C"] + h * M["H"] + o * M["O"])
    nc_expected = c * n_motif - m_char / M["C"]
    nc_gc = n_raw["C"]
    heavy = ("C6H6O", "C7H8O", "C8H10O", "C9H12O", "C6H6", "C7H8", "C8H10")
    nc_heavy = sum(FORMULA[k][0] * WONG_TOTALS[k] for k in heavy)
    n_heavy = sum(WONG_TOTALS[k] for k in heavy)

    print(f"""  (b) CARBONE MANQUANT. Le GC est limité à M < 400 g/mol.
      résine = {m_resin:.0f} mg, char = {m_char:.0f} mg  ->  Y(PICA) = {y_pica:.3f}
      C attendu par fermeture : {nc_expected:.3f} mmol
      C mesuré au GC          : {nc_gc:.3f} mmol  ({100 * nc_gc / nc_expected:.0f} %)
      -> {100 * (1 - nc_gc / nc_expected):.0f} % du carbone échappe à la mesure.
      Cohérent avec la structure du gaz : les aromatiques ne pèsent que
      {100 * n_heavy / n_mol:.1f} % des moles mais {100 * nc_heavy / nc_gc:.0f} % du carbone mesuré — c'est la
      queue lourde qui manque, et elle est presque tout le carbone.

  (c) PICA n'est pas la résine nue : la préforme capte des radicaux ([Wo]
      mesure H2 divisé par ~5 et des hydrocarbures légers abaissés par
      rapport à la résole nue), ce qui déplace aussi la spéciation.

  Les corrections vont dans le bon sens (H/O {x_raw['H'] / x_raw['O']:.2f} -> {x_b['H'] / x_b['O']:.2f}) mais ne suffisent
  pas : le carbone manquant explique le reste. La mesure de [Wo] est donc
  utilisable pour la SPÉCIATION et pour le rendement en char, PAS comme
  composition élémentaire du gaz.
""")
    ok = x_b["H"] / x_b["O"] > x_raw["H"] / x_raw["O"]
    print(f"  [{'OK' if ok else 'ÉCHEC'}] les corrections d'humidité "
          f"rapprochent bien H/O de la fermeture")
    return y_pica, ok


# ---------------------------------------------------------------------------
# 4. L'invariant H/O
# ---------------------------------------------------------------------------

def read_mixture_compositions(root):
    """Extrait {fichier: {nom: {élément: valeur}}} de data/mixtures/*.xml."""
    out = {}
    d = os.path.join(root, "data", "mixtures")
    if not os.path.isdir(d):
        return out
    pat = re.compile(r'<composition\s+name="([^"]+)"[^>]*>([^<]*)</composition>')
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".xml"):
            continue
        with open(os.path.join(d, fn)) as f:
            for name, body in pat.findall(f.read()):
                comp = {}
                for tok in body.split(","):
                    if ":" not in tok:
                        continue
                    el, val = tok.split(":", 1)
                    comp[el.strip()] = float(val)
                out.setdefault(fn, {})[name] = comp
    return out


def check_ho_invariant(root, x_ts):
    print()
    print("=" * 78)
    print("  4. L'INVARIANT H/O")
    print("=" * 78)
    print("""
  Avec un char de CARBONE PUR, H et O traversent la pyrolyse sans être
  touchés : seul le carbone est prélevé. Donc

      H/O (gaz)  =  H/O (résine)        indépendamment du rendement en char

  C'est un test de cohérence gratuit : il ne dépend d'AUCUNE des hypothèses
  fragiles (motif exact, rendement en char, degré de cuisson).
""")
    print("  Démonstration numérique, motif C7.5H6O :")
    for y in (0.30, 0.50, 0.575, 0.70):
        x = normalize(gas_from_resin(MOTIF, y))
        print(f"      Y = {y:.3f}   {fmt(x, x['H'] / x['O'])}")

    print("""
  Bornes des motifs phénoliques concevables — aucun ne descend sous 3.6 :""")
    for label, (c, h, o) in (
            ("résol prépolymère non cuit C7.5H9O2.5", (7.5, 9, 2.5)),
            ("résol mi-cuit             C7.5H6.68O1.34", (7.5, 6.68, 1.34)),
            ("résol cuit                C7.5H6O", (7.5, 6, 1)),
            ("novolac linéaire          C7H6O", (7, 6, 1)),
            ("TACOT mesuré (Sykes)      C6.32H6.10O", (6.32, 6.10, 1))):
        print(f"      {label:42s} H/O = {h / o:5.2f}")
    print("""
  Et l'ajout d'eau ne sauve rien : C7.5H6O.(H2O)_w donne H/O = (6+2w)/(1+w),
  qui tend vers 2.0 PAR LE HAUT. La valeur 2.0 est une borne asymptotique
  inatteignable.""")
    for w in (0, 1, 2, 5, 20, 1000):
        print(f"      w = {w:5d}   H/O = {(6 + 2 * w) / (1 + w):.3f}")

    print("\n  Application aux mélanges du dépôt :\n")
    print(f"      {'fichier':26s} {'composition':22s} {'H/O':>6s}   verdict")
    print("      " + "-" * 70)
    ok = True
    mixes = read_mixture_compositions(root)
    for fn, comps in mixes.items():
        for name, comp in comps.items():
            if "H" not in comp or "O" not in comp:
                continue
            ho = comp["H"] / comp["O"]
            # Le test ne s'applique qu'aux matrices PHÉNOLIQUES. Le liège
            # (subérine, lignine, polysaccharides) est nativement bien plus
            # oxygéné : un H/O bas y est attendu, pas suspect.
            phenolic = "cork" not in fn
            if not phenolic:
                verdict = "hors test (liège, pas phénolique)"
            elif 5.0 <= ho <= 7.5:
                verdict = "cohérent phénolique"
            elif ho < 3.6:
                verdict = "*** INCOHÉRENT ***"
                ok = False
            else:
                verdict = "à examiner"
            print(f"      {fn:26s} {name:22s} {ho:6.2f}   {verdict}")

    print(f"""
      {'[TS] reconstruit':26s} {'(mesure, résine nue)':22s} """
          f"""{x_ts['H'] / x_ts['O']:6.2f}   référence expérimentale""")

    print(f"\n  [{'OK' if ok else 'ÉCHEC'}] aucun jeu phénolique sous la borne 3.6")
    return ok


# ---------------------------------------------------------------------------
# 5. Convergence du rendement en char
# ---------------------------------------------------------------------------

def check_char_yield(y_ts, y_pica):
    print()
    print("=" * 78)
    print("  5. RENDEMENT EN CHAR — TROIS SOURCES INDÉPENDANTES")
    print("=" * 78)
    rows = [
        ("[TS], via la composition du gaz reconstruite", y_ts),
        ("ATG SC-1008 nue (55-60 % > 650 °C, fiche)", 0.575),
        ("[Wo], PICA : 19 % de perte / 50 % de résine", y_pica),
    ]
    print()
    for label, y in rows:
        print(f"      {label:46s} Y = {y:.3f}")
    lo = min(y for _, y in rows)
    hi = max(y for _, y in rows)
    x_ret = normalize(gas_from_resin(MOTIF, Y_RETAINED))
    print(f"""
  Enveloppe : {lo:.3f} - {hi:.3f}. Les trois voies — composition du gaz,
  thermogravimétrie de la résine nue, bilan de masse du composite — sont
  méthodologiquement disjointes.

  RETENU pour le XML : Y = {Y_RETAINED:.3f}.

  C'est la seule valeur qui satisfasse les deux contraintes fortes à la
  fois : elle tombe dans la fourchette ATG (0.55-0.60) ET elle reproduit la
  composition de [TS] à 1.5 % (§ 2). Le milieu arithmétique 0.575 s'écarte
  de [TS] de 8.9 % : la moyenne des trois sources n'est pas le bon
  estimateur, parce qu'elles ne mesurent pas la même chose.

  La borne haute {hi:.3f} est propre au PICA — la préforme carbone capte des
  radicaux et relève le rendement — et ne vaut pas pour la résine seule.

  Composition portée dans le XML :
      {fmt(x_ret, x_ret['H'] / x_ret['O'])}
""")
    ok = (hi - lo) < 0.12
    print(f"  [{'OK' if ok else 'ÉCHEC'}] dispersion {hi - lo:.3f} < 0.12")
    return ok


# ---------------------------------------------------------------------------

def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print()
    print("  TRAÇABILITÉ DE LA RÉSINE SC-1008")
    print("  Trick & Saliba 1995 + Wong et al. 2016  ->  "
          "data/mixtures/sc1008-*.xml")
    print()
    print("  Résine : Durite SC-1008 (Hexion / Bakelite Synthetics), résol")
    print("  phénol-formaldéhyde monocomposant, réticulant par chauffage seul.")
    print("  Matrice du PICA, sur préforme FiberForm.  ([Wo] § 2.2)")
    print()

    regions, ok1 = check_ts_consistency()
    x_ts, y_ts, ok2 = check_ts_vs_closure(regions)
    y_pica, ok3 = check_wong()
    ok4 = check_ho_invariant(root, x_ts)
    ok5 = check_char_yield(y_ts, y_pica)

    print("\n" + "=" * 78)
    results = [("cohérence interne de Trick & Saliba", ok1),
               ("fermeture élémentaire vs Trick & Saliba", ok2),
               ("diagnostic des écarts de Wong et al.", ok3),
               ("invariant H/O sur les mélanges du dépôt", ok4),
               ("convergence du rendement en char", ok5)]
    for label, ok in results:
        print(f"  [{'OK' if ok else 'ÉCHEC'}] {label}")
    print("=" * 78)
    return 0 if all(ok for _, ok in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
