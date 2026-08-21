# La résine du TACOT : choix, type, composition

Document de traçabilité complet. Il répond à trois questions — **quelle résine,
pourquoi celle-là, et quelle composition** — et rattache chaque nombre soit à
une cellule du classeur `TACOT_3.0.xls`, soit à un tableau de la référence
primaire :

> G. F. Sykes, Jr., *Decomposition characteristics of a char-forming phenolic
> polymer used for ablative composites*, NASA TN D-3810, février 1967.

Toutes les vérifications numériques de ce document sont rejouables :

```bash
python tacot_bprime/resine_tacot_verification.py
```

---

## 0. Réponse courte

| Question | Réponse |
|---|---|
| **Type** | Phénol-formaldéhyde **novolac** (stade « B »), réticulée à l'hexaméthylènetétramine |
| **Référence** | **Union Carbide BRP 5549** |
| **Composition mesurée** | C 75.60 / H 6.12 / N 2.35 / O 15.93 % masse → **C₆.₃₂H₆.₁₀O** hors azote |
| **Motif équivalent** | le **phénol C₆H₆O**, à 0.85 point de % près |
| **Pourquoi celle-là** | c'est une résine **industrielle** réellement utilisée sur des modèles arc-jet, pas une résine de laboratoire — et la seule dont on ait le gaz de pyrolyse ET le résidu mesurés sur le même échantillon |

Ce qui finit dans `data/mixtures/tacot-*.xml` n'est pas la résine mais le **gaz**
qu'elle produit : `C:0.206, H:0.679, O:0.115`. Le § 4 montre que ces trois
nombres se reconstruisent **exactement** depuis la table I de Sykes.

---

## 1. Identité de la résine

| | Valeur | Source |
|---|---|---|
| Nature | phénol-formaldéhyde, **novolac** (stade « B ») | Sykes, § MATERIAL, p. 2 |
| Référence commerciale | **Union Carbide BRP 5549** | Sykes, § MATERIAL, p. 2 |
| Durcisseur | hexaméthylènetétramine (HMTA) prémélangée | Sykes, § MATERIAL, p. 2 |
| Cuisson | 165 °C 1 h sous vide (5 torr), montée 210 °C/h | Sykes, p. 3 |
| Post-cuisson | paliers 38 → 94 → 116 → 149 °C, 10 h chacun | Sykes, p. 3 |
| Densité vierge | 1200 kg/m³ | `Description!A12` |

Sykes, § MATERIAL :

> *« The material characterized in this investigation was a phenol-formaldehyde
> resin of the Union Carbide Corporation (type designation, **BRP 5549**). The
> "B" stage or **novolac**-type resin premixed with **hexamethylenetetramine**
> is supplied by the manufacturer. »*

Classeur, `Description!A12` :

> *« Matrix: ex-novolac/formaldehyde polymer cured at 300 K, virgin density
> 1200 kg/m³ »*

### Chimie de la réticulation

Le novolac est formé **en défaut de formaldéhyde** : il reste donc fusible et
soluble tant qu'on ne le durcit pas. La HMTA prémélangée joue le rôle de
donneur : à la cuisson elle se décompose, libère le formaldéhyde qui manquait
pour achever la condensation, et de l'ammoniac qui sert de catalyseur.

Le réseau obtenu (structure p. 3 de Sykes) est un empilement de **cycles
phénoliques pontés par des méthylènes –CH₂–**, chaque cycle portant son
hydroxyle. C'est ce motif qui gouverne toute la suite : le gaz de pyrolyse est
dominé par le phénol lui-même (§ 4) et le résidu est un carbone aromatique.

### Pourquoi cette résine

Ce n'est **pas** un choix de modélisation, et ce n'est pas non plus un choix
des auteurs du TACOT : c'est le matériau que Sykes a caractérisé, et le TACOT
en hérite.

Sykes justifie explicitement son propre choix en introduction (p. 2) : les
études antérieures sur la décomposition des phénoliques avaient été menées sur
des résines *préparées en conditions de laboratoire*, donc

> *« quite different materials from those actually fabricated in large
> quantities for evaluation as reentry ablation materials »*

Il a donc pris la résine de qualité industrielle des modèles arc-jet des
références 1–3, et l'a soumise **au même cycle de cuisson** que dans le
composite, en moulant simplement un disque plus petit (2.5 cm × 10 cm) pour
obtenir un degré de cuisson uniforme.

Conséquence pour le TACOT : c'est la seule phénolique pour laquelle on dispose,
**mesurés sur le même échantillon et le même cycle thermique**, à la fois de la
composition du gaz de pyrolyse (table I) et de l'analyse élémentaire du résidu
(table II). Cette redondance permet le test de fermeture croisée du § 6 — c'est
ce qui fait la valeur de la donnée, bien plus que la nature exacte du novolac.

---

## 2. Composition élémentaire mesurée (Sykes, table II)

Analyse élémentaire du polymère avant chauffage, puis du résidu à chaque
palier, avec la fraction massique résiduelle mesurée en DTA :

| T (°C) | C | H | N | O | m/m₀ |
|---|---|---|---|---|---|
| **ambiant** | **75.60** | **6.12** | **2.35** | **15.93** | **1.000** |
| 150 | 76.08 | 5.81 | 2.02 | 16.09 | 0.994 |
| 200 | 76.98 | 5.82 | 1.72 | 15.48 | 0.979 |
| 300 | 77.14 | 5.60 | 0.78 | 16.48 | 0.966 |
| 400 | 78.45 | 5.42 | 0.42 | 15.71 | 0.930 |
| 450 | 79.40 | 5.36 | 0.40 | 14.84 | 0.902 |
| 500 | 81.53 | 4.89 | — | 13.58 | 0.855 |
| 550 | 84.22 | 4.47 | — | 11.31 | 0.800 |
| 600 | 88.00 | 3.59 | — | 8.41 | 0.730 |
| 650 | 90.81 | 2.85 | — | 6.35 | 0.680 |
| 750 | 92.31 | 1.54 | — | 6.15 | 0.590 |
| **850** | **92.65** | **0.90** | — | **6.45** | **0.540** |

Trois précautions de lecture, toutes signalées par Sykes :

- L'analyse ne dose que **C, H et N** ; **l'oxygène est le complément à 100 %**.
  C'est licite ici : la teneur en cendres de la résine non cuite a été mesurée à
  0.04 %, donc négligeable.
- L'azote (2.35 %) provient de la HMTA et **n'est pas dans le squelette** du
  polymère. Sa fraction résiduelle s'effondre dès les premiers paliers
  (2.35 → 0.42 % à 400 °C, et zéro à 500 °C) : signature d'un composé **piégé**,
  pas d'un azote de structure. C'est ce qui justifie que le gaz du modèle TACOT
  soit déclaré **sans azote** dans les XML.
- De l'eau est elle aussi piégée au moulage : c'est la seule espèce évoluée en
  dessous de 300 °C (fig. 6), et une post-cuisson prolongée à 225 °C pendant
  24 h l'élimine (fig. 8).

### Formule brute

En moles, normalisé sur l'oxygène :

```
C6.32 H6.10 O1   (N0.17)
```

### Comparaison aux motifs idéalisés du réseau novolac

| Motif | C % | H % | O % | écart max |
|---|---|---|---|---|
| **mesuré, renormalisé hors azote** | **77.42** | **6.27** | **16.31** | — |
| **phénol C₆H₆O** | 76.57 | 6.43 | 17.00 | **0.85** |
| novolac linéaire C₇H₆O (1 pont CH₂/cycle) | 79.23 | 5.70 | 15.08 | 1.81 |
| novolac réticulé C₇.₅H₆O (1.5 pont CH₂/cycle) | 80.34 | 5.39 | 14.27 | 2.92 |

Le polymère mesuré est **moins carboné que le réseau idéal**, et tombe très
près du **phénol C₆H₆O** — c'est-à-dire du monomère, pas du réseau ponté.

Sykes explique l'écart (p. 13) : la différence entre l'analyse expérimentale et
l'analyse théorique de la structure supposée est précisément l'argument qui
établit l'existence de composés piégés (eau, dérivés de la HMTA) dans le
polymère moulé. L'eau et l'ammoniac piégés diluent le carbone et enrichissent
artificiellement H et O — d'où le glissement du réseau C₇H₆O vers C₆H₆O.

> Le script `pyrolysis_gas_from_resin.py` était déjà arrivé à cette conclusion
> par **fermeture inverse** (partir du gaz, remonter à la résine). Elle est
> désormais **confirmée par la mesure directe**, et non plus seulement inférée.

---

## 3. Le matériau complet : où se place la résine

| | Valeur | Source |
|---|---|---|
| Fibres | carbone **ex-cellulose**, traitées à 2000 K, 1600 kg/m³, L 1 mm, ⌀ 10 µm | `Description!A11` |
| Matrice | novolac/formaldéhyde, 1200 kg/m³ | `Description!A12` |
| Microstructure | eps_f 0.1 / eps_m 0.1 / **porosité 0.8** | `Description!A15:A17` |
| ρ vierge | 280 kg/m³ | `'Pyrolysis model'!A32` |
| ρ char | 220 kg/m³ | `'Pyrolysis model'!A33` |

La résine ne représente donc que **10 % du volume** du TACOT (et 42.9 % de sa
masse), mais c'est **elle seule** qui produit le gaz de pyrolyse : les fibres,
déjà traitées à 2000 K, ne pyrolysent pas — le classeur le dit explicitement en
`'Pyrolysis model'!H20` : *« The carbon fibers do not pyrolyse »*.

D'où la règle qui structure toute la mise en données (cf. `mise_en_donnees_xml.md`) :

- **gaz de pyrolyse** (`-py`) ← la **résine seule**
- **char** (`-char`) ← fibres + résine carbonisée, toutes deux du carbone pur

---

## 4. Le gaz de pyrolyse : reconstruction exacte

### 4.1 La donnée brute — Sykes, table I

Distribution des produits de décomposition, par chromatographie en phase
gazeuse, paliers de 50 °C entre 100 et 1000 °C, en pourcentage des moles
totales observées :

| T (°C) | CO₂ | CO | C₆H₆ | C₇H₈ | C₆H₅OH | (CH₃)₂C₆H₃OH | CH₄ | H₂O | H₂ |
|---|---|---|---|---|---|---|---|---|---|
| 100 | | | | | | | | 1.47 | |
| 150 | | | | | | | | 0.75 | |
| 200 | | | | | | | | 0.48 | |
| 250 | | | | | | | | 0.57 | |
| 300 | | | | | | | | 1.28 | |
| 350 | | | | | | | | 3.44 | |
| 400 | | | | | 0.46 | 0.13 | 0.05 | 3.44 | |
| 450 | | 0.21 | | 0.08 | 0.81 | 0.25 | 0.15 | 5.42 | 0.76 |
| 500 | 0.09 | 0.44 | 0.02 | 0.13 | 2.72 | 0.75 | 0.75 | 3.35 | 1.47 |
| 550 | 0.11 | 0.87 | 0.06 | 0.05 | 1.62 | 0.38 | 1.29 | 2.44 | 2.18 |
| 600 | 0.32 | 1.30 | 0.06 | 0.04 | 0.79 | 0.14 | 2.61 | 0.40 | 3.65 |
| 650 | 0.51 | 1.19 | 0.03 | 0.01 | 0.44 | 0.10 | 2.35 | 0.26 | 5.17 |
| 700 | 0.26 | 0.77 | | | 0.21 | 0.05 | 1.32 | 0.13 | 5.88 |
| 750 | 0.17 | 0.54 | | | 0.09 | | 0.83 | | 6.64 |
| 800 | 0.11 | 0.26 | | | | | 0.40 | | 7.35 |
| 850 | | 0.20 | | | | | 0.20 | | 5.88 |
| 900 | | | | | | | 0.08 | | 4.50 |
| 950 | | | | | | | | | 3.65 |
| 1000 | | | | | | | | | 2.94 |
| **somme réelle** | **1.57** | **5.78** | **0.17** | **0.31** | **7.14** | **1.80** | **10.03** | **23.43** | **50.07** |
| *total imprimé* | *1.6* | *5.5* | *0.2* | *0.3* | *7.1* | *1.8* | *10.0* | *23.4* | *50.1* |

**Détail décisif** : les totaux imprimés au bas de la table sont arrondis et
somment à 100.0, alors que les colonnes elles-mêmes somment à **100.30**. C'est
cette dernière valeur que le classeur utilise pour renormaliser — c'est ce qui
permet la reproduction exacte ci-dessous.

### 4.2 La réduction TACOT

`'Pyrolysis-gas chemistry'!A8` :

> *« Toluene [C7H8] and Xylenol [C8H10O], present in Sykes but not modeled in
> the finite-rate mechanism of April have been transferred to C6H6 and C6H6O,
> respectively, for consistency of the model. »*

soit :

- C₇H₈ → C₆H₆ : 0.17 + 0.31 = **0.48**
- 2,4-xylénol → C₆H₅OH : 7.14 + 1.80 = **8.94**

puis renormalisation par **100.30**.

> Cette réduction est faite pour rester compatible avec le mécanisme cinétique
> d'April (`'Pyrolysis-gas chemistry'!A11`, 10 réactions, Chemkin), qui ne
> connaît ni le toluène ni le xylénol. Elle est **conservative en atomes à peu
> près, pas exactement** : C₇H₈ → C₆H₆ perd un CH₂, C₈H₁₀O → C₆H₆O en perd
> deux. L'effet reste sous le pourcent vu les quantités en jeu.

### 4.3 Résultat — reproduction exacte

Mole fractions comparées à `'Pyrolysis-gas chemistry'!B5:H5` :

| espèce | reconstruit | classeur TACOT 3.0 | écart |
|---|---|---|---|
| CO₂ | 0.0156530408773679 | 0.0156530408773679 | 7e-18 |
| CO | 0.0576271186440678 | 0.0576271186440678 | 1e-17 |
| C₆H₆ | 0.0047856430707876 | 0.0047856430707876 | 4e-18 |
| C₆H₅OH | 0.0891326021934197 | 0.0891326021934197 | 4e-17 |
| CH₄ | 0.1000000000000000 | 0.1000000000000000 | 3e-17 |
| H₂O | 0.2335992023928215 | 0.2335992023928220 | 5e-16 |
| H₂ | 0.4992023928215353 | 0.4992023928215351 | 3e-16 |

**Écart maximal 5×10⁻¹⁶** — la précision de la double précision. La chaîne est
donc entièrement élucidée.

### 4.4 Du moléculaire à l'élémentaire

Sommation des atomes, `n_E = Σ_j ν_{E,j} · x_j`, puis normalisation :

```
C = 0.7367   H = 2.4290   O = 0.4117  mol   (total 3.5774)
-> C:0.2059   H:0.6790   O:0.1151
```

arrondi en `'Pyrolysis model'!B4:B6` en **C:0.206, H:0.679, O:0.115** — la
valeur portée par les trois XML `data/mixtures/tacot-*.xml`.

Accessoirement : M_gaz ≈ **5.0 g/mol**, et en masse C 49.5 % / H 13.7 % /
O 36.8 %.

`'Pyrolysis model'!A8` porte d'ailleurs l'attribution explicite :

> *« [from table 1 of DECOMPOSITION CHARACTERISTICS OF A CHAR-FORMING PHENOLIC
> POLYMER USED FOR ABLATIVE COMPOSITES; George F. Sykes; NASA TN D 3810;
> 1967.] »*

---

## 5. Le rendement en char de 50 %

### 5.1 Il n'y a aucune cellule « char yield »

La valeur de 50 % est **implicite** dans le classeur. Elle y apparaît deux fois,
numériquement concordantes mais **microstructuralement contradictoires** :

**(a) Bilan de microstructure**, `'Pyrolysis model'!A32:A34` :

```
rho_v = eps_f*rho_f + eps_m*rho_m     = 0.1*1600 + 0.1*1200  = 280 kg/m³
rho_c = eps_f*rho_f + eps_m_c*rho_m_c = 0.1*1600 + 0.05*1200 = 220 kg/m³
hyp: intrinsic densities of virgin and charred matrix are equal
     (only volume loss due to shrinkage during the pyrolysis)
```

→ la matrice passe de `eps_m = 0.10` à `0.05` **à densité intrinsèque
constante** : rendement massique `0.05×1200 / (0.10×1200)` = **0.50**.

**(b) Cinétique de Goldstein**, `'Pyrolysis model'!B18:C19` :

| réaction | ρ_v (kg/m³) | ρ_c (kg/m³) | A (s⁻¹) | E/R (K) | ψ | T_reac (K) |
|---|---|---|---|---|---|---|
| A | 300 | 0 | 1.2e4 | 8 556 | 3 | 333 |
| B | 900 | 600 | 4.48e9 | 20 444 | 3 | 556 |
| C (fibres) | 1600 | 1600 | 0 | 0 | 0 | 5 556 |

→ la matrice passe de `300 + 900 = 1200` à `0 + 600 = 600` kg/m³, **à volume
constant** : rendement **0.50** également.

Les deux lectures décrivent des physiques différentes (perte de volume à densité
constante *vs* perte de densité à volume constant) mais donnent les mêmes
ρ_v = 280 et ρ_c = 220 kg/m³, seules grandeurs qui comptent en pratique :

```
(1-φ)·[γ·ρ_matrice + (1-γ)·ρ_fibre]   avec φ = 0.8, γ = 0.5
  vierge : 0.2·(0.5·1200 + 0.5·1600) = 280
  char   : 0.2·(0.5· 600 + 0.5·1600) = 220
```

### 5.2 Mais Sykes mesure 0.54

Table II, dernière ligne : fraction massique résiduelle **0.540** à 850 °C.

Test de fermeture indépendant — partir de la résine **mesurée** (table II,
ambiant), retrancher un char de carbone pur, comparer au gaz de la **table I** :

| rendement | origine | gaz obtenu | écart max / table I |
|---|---|---|---|
| 0.500 | TACOT : microstructure + Goldstein | C:0.2317 H:0.6600 O:0.1082 | 12.5 % |
| **0.540** | **Sykes table II, mesuré à 850 °C** | **C:0.2028 H:0.6848 O:0.1123** | **2.4 %** |
| 0.534 | ajustement optimal | C:0.2072 H:0.6811 O:0.1117 | 2.9 % |
| — | *cible : gaz de la table I* | *C:0.2059 H:0.6790 O:0.1151* | — |

**Conclusion.** Le 0.54 mesuré ferme le bilan à 2.4 % près entre **deux
techniques expérimentales totalement indépendantes** — analyse élémentaire du
résidu (DTA) d'un côté, chromatographie des effluents de l'autre. C'est une
validation forte de la cohérence interne de Sykes, et accessoirement du fait
que le char soit assimilable à du carbone pur.

Le **0.50 du TACOT est une valeur de modèle**, arrondie pour produire des
densités rondes (280 / 220 kg/m³) et une fraction volumique de matrice qui se
divise proprement par deux. Ce n'est pas la mesure.

### 5.3 Pourquoi c'est sans conséquence sur la table B'

`bprime` ne consomme que **trois compositions élémentaires** (bord de couche
limite, gaz de pyrolyse, char) plus T, P et B'g — cf. `src/apps/bprime.cpp:304-330`.

La composition du gaz est prise **directement de la table I**, elle n'est
jamais reconstruite à partir du rendement en char. Ce dernier n'intervient que
dans la **réponse matériau**, via le couplage

```
k = B'g / B'c = (ρ_v - ρ_c) / ρ_c = (280 - 220) / 220 = 0.2727
```

c'est-à-dire dans le choix du point de fonctionnement sur la table, pas dans la
table elle-même. Cf. `autre_materiau.md` § 4.

### 5.4 Le char n'est pas du carbone pur à 850 °C

Sykes mesure **92.65 % C / 0.90 % H / 6.45 % O** à 850 °C. L'hypothèse `C:1.0`
des XML reste néanmoins légitime :

- en ablation, la surface dépasse largement 2000 K, où la carbonisation est
  complète (les 6.45 % d'oxygène résiduels partent en CO bien avant) ;
- les fibres ex-cellulose sont déjà traitées à 2000 K (`Description!A11`), donc
  du carbone pur, et elles représentent 57 % de la masse vierge.

---

## 6. Récapitulatif des sources

| Donnée | Valeur | Localisation exacte |
|---|---|---|
| Type de résine | novolac phénol-formaldéhyde + HMTA | Sykes p. 2 ; `Description!A12` |
| Référence commerciale | Union Carbide BRP 5549 | Sykes p. 2 |
| Cycle de cuisson | 165 °C, post-cuisson jusqu'à 149 °C | Sykes p. 3 |
| Composition résine vierge | C 75.60 / H 6.12 / N 2.35 / O 15.93 (% masse) | Sykes table II |
| Formule brute équivalente | C₆.₃₂H₆.₁₀O ≈ phénol C₆H₆O | dérivée de la table II |
| Densité résine vierge | 1200 kg/m³ | `Description!A12` |
| Fibres | carbone ex-cellulose, 2000 K, 1600 kg/m³, L 1 mm, ⌀ 10 µm | `Description!A11` |
| Microstructure | eps_f 0.1 / eps_m 0.1 / porosité 0.8 | `Description!A15:A17` |
| ρ vierge / ρ char | 280 / 220 kg/m³ | `'Pyrolysis model'!A32:A33` |
| Gaz, composition moléculaire | 7 espèces | `'Pyrolysis-gas chemistry'!B5:H5` ← Sykes table I |
| Gaz, réduction toluène/xylénol | → C₆H₆ / C₆H₅OH | `'Pyrolysis-gas chemistry'!A8` |
| **Gaz, composition élémentaire** | **C:0.206 H:0.679 O:0.115** | `'Pyrolysis model'!B4:B6` = les XML |
| Cinétique de décomposition | Goldstein 1965, 2 phases A/B | `'Pyrolysis model'!A11:G20` |
| Mécanisme cinétique du gaz | April 1969, 10 réactions | `'Pyrolysis-gas chemistry'!A11:A34` |
| Rendement en char (modèle) | 0.50 | implicite : `'Pyrolysis model'!A32:A33` **et** `B18:C19` |
| Rendement en char (mesuré) | 0.54 à 850 °C | Sykes table II |
| Char, composition à 850 °C | C 92.65 / H 0.90 / O 6.45 | Sykes table II |
| Chaleur de pyrolyse | 293 kJ/kg entre 350 et 850 °C | Sykes, résumé |

---

## 7. Fichiers liés

| Fichier | Rôle |
|---|---|
| `resine_tacot_verification.py` | rejoue les 4 vérifications de ce document |
| `pyrolysis_gas_from_resin.py` | fermeture résine → gaz, sensibilité au rendement en char |
| `resin_ranges_to_pyro.py` | propagation d'incertitude quand la résine est connue par plages |
| `material_response.py` | bilan de phase solide, ligne de fonctionnement B'g/B'c |
| `mise_en_donnees_xml.md` | ce qu'il faut renseigner dans un XML de mélange |
| `../zuram_bprime/resine_zuram.md` | même travail pour le ZURAM (novolac + hexamine) |
| `autre_materiau.md` | ce que change (et ne change pas) le rapport fibres/résine |
| `data/mixtures/tacot-air_25.xml` | 25 espèces de la table B' de référence + C(gr) |
| `data/mixtures/tacot-air_35.xml` | liste étendue |
| `data/mixtures/tacot-pyrogas.xml` | gaz de pyrolyse seul, pour h_g(T, P) |
