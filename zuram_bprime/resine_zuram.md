# La résine du ZURAM : choix, type, composition

Document de traçabilité, pendant de `tacot_bprime/resine_tacot.md`. Il rattache
chaque nombre soit à une cellule du classeur, soit à la publication primaire :

> [ODS] `Tacot_Zuram_Calcarb_database_v4.3.1.ods`, onglet **`ZURAM_official`**
> — A. Turchi, F. Torres-Herrador, J. Rico Orero (VKI), 2021.
>
> [THo] F. Torres-Herrador, A. Eschenbacher, J. Coheur, J. Blondeau,
> T. E. Magin, K. M. Van Geem, *Decomposition of carbon/phenolic composites for
> aerospace heatshields: Detailed speciation of phenolic resin pyrolysis
> products*, **Aerospace Science and Technology 119 (2021) 107079**.

Vérifications rejouables :

```bash
python zuram_bprime/resine_zuram_verification.py
```

---

## 0. Réponse courte

| Question | Réponse |
|---|---|
| **Matériau** | ZURAM® 18/50, ablateur carbone/phénolique du **DLR** |
| **Renfort** | préforme **Mersen Calcarb CBCF 18/2000** (fibres de carbone) |
| **Résine** | phénolique, **catalysée à l'hexamine** (hexaméthylènetétramine) |
| **Type déduit** | **novolac** — l'hexamine est le durcisseur des systèmes novolac deux étages, et la composition mesurée tombe sur le motif C₇H₆O |
| **Composition mesurée** | C 75.2 / H 5.8 / N 1.4 / O 14.1 % masse → **C₇.₀₉H₆.₄₆O** |
| **Motif équivalent** | **novolac linéaire C₇H₆O**, à **0.35** point de % près |

> **Réserve importante.** Contrairement au TACOT — où Sykes nomme la résine
> (*Union Carbide BRP 5549, novolac*) — la formulation exacte du ZURAM **n'est
> pas publique**. [THo] l'écrit noir sur blanc : *« the preparation of the
> material is often not publicly available »*. « Novolac » est donc ici une
> **déduction** (hexamine + composition élémentaire), pas une donnée déclarée.

Ce qui finit dans `data/mixtures/zuram-*.xml` n'est pas la résine mais le
**gaz** : `C:0.171, H:0.722, O:0.107`. Le § 4 montre d'où viennent ces trois
nombres — et la réponse diffère radicalement du TACOT.

---

## 1. Identité du matériau et de la résine

`'ZURAM_official'!A4` :

> *« Carbon-phenolic ablator produced by DLR. Composed by **Mersen Calcarb
> Preform (CBCF 18/2000)** impregnated with phenolic resin. A mix of available
> data and data obtained through the **AblaNTIS** characterization campaign is
> used »*

[THo] § 4.2, sur le catalyseur :

> *« The production of NH₃ and other nitrogenated compounds are due to the
> catalyst used in the production of the material (**hexamine** in the case of
> Zuram®). This is confirmed by a relatively constant mass yield of 1.5 % of
> NH₃ for temperatures above 500 °C. »*

et, plus prudemment, sur la formulation en général :

> *« Even though the preparation of the material is often not publicly
> available, some authors mention that nitrogen containing compounds
> (hexamine, NH₄OH) are often used as catalyst in the synthesis of phenolic
> resins »*

### Pourquoi « novolac »

Deux arguments convergents, aucun décisif seul :

1. **L'hexamine.** C'est le durcisseur du système **novolac deux étages** :
   un novolac formé en défaut de formaldéhyde reste fusible, et l'hexamine
   apporte à la cuisson le formaldéhyde manquant (plus l'ammoniac catalyseur).
   Un résol, lui, est formé en excès de formaldéhyde et réticule seul, sans
   durcisseur. La présence d'hexamine oriente donc vers un novolac.
   [THo] mesure d'ailleurs 1.5 % de NH₃ en masse au-dessus de 500 °C — c'est
   la signature de l'hexamine.
2. **La composition élémentaire** (§ 2) tombe à 0.35 point du motif novolac
   linéaire C₇H₆O.

### Que signifie « 18/50 » ?

Le classeur documente le **18** et reste muet sur le **50**.

**18 — la préforme.** Il vient directement de la désignation Mersen
**Calcarb® CBCF 18/2000**, où 18 est la **masse volumique apparente en
0.01 g/cm³**, soit **180 kg/m³**. Le classeur le confirme deux fois :

- `Calcarb_official!D15` : *Average density* = **180**, et
  `Calcarb_official!D12` donne une fraction volumique de fibres de 0.114141
  qui, multipliée par la densité intrinsèque 1577, redonne exactement 180.
- Note de version 4.3.0 : *« Updated volume fraction in Calcarb_official,
  computed now from fiber intrinsic density and **nominal average density of
  180 kg/m³** »*.

(Le **2000** de `CBCF 18/2000` est, dans la nomenclature Mersen, la
**température de traitement thermique** de la préforme en °C. Ce point n'est
pas énoncé dans les documents dont nous disposons.)

**50 — non documenté.** Aucun des fichiers fournis n'explique ce second nombre.
La lecture qui colle aux chiffres est la **teneur en résine, ~50 % en masse** :

| Lecture | Résine, % masse |
|---|---|
| `ZURAM_official`, microstructure (205 kg/m³ fibres + 214 résine) | **51.07** |
| `ZURAM_official!B76`, échantillon TGA | **54.18** |
| version « V12 » à 0.366 g/cm³, préforme nominale 180 | **50.82** |
| version livrée à 0.43 g/cm³, préforme mesurée 205 | 52.33 |
| version livrée à 0.43 g/cm³, préforme nominale 180 | 58.14 |

L'hypothèse concurrente — 50 = masse volumique visée de 0.50 g/cm³ — est
**exclue** par le courriel du DLR cité en `Versions_and_issues!G25` :

> *« The final density of the provided virgin ZURAM® **(18/50)** is
> **0.43 g/cm³**. However, in the thermal conductivity data files you provided,
> the density indicated (ZURAM® V12) is 0.366 g/cm³. »*

Deux matériaux portant le même nom 18/50 y ont des masses volumiques de 0.43 et
0.366 g/cm³ : le nombre ne peut donc pas désigner la densité finale. Le DLR
explique l'écart par l'amélioration du procédé d'infiltration (passage d'un
dessiccateur à vide manuel à un moule RTM à vide asservi), *« which led to an
increase in density »* — la préforme, elle, *« was all the time CBCF 18-2000 »*.

En résumé : **18 = la préforme Calcarb à 180 kg/m³** (certain), **50 ≈ la part
de résine en masse** (déduction cohérente avec toutes les mesures, mais non
énoncée dans les sources disponibles).

### Comment la résine a été prélevée

[THo] § 3.1 : les essais ont porté sur **la résine seule**, prélevée sur une
plaque de ZURAM à un endroit *« where an excess of resin did not impregnate the
fibers, but it underwent the same preparation processes »*, puis broyée. Des
essais préliminaires avaient montré que les fibres n'affectent pas la
distribution des produits de pyrolyse.

C'est exactement la précaution qu'avait prise Sykes pour le TACOT (même cycle
de cuisson, échantillon de résine pure) — et c'est ce qui rend les deux jeux de
données comparables.

---

## 2. Composition élémentaire mesurée

[THo] fig. 7 compare la composition élémentaire de la résine vierge obtenue par
**analyseur élémentaire** (mesure directe) et par **micropyrolyse** (somme de
tous les produits quantifiés + char).

Relevé graphique des barres, en % masse :

| | C | H | N | O | total |
|---|---|---|---|---|---|
| **analyseur élémentaire** | **75.21** | **5.75** | **1.44** | **14.13** | 96.5 |
| micropyrolyse | 74.17 | 5.27 | 1.04 | 19.32 | **99.8** |

> **Contrôle de la calibration.** [THo] précise que la colonne micropyrolyse a
> été *« re-scaled to 100 % from a mass balance of ~90 % »*. Notre relevé la
> donne à **99.8 %** : la lecture graphique est donc juste à 0.2 point près.

[THo] commente lui-même l'écart : bon accord sur C, H et N, mais discordance
sur l'oxygène, attribuée aux composés de haute masse molaire non détectés (qui
contenaient probablement de l'azote, réduisant d'autant la part attribuée à O).

### Formule brute

Renormalisée sur C+H+N+O = 100 : C 77.91 / H 5.96 / N 1.49 / O 14.64, soit

```
C7.09 H6.46 O1   (N0.12)
```

| Motif | C % | H % | O % | écart max |
|---|---|---|---|---|
| **ZURAM mesuré, hors azote** | **79.09** | **6.05** | **14.86** | — |
| phénol C₆H₆O | 76.57 | 6.43 | 17.00 | 2.52 |
| **novolac linéaire C₇H₆O** | 79.23 | 5.70 | 15.08 | **0.35** |
| novolac réticulé C₇.₅H₆O | 80.34 | 5.39 | 14.27 | 1.24 |

### Comparaison au TACOT

| | C | H | N | O | motif | écart |
|---|---|---|---|---|---|---|
| **ZURAM** ([THo] fig. 7) | 75.21 | 5.75 | 1.44 | 14.13 | **novolac C₇H₆O** | 0.35 |
| **TACOT** (Sykes, table II) | 75.60 | 6.12 | 2.35 | 15.93 | **phénol C₆H₆O** | 0.85 |

Les deux résines sont **très proches** en composition brute, mais tombent sur
des motifs différents, et l'écart s'explique :

- Le TACOT porte **2.35 % d'azote** contre **1.44 %** pour le ZURAM, plus de
  l'eau piégée au moulage. Ces composés piégés diluent le carbone et
  enrichissent H et O — ce qui fait glisser le TACOT du réseau C₇H₆O vers le
  monomère C₆H₆O. Sykes l'écrit explicitement (p. 13).
- Le ZURAM, moins pollué, tombe donc plus près du **réseau réticulé réel**.

Autrement dit : les deux matériaux emploient très vraisemblablement la **même
chimie** (novolac + hexamine) ; c'est la qualité de la cuisson et la quantité
de composés piégés qui distingue les deux analyses, pas la nature du polymère.

---

## 3. Le matériau complet

`'ZURAM_official'!D12:F21` — fractions volumiques et densités intrinsèques :

| | gaz (pores) | fibres | résine | somme |
|---|---|---|---|---|
| vierge | 0.707299 | 0.129973 | 0.162729 | 1.0000 |
| char | 0.769174 | 0.129973 | 0.100854 | 1.0000 |

| | fibres | matrice vierge | matrice charbonnée |
|---|---|---|---|
| densité intrinsèque | 1.5772 | 1.3151 | **1.3151** |

Densités moyennes (`D17:E17`), reproduites exactement :

```
rho_v = 0.129973*1.5772 + 0.162729*1.3151 = 0.418998084502609
rho_c = 0.129973*1.5772 + 0.100854*1.3151 = 0.337628084504206
```

> **Piège d'unités.** L'en-tête annonce `[kg/m³]` mais les valeurs sont en
> **g/cm³** : le ZURAM fait **419.0 kg/m³** vierge et **337.6 kg/m³** charbonné.

Couplage matériau : `k = B'g/B'c = (ρ_v − ρ_c)/ρ_c = **0.2410**`
(TACOT : 0.2727). Le ZURAM libère 81.4 kg/m³ de gaz, soit **19.4 %** de sa
masse vierge (TACOT : 21.4 %).

La résine représente **51.07 % de la masse vierge** (TACOT : 42.9 %) pour
16.3 % du volume.

### Enthalpie de formation — contrôle de cohérence

`'ZURAM_official'!D69:D71` donne char 0, résine phénolique −2 143 700 J/kg,
vierge −1 094 878.036 J/kg. Or :

```
0.510742 × (−2 143 700) + 0.489258 × 0 = −1 094 878.036 J/kg      exact
```

Les fibres de carbone ont `h_f = 0` par convention, et la fraction massique de
résine se déduit de la microstructure. La feuille est donc interne­ment
cohérente à la dernière décimale.

---

## 4. Le gaz de pyrolyse — une chaîne différente de celle du TACOT

C'est le point le plus important pour qui vient du TACOT.

| | TACOT | ZURAM |
|---|---|---|
| donnée primaire | **spéciation moléculaire** (7 espèces, table I de Sykes) | **composition élémentaire massique**, mesurée |
| passage à l'élémentaire | somme des atomes `n_E = Σ ν_{E,j} x_j` | simple conversion masse → mole |
| source | Sykes 1967, chromatographie | VKI, campagne AblaNTIS ([4] du classeur) |

`'ZURAM_official'!C86:F90` :

| Élément | Fraction molaire (D) | Fraction massique (E) | Masse molaire (F) |
|---|---|---|---|
| C | 0.170941536278466 | 0.457 | 12.0107 |
| H | 0.722071254424386 | 0.162 | 1.00794 |
| O | 0.106987209297148 | 0.381 | 15.999 |

La colonne D se déduit **exactement** de la colonne E avec les masses molaires
de la colonne F :

```
x_E = (y_E / M_E) / Σ_j (y_j / M_j)
```

| | y (masse) | x calculé | x feuille | écart |
|---|---|---|---|---|
| C | 0.457 | 0.1709415362784664 | 0.1709415362784660 | 4e-16 |
| H | 0.162 | 0.7220712544243860 | 0.7220712544243860 | 0 |
| O | 0.381 | 0.1069872092971476 | 0.1069872092971480 | 4e-16 |

puis arrondi à trois décimales → **C:0.171, H:0.722, O:0.107**, la valeur des
XML `zuram-air.xml` et `zuram-pyrogas.xml`, et du fichier de référence AblaNTIS
`carbonPhenolInAir_AblaNTIS.xml` (`VKIZuramPyroGas`).

**Masse molaire du gaz : 4.4926 g/mol** (TACOT : 4.9985). Le gaz du ZURAM est
plus léger — plus riche en hydrogène (0.722 contre 0.679) et nettement plus
pauvre en carbone (0.171 contre 0.206).

> **Limite de la traçabilité.** Les fractions massiques 0.457 / 0.162 / 0.381
> ne sont données qu'à **trois décimales**, et leur source ([4]
> `VKI_ZURAM_characterization_v1-5.ods`, reprise dans AblaNTIS TN-2.2) n'est pas
> dans les documents fournis. La chaîne s'arrête donc là — d'où l'intérêt du
> contrôle indépendant du § 6.
>
> Une valeur *préliminaire* figure d'ailleurs dans l'onglet
> `ZURAM_preliminary` : C:0.1246, H:0.7506, O:0.1249 (masse 0.352/0.178/0.470),
> issue d'une autre voie. **Ce n'est pas** la valeur du XML. L'écart entre les
> deux jeux est considérable (37 % sur le carbone) : bien vérifier de quel
> onglet on part.

---

## 5. Le rendement en char

Comme pour le TACOT, **aucune cellule ne porte le rendement en char**. Il est
implicite, et se retrouve par deux voies indépendantes qui coïncident :

**(a) Fractions volumiques**, à densité intrinsèque de matrice inchangée
(`D21` = `F21` = 1.3151) — donc pure perte de volume :

```
0.100853796331154 / 0.1627287372 = 0.619766
```

**(b) Cinétique de pyrolyse**, `'ZURAM_official'!B80:B83`. Le classeur note en
`B77` que les `f` de la source ont été **divisées par la fraction massique de
résine** pour porter sur la résine seule et non sur le composite :

| Réaction | f | log₁₀A | E/R (K) | m |
|---|---|---|---|---|
| 1 | 0.035070 | 5.33 | 8 178.5 | 4.30 |
| 2 | 0.027687 | 8.69 | 16 068.4 | 3.70 |
| 3 | 0.095981 | 10.60 | 21 612.9 | 2.57 |
| 4 | 0.221495 | 11.67 | 26 423.8 | 4.63 |
| **Σ** | **0.380234** | | | |

```
1 − 0.380234 = 0.619766       (écart avec (a) : 4e-16)
```

**Rendement en char de la résine du ZURAM : 0.6198**, contre **0.50** pour le
TACOT. Le ZURAM charbonne nettement plus — cohérent avec une résine plus
réticulée (motif C₇H₆O plutôt que C₆H₆O).

> Les deux classeurs font la même hypothèse microstructurale : densité
> intrinsèque de matrice **inchangée**, la pyrolyse ne retirant que du volume.
> Le classeur en fait d'ailleurs une **anomalie ouverte** (`Versions_and_issues`,
> issue n° 3, statut *Open*) : *« Intrinsic density of charred resin not well
> evaluated by He pycnometry. Hypothesis of constant resin density and only
> volume shrinkage used. »*
> Pour le TACOT, cette lecture coexistait avec une seconde (Goldstein :
> densité divisée par deux à volume constant) ; ici il n'y en a qu'une.

---

## 6. Validation croisée : résine − char → gaz

La composition du gaz vient d'une source VKI ([4]) ; celles de la résine et du
char viennent de la campagne de [THo]. Deux campagnes distinctes : la fermeture
est un vrai contrôle, pas une tautologie.

Partir de la résine (fig. 7, analyseur élémentaire), retrancher le char au
rendement 0.6198 :

| char retranché | C | H | O | écart max / XML |
|---|---|---|---|---|
| **carbone pur** | **0.1628** | **0.7250** | **0.1122** | **4.9 %** |
| mesuré à 800 °C | 0.1857 | 0.7184 | 0.0959 | 10.4 % |
| *cible : le XML* | *0.1709* | *0.7221* | *0.1070* | — |

**Accord à 4.9 %** avec l'hypothèse de char de carbone pur — bon, compte tenu
que la composition de la résine est relevée sur un graphique (±0.5 point) et
provient d'une campagne différente de celle du gaz.

### Le char n'est pas du carbone pur

[THo] § 4.2.1 mesure le char à 800 °C et le souligne explicitement :

> *« the char, even at this high pyrolysis temperature contained relevant
> percentages of the different elements (**C: 96.36 %, N: 0.13 %, H: 0.08 %,
> O: 3.43 %**), in contrast to the typical assumption of 100 % carbon in the
> char »*

L'hypothèse `Char = C:1.0` des XML reste néanmoins licite en ablation : la
surface dépasse largement 2000 K, bien au-delà des 800 °C de l'essai, et les
fibres Calcarb sont déjà du carbone pur. C'est d'ailleurs le cas de char pur
qui ferme le mieux (4.9 % contre 10.4 %) — ce qui suggère que la donnée [4] a
elle aussi été établie sous cette hypothèse.

Sykes faisait la même observation sur le TACOT (char à 850 °C : C 92.65 /
H 0.90 / O 6.45 %) — les deux matériaux se comportent pareillement.

---

## 7. Récapitulatif des sources

| Donnée | Valeur | Localisation exacte |
|---|---|---|
| Matériau | ZURAM® 18/50, carbone/phénolique DLR | `'ZURAM_official'!A4` |
| « 18 » du nom | préforme Calcarb CBCF **18**/2000, 180 kg/m³ nominal | `'Calcarb_official'!D15`, note v4.3.0 |
| « 50 » du nom | ≈ % massique de résine (**déduit**, non énoncé) | déduction — cf. § 1 |
| Renfort | préforme Mersen Calcarb CBCF 18/2000 | `'ZURAM_official'!A4` |
| Catalyseur de la résine | hexamine (HMTA) | [THo] § 4.2 |
| Type de résine | novolac (**déduit**, non déclaré) | déduction — cf. § 1 |
| Composition résine vierge | C 75.21 / H 5.75 / N 1.44 / O 14.13 (% masse) | [THo] fig. 7 |
| Formule brute équivalente | C₇.₀₉H₆.₄₆O ≈ novolac C₇H₆O | dérivée de fig. 7 |
| Fractions volumiques | vierge 0.7073 / 0.1300 / 0.1627 | `'ZURAM_official'!D12:F12` |
| | char 0.7692 / 0.1300 / 0.1009 | `'ZURAM_official'!D13:F13` |
| Densités intrinsèques | fibres 1577.2, matrice 1315.1 kg/m³ | `'ZURAM_official'!D21:F21` |
| ρ vierge / ρ char | 419.0 / 337.6 kg/m³ | `'ZURAM_official'!D17:E17` |
| h_f résine / char / vierge | −2 143 700 / 0 / −1 094 878 J/kg | `'ZURAM_official'!D69:D71` |
| Cinétique de pyrolyse | 2 pseudo-phases, 4 réactions | `'ZURAM_official'!A74:F83` |
| Rendement en char | 0.6198 | implicite : `D12:F13` **et** `B80:B83` |
| Gaz, fractions massiques | C 0.457 / H 0.162 / O 0.381 | `'ZURAM_official'!E87:E89` |
| **Gaz, fractions molaires** | **C:0.171 H:0.722 O:0.107** | `'ZURAM_official'!D87:D89` = les XML |
| Char, composition à 800 °C | C 96.36 / N 0.13 / H 0.08 / O 3.43 | [THo] § 4.2.1 |
| Émissivité vierge / char | 0.8 / 0.9 | `'ZURAM_official'!D95:D96` |

---

## 8. Fichiers liés

| Fichier | Rôle |
|---|---|
| `resine_zuram_verification.py` | rejoue les 5 vérifications de ce document |
| `README.md` | table B', validation AblaNTIS, propriétés du gaz |
| `../tacot_bprime/resine_tacot.md` | même travail pour le TACOT |
| `../tacot_bprime/mise_en_donnees_xml.md` | ce qu'il faut renseigner dans un XML |
| `../data/mixtures/zuram-air.xml` | 20 espèces + C(gr), cas AblaNTIS |
| `../data/mixtures/zuram-pyrogas.xml` | gaz de pyrolyse seul, pour h_g(T, P) |
