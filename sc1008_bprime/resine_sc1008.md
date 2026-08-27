# La résine SC-1008 : choix, type, composition

Document de traçabilité, pendant de `tacot_bprime/resine_tacot.md` et de
`zuram_bprime/resine_zuram.md`. Il rattache chaque nombre à l'une des deux
publications primaires :

> [TS] K. A. Trick, T. E. Saliba, *Mechanisms of the pyrolysis of phenolic resin
> in a carbon/phenolic composite*, **Carbon 33 (1995) 1509–1515**.
>
> [Wo] H.-W. Wong, J. Peck, J. Assif, F. Panerai, J. Lachaud, N. N. Mansour,
> *Detailed analysis of species production from the pyrolysis of the Phenolic
> Impregnated Carbon Ablator*, **J. Anal. Appl. Pyrolysis 122 (2016) 258–267**.

Toutes les vérifications numériques de ce document sont rejouables :

```bash
python sc1008_bprime/resine_sc1008_verification.py
```

---

## 0. Réponse courte

| Question | Réponse |
|---|---|
| **Type** | **résol** (resole) phénol-formaldéhyde, monocomposant, réticulant par chauffage seul |
| **Référence** | **Durite™ SC-1008**, Hexion / Bakelite Synthetics (ex-Borden) |
| **Livraison** | solution : **20–25 % d'isopropanol**, **11–18 % de phénol libre**, **0.6–2 % de formaldéhyde libre** ([Wo] § 2.2) |
| **Emploi** | matrice du **PICA**, sur préforme carbone FiberForm |
| **Motif retenu** | **C₇.₅H₆O** — 80.3 % C / 5.4 % H / 14.3 % O en masse |
| **Azote** | **aucun** — un résol n'a pas de durcisseur, donc pas d'hexamine |
| **Rendement en char** | **0.55** (enveloppe 0.545 – 0.62) |
| **Ce qui entre dans le XML** | `C:0.2526, H:0.6407, O:0.1068` |

Ce qui finit dans `data/mixtures/sc1008-*.xml` n'est pas la résine mais le
**gaz** qu'elle produit. Le § 4 montre que ces trois nombres reproduisent la
composition mesurée de [TS] à **1.5 % près sur C, 0.3 % sur H, 1.9 % sur O**.

---

## 1. Identité de la résine

| | Valeur | Source |
|---|---|---|
| Nature | résol phénol-formaldéhyde | [Wo] § 2.2, « resole PF resin » partout |
| Référence commerciale | **Durite SC-1008** | [Wo] § 2.2 |
| Fournisseur | Hexion, Batesville (AR), USA | [Wo] § 2.2 |
| Isopropanol (solvant) | **20–25 %** masse | [Wo] § 2.2 |
| Phénol libre | **11–18 %** masse | [Wo] § 2.2 |
| Formaldéhyde libre | **0.6–2 %** masse | [Wo] § 2.2 |
| Mise en œuvre PICA | imprégnation, cuisson à chaud, séchage en étuve sous vide (conditions propriétaires) | [Wo] § 2.2 |
| Densité PICA | 0.274 g/cm³ ± 10 %, porosité 0.8 ± 10 % | [Wo] § 2.2 |

[Wo] § 2.2 :

> *« A typical SC-1008 resin contains 0.6–2 weight percent of formaldehyde,
> 11–18 weight percent of phenol, and 20–25 weight percent of isopropyl alcohol
> as a solvent. »*

### Résol, et ce que ça implique

Un résol est formé **en excès de formaldéhyde** (F/P > 1) sous catalyse basique.
Il porte donc déjà les méthylols nécessaires à sa propre réticulation : il durcit
à la chaleur **sans durcisseur ajouté**.

C'est la différence de fond avec le TACOT et le ZURAM, qui sont des **novolac**
durcis à l'hexaméthylènetétramine. Conséquence directe pour la mise en données :

- le TACOT mesure **2.35 % d'azote** (Sykes table II), le ZURAM **1.4 %** ;
- **le SC-1008 n'en a aucun** — il n'y a pas de HMTA dans la formulation.

L'argument que le TACOT doit construire (l'azote est piégé, il s'effondre dès
400 °C, on peut donc le retirer du gaz) est ici **sans objet** : le gaz est
nativement C/H/O.

### Les 11–18 % de phénol libre

C'est une particularité du SC-1008 qu'aucune des deux autres résines du dépôt ne
partage, et elle a une conséquence mesurable : ce phénol part **en région 1**
(300–550 °C), sans provenir du réseau. Il explique en partie que [TS] trouve
**14 % de phénol/crésol** dans le gaz global (§ 4).

---

## 2. Le motif — base de calcul

Une seule hypothèse : **une résine résol cuite = un cycle phénolique + ses ponts
méthylène**.

```
phénol       C6H6O
+ 1.5 CH2O            (formaldéhyde)
− 1.5 H2O             (eau de condensation, libérée à la cuisson)
──────────────────────
   C7.5 H6 O
```

Le 1.5 n'est pas un réglage libre : le phénol a **3 sites réactifs** (2 ortho +
1 para), un pont est partagé entre 2 cycles, donc **1.5 pont par cycle** à réseau
saturé.

| | moles | × M | masse |
|---|---|---|---|
| C | 7.5 | 12.011 | 90.08 |
| H | 6 | 1.008 | 6.05 |
| O | 1 | 15.999 | 16.00 |
| | | | **112.13** |

> **C 80.3 % · H 5.4 % · O 14.3 %**

### Généralisation, et pourquoi le choix résol/novolac importe peu

Le motif cuit général est `C(6+b)H6O`, où `b` est le nombre de ponts par cycle.
Résol et novolac ne sont pas deux points distincts : ils **glissent sur la même
courbe**, seule la valeur de `b` accessible diffère (≈ 1 pour un novolac
linéaire, ≈ 1.5 pour un résol saturé).

| Motif | C % | H % | O % |
|---|---|---|---|
| novolac linéaire C₇H₆O (*ZURAM*) | 79.2 | 5.7 | 15.1 |
| **résol saturé C₇.₅H₆O (*retenu*)** | **80.3** | **5.4** | **14.3** |
| *dispersion entre deux novolacs **mesurés*** | *77.4 → 79.1* | | |

L'écart résol/novolac (**1.1 point sur C**) est **plus petit que la dispersion
entre deux novolacs réellement mesurés** (TACOT 77.4 vs ZURAM 79.1, soit 1.7
point). L'analyse élémentaire C/H/O ne permet donc pas d'identifier le type de
résine — seul l'azote le fait.

---

## 3. L'invariant H/O

C'est l'outil central de ce document, et il ne coûte rien.

Avec un **char de carbone pur**, seul le carbone est prélevé à la résine : H et O
traversent la pyrolyse sans être touchés. Donc

```
    H/O (gaz)  =  H/O (résine)          indépendamment du rendement en char
```

| Y | C | H | O | H/O |
|---|---|---|---|---|
| 0.30 | 0.402 | 0.513 | 0.086 | **6.00** |
| 0.50 | 0.288 | 0.610 | 0.102 | **6.00** |
| 0.575 | 0.234 | 0.657 | 0.110 | **6.00** |
| 0.70 | 0.121 | 0.753 | 0.126 | **6.00** |

Le rendement en char — le paramètre le plus incertain de toute la chaîne — n'a
**aucun effet** sur H/O. C'est donc un test de cohérence qui ne dépend d'aucune
des hypothèses fragiles.

### Les bornes

Aucun motif phénolique concevable ne descend sous **3.6** :

| Motif | H/O |
|---|---|
| résol prépolymère non cuit C₇.₅H₉O₂.₅ | 3.60 |
| résol mi-cuit C₇.₅H₆.₆₈O₁.₃₄ | 4.99 |
| résol / novolac cuit | 6.00 |
| TACOT mesuré C6.32H6.10O (Sykes) | 6.10 |

Et **ajouter de l'eau ne sauve rien** : `C7.5H6O·(H2O)_w` donne
`H/O = (6+2w)/(1+w)`, qui tend vers **2.0 par le haut**. La valeur 2.0 est une
borne asymptotique inatteignable.

### Application au dépôt

| Jeu | H/O | Verdict |
|---|---|---|
| `tacot_pyro` | 5.90 | cohérent phénolique |
| `cph70_pyro` | 5.90 | cohérent phénolique |
| `VKIZuramPyroGas` | 6.75 | cohérent phénolique |
| **SC-1008 (ce document)** | **6.00** | par construction |
| **[TS], mesuré** | **6.10** | **référence expérimentale** |
| `cork_pyro` | 4.89 | hors test — le liège (subérine, lignine, polysaccharides) est nativement bien plus oxygéné qu'un phénolique |
| ~~`pica_pyro` (Goldstein)~~ | **1.93** | **incohérent — supprimé du dépôt** |

Le jeu PICA hérité de Goldstein tombait sous la borne asymptotique 2.0 : il
n'était reproductible par **aucune** résine phénolique, à **aucun** degré de
cuisson, avec **aucune** quantité d'eau ajoutée. C'est ce qui a motivé sa
suppression.

---

## 4. Le gaz de pyrolyse : reconstruction et validation

### 4.1 Ce que [TS] donne réellement

**Attention à la nature de la donnée.** Le travail expérimental propre de [TS]
est la **FTIR** sur un préimprégné carbone/phénolique Amoco T300 — la résine n'y
est pas nommée. Les données d'évolution gazeuse **ne sont pas leur mesure** :

> *« Gas evolution data from the pyrolysis of a neat phenolic resin is available
> in the literature [2,3]. This data was used to determine the evolution rates of
> seven identified pyrolysis products »*

[TS] les **déconvolue** en trois régions réactionnelles. Ce n'est donc ni du
SC-1008, ni une mesure primaire : c'est une **ré-analyse de données publiées sur
résine phénolique nue**. Pour une traçabilité de même niveau que Sykes pour le
TACOT, il faudrait remonter aux réfs [2,3] de l'article.

**Table 2** — % molaire *dans* chaque région :

| Région 1 (300–550 °C) | Région 2 (400–800 °C) | Région 3 (560–900 °C) |
|---|---|---|
| H₂O 49.8 | H₂ 59.4 | H₂ 85.7 |
| LMS 50.1 | CH₄ 14.9 | CO 9.5 |
| CO₂ 0.1 | CO 12.7 | H₂O 4.7 |
| | H₂O 12.7 | CO₂ 0.1 |
| | CO₂ 0.2 · C₂H₆ 0.1 | |

**Table 4** — distribution de chaque espèce *sur* les régions (% du total) :

| Espèce | R1 | R2 | R3 |
|---|---|---|---|
| H₂ | — | 57.1 | 42.9 |
| H₂O | 66.3 | 28.2 | 5.5 |
| Phénol/crésol (LMS) | 100 | — | — |
| CO | — | 71.9 | 28.1 |
| CH₄ | — | 100 | — |
| CO₂ | 20.5 | 66.8 | 12.7 |
| C₂H₆ | — | 100 | — |

### 4.2 Reconstruction de la composition globale

Aucune des deux tables ne donne la composition globale. Mais **elles se
recoupent** : pour une espèce *j* présente dans les régions *i* et *k*,

```
    f_ij · N_i / (f_kj · N_k) = d_ij / d_kj      →      N_i / N_k
```

Chaque espèce commune fournit une estimation **indépendante** des tailles
relatives des régions :

| | N₂/N₃ |
|---|---|
| via H₂O | 1.8975 |
| via H₂ | 1.9203 |
| via CO | 1.9140 |
| *via CO₂* | *2.6299 — écarté* |

Les trois premières s'accordent à **1.2 %**. C'est le test : deux tableaux
construits indépendamment redonnent le même découpage. CO₂ est écarté parce qu'il
pèse 0.1–0.2 % dans la table 2 — l'arrondi à la première décimale y domine le
signal.

H₂O est la seule espèce commune aux régions 1 et 2 : N₁/N₂ = 0.5996.

```
    N1 : N2 : N3  =  1.146 : 1.911 : 1.000
```

**Composition moléculaire globale reconstruite :**

| H₂ | H₂O | LMS | CO | CH₄ | CO₂ | C₂H₆ |
|---|---|---|---|---|---|---|
| 0.4911 | 0.2121 | 0.1415 | 0.0832 | 0.0702 | 0.0015 | 0.0005 |

### 4.3 Du moléculaire à l'élémentaire

LMS (« low molecular weight substances ») = phénol + crésol, non résolus par
[TS]. **C'est la seule hypothèse vraiment libre du calcul.** Sensibilité :

| Hypothèse LMS | C | H | O | H/O |
|---|---|---|---|---|
| 100 % phénol C₆H₆O | 0.2523 | 0.6374 | 0.1104 | 5.77 |
| **50/50 C₆.₅H₇O** | **0.2564** | **0.6388** | **0.1048** | **6.10** |
| 100 % crésol C₇H₈O | 0.2601 | 0.6402 | 0.0998 | 6.42 |

**L'invariant H/O ≈ 6 tient sur toute la plage.**

### 4.4 Confrontation à la fermeture élémentaire

Le calcul du § 4.2–4.3 ne fait intervenir **ni le motif de la résine, ni le
rendement en char** : il ne sort que des deux tableaux publiés. La fermeture est
donc une **prédiction indépendante**.

| | C | H | O | H/O |
|---|---|---|---|---|
| **[TS] reconstruit (mesure)** | **0.2564** | **0.6388** | **0.1048** | **6.10** |
| **fermeture C₇.₅H₆O, Y = 0.55** | **0.2526** | **0.6407** | **0.1068** | **6.00** |
| **écart** | **1.5 %** | **0.3 %** | **1.9 %** | **1.6 %** |

Deux chaînes de raisonnement entièrement disjointes — l'une partant d'une mesure
de gaz déconvoluée en régions, l'autre d'un motif chimique idéalisé et d'un bilan
de masse — convergent à **2 %**.

Sensibilité au rendement en char, à motif fixé :

| Y | écart sur C vs [TS] |
|---|---|
| 0.500 | 12.4 % |
| **0.550** | **1.5 %** |
| 0.575 | 8.9 % |
| 0.600 | 16.8 % |

---

## 5. Le rendement en char

### 5.1 Trois sources indépendantes

| Source | Y | Nature |
|---|---|---|
| [TS], via la composition du gaz reconstruite | **0.545** | composition élémentaire |
| ATG SC-1008 nue (55–60 % au-delà de 650 °C) | 0.55 – 0.60 | thermogravimétrie |
| [Wo], PICA : 19 % de perte / 50 % de résine | **0.620** | bilan de masse du composite |

[Wo] § 3.2 :

> *« the cumulative mass loss of PICA reaches ≈19 % after 1250 K »*

et, sur la fraction de résine :

> *« phenolic resin consists of approximately half of the PICA mass, given that
> the carbon preform does not undergo mass loss during pyrolysis »*

d'où Y(PICA) = (50 − 19)/50 = **0.62**.

### 5.2 Valeur retenue : 0.55

**Ce n'est pas la moyenne des trois.** Le milieu arithmétique (0.575) s'écarte de
[TS] de 8.9 % (§ 4.4). La valeur 0.55 est la seule qui satisfasse les **deux
contraintes fortes simultanément** :

- elle tombe dans la fourchette ATG (0.55–0.60) ;
- elle reproduit la composition de [TS] à 1.5 %.

Les trois sources ne mesurent pas la même chose : moyenner n'a pas de sens.

La borne haute **0.62 est propre au PICA** et ne vaut pas pour la résine seule :
la préforme carbone capte des radicaux pendant la pyrolyse, ce qui relève le
rendement. [Wo] l'observe directement — H₂ divisé par ~5 et hydrocarbures légers
abaissés par rapport à la résole nue :

> *« suggesting that carbon preform in PICA may scavenge the hydrocarbon radicals
> produced during the pyrolysis and charring process »*

---

## 6. Pourquoi [Wo] n'est pas utilisé comme composition élémentaire

C'est pourtant la mesure la plus directe : SC-1008 réel, spéciation GC détaillée,
16 espèces, 20 paliers de 323 à 1252 K. Mais sa composition élémentaire intégrée
est inexploitable telle quelle.

| | C | H | O | H/O | masse |
|---|---|---|---|---|---|
| **brut, table 1 intégrée** | 0.120 | 0.642 | **0.238** | **2.70** | 24.5 mg/100 mg |
| A : − eau T ≤ 664 K | 0.138 | 0.639 | 0.223 | 2.86 | 21.2 mg |
| B : + eau > 850 K (attribution de [Wo]) | 0.222 | 0.622 | 0.156 | **3.98** | 13.0 mg |
| C : calage sur la balance | 0.154 | 0.635 | 0.211 | **3.02** | 19.0 mg |
| *fermeture, Y = 0.62* | *0.197* | *0.689* | *0.115* | *6.00* | *19.0 mg* |

Trois causes, **toutes documentées par les auteurs eux-mêmes** :

**(a) Humidité.** PICA est fortement hygroscopique :

> *« The highly porous architecture of the nano-dispersed phenolic is responsible
> for adsorption of atmospheric moisture »*

> *« despite our best effort to keep the samples in vacuum environment, there is
> unavoidable atmospheric moisture leaking into the system at room temperature
> between runs »*

L'eau au-delà de 850 K est explicitement attribuée à de la désorption. Le GC
totalise **24.5 mg** contre **19 mg** à la balance : l'excédent de 5.5 mg vaut
exactement **0.306 mmol d'eau**.

**(b) Carbone manquant.** Le GC est limité aux produits de **M < 400 g/mol**.

```
    C attendu par fermeture : 0.763 mmol
    C mesuré au GC          : 0.499 mmol   (65 %)
    →  35 % du carbone échappe à la mesure
```

Cohérent avec la structure du gaz : les aromatiques ne pèsent que **3.0 % des
moles** mais **55 % du carbone mesuré**. C'est la queue lourde qui manque, et
elle est presque tout le carbone.

**(c) Composite ≠ résine.** La préforme déplace la spéciation (§ 5.2).

Les corrections d'humidité vont dans le bon sens (H/O 2.70 → 3.98) mais ne
suffisent pas ; le carbone manquant explique le reste.

> **Conclusion.** [Wo] est la meilleure source pour la **spéciation** (quelles
> espèces, à quelle température) et pour le **rendement en char**. Ce n'est pas
> une source de composition élémentaire : sa fermeture n'est pas assurée.

---

## 7. Ce qui entre dans le XML

```xml
<element_compositions default="air">
    <composition name="air">N:0.79, O:0.21</composition>
    <composition name="sc1008_pyro">C:0.2526, H:0.6407, O:0.1068</composition>
    <composition name="sc1008_char">C:1.0</composition>
</element_compositions>
```

```bash
bprime -T 300:25:4000 -P 101325 -b 0.5 \
       -m sc1008-air -bl air -py sc1008_pyro -char sc1008_char -char-elem C
```

Équivalent en fractions massiques (Mutation++ normalise seul, mais `type="mass"`
est alors obligatoire — sans lui le parseur suppose du molaire) :

```xml
<composition name="sc1008_pyro" type="mass">C:0.563, H:0.120, O:0.317</composition>
```

### Vérification bout en bout

Les deux mélanges ont été chargés et exercés avec les trois outils :

```
checkmix sc1008-air        27 espèces, 4 éléments        OK
checkmix sc1008-pyrogas    18 espèces, 3 éléments        OK
bprime   sc1008-air        table B'c(T) complète          OK
mppequil sc1008-pyrogas    h_g(T, P), M, Cp, gamma        OK
```

Ordres de grandeur du gaz de pyrolyse, à 1 atm, comparés au TACOT :

| T (K) | M SC-1008 (g/mol) | M TACOT | h SC-1008 (MJ/kg) | h TACOT |
|---|---|---|---|---|
| 500 | 25.69 | 21.95 | −5.45 | −6.71 |
| 1000 | 18.69 | 13.94 | −2.43 | −2.17 |
| 2000 | 12.64 | 11.02 | +6.01 | +5.00 |
| 3000 | 11.84 | 10.32 | +11.98 | +11.51 |

Le gaz du SC-1008 est plus lourd et moins exothermique à basse température : il
est plus carboné (C:0.253 contre 0.206) et moins riche en eau que celui du
TACOT — cohérent avec un résol saturé face à une résine diluée par l'eau piégée
au moulage (cf. `../tacot_bprime/resine_tacot.md` § 2).

---

## 8. Récapitulatif des sources

| Donnée | Valeur | Localisation exacte |
|---|---|---|
| Type de résine | résol phénol-formaldéhyde | [Wo] § 2.2 |
| Référence commerciale | Durite SC-1008, Hexion | [Wo] § 2.2 |
| Isopropanol / phénol libre / formaldéhyde libre | 20–25 / 11–18 / 0.6–2 % masse | [Wo] § 2.2 |
| Densité PICA / porosité | 0.274 g/cm³ ± 10 % / 0.8 ± 10 % | [Wo] § 2.2 |
| Fraction massique de résine dans PICA | ≈ 50 % | [Wo] § 3.2 |
| Perte de masse PICA à 1250 K | 19 % (balance), 17 % (ATG) | [Wo] § 3.2 |
| Spéciation du gaz, 16 espèces, 20 paliers | mmol/100 mg PICA | [Wo] table 1 |
| Régions réactionnelles et températures | 300–550 / 400–800 / 560–900 °C | [TS] table 3 |
| Composition par région | % molaire | [TS] table 2 |
| Distribution des espèces sur les régions | % du total | [TS] table 4 |
| Motif de résine retenu | C₇.₅H₆O | dérivé, § 2 |
| Rendement en char retenu | 0.55 | recoupement [TS] / ATG, § 5.2 |
| **Gaz, composition élémentaire** | **C:0.2526 H:0.6407 O:0.1068** | fermeture, § 4.4 = les XML |

---

## 9. Fichiers liés

| Fichier | Rôle |
|---|---|
| `resine_sc1008_verification.py` | rejoue les 5 vérifications de ce document |
| `../tacot_bprime/resine_tacot.md` | même travail pour le TACOT (novolac, Sykes 1967) |
| `../zuram_bprime/resine_zuram.md` | même travail pour le ZURAM (novolac + hexamine) |
| `../tacot_bprime/mise_en_donnees_xml.md` | ce qu'il faut renseigner dans un XML de mélange |
| `../tacot_bprime/pyrolysis_gas_from_resin.py` | fermeture résine → gaz, sensibilité au rendement |
| `../../data/mixtures/sc1008-air.xml` | ablation dans l'air, avec C(gr) |
| `../../data/mixtures/sc1008-pyrogas.xml` | gaz de pyrolyse seul, pour h_g(T, P) |
