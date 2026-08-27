# Mise en données XML d'un liège/phénolique (cork phenolic)

**Matériau traité** : 80 % liège / 20 % résine phénolique (fractions
**massiques** du solide), résine de rendement en char **50 %**, et rendement
en char **du composite 20 %** — valeur mesurée par TGA sur le cork P50, dont
se déduit le rendement du liège, **12.5 %** (§3.1, §7).

Objectif : produire les trois données que consomment `bprime` et `mppequil`
pour la table **B'c(T, P, B'g)**, l'enthalpie de paroi **h_w** et l'enthalpie
du gaz de pyrolyse **h_g**.

Pour la mécanique générale du fichier XML (balises, parseur, liste d'espèces,
pièges), voir `../tacot_bprime/mise_en_donnees_xml.md`. Le présent document
traite ce qui est **spécifique au liège** : *le renfort pyrolyse aussi*.

---

## 1. Ce qui change par rapport à un carbone/phénolique

Pour un TACOT, un PICA, un Zuram ou un CPh70, le renfort est de la fibre de
carbone : **inerte à la pyrolyse**. D'où deux conséquences bien connues
(`../tacot_bprime/cph70_vs_tacot.md`) :

- le gaz de pyrolyse est produit par la **résine seule** ;
- fibres et résine carbonisée donnent le **même** char (C pur) ;
- donc le **rapport renfort/résine n'entre pas** dans la table B'.

Pour un liège/phénolique, le liège pyrolyse — c'est même l'essentiel du gaz.
Il faut donc :

| | carbone/phénolique | liège/phénolique |
|--|--|--|
| gaz de pyrolyse | résine seule | **liège + résine**, au prorata des masses dégazées |
| rapport renfort/résine dans la table B' | non | **oui** (il pondère le mélange des deux gaz) |
| rendement en char du renfort | sans objet (0 % de gaz) | **donnée d'entrée indispensable** (ici déduite de la TGA du composite) |
| couplage stationnaire k = B'g/B'c | 0.27 (TACOT) | **4.0** — quinze fois plus de soufflage |

La dernière ligne est la plus lourde de conséquences : le matériau ne lit pas
la table au même endroit du tout (§6).

---

## 2. Les trois compositions élémentaires à déterminer

`bprime` ne consomme que des **fractions massiques élémentaires**
(`src/apps/bprime.cpp:304-330`) :

| # | Composition | Option | Ici |
|---|-------------|--------|-----|
| 1 | bord de couche limite | `-bl` | air `N:0.79, O:0.21` |
| 2 | gaz de pyrolyse | `-py` | à calculer — §3 |
| 3 | char | `-char` + `-char-elem` | `C:1.0` (§5 pour le cas multi-élément) |

Aucune grandeur volumique (ρ, porosité, épaisseur) n'entre dans ce bilan.

---

## 3. Le gaz de pyrolyse par fermeture élémentaire

La pyrolyse ne crée ni ne détruit d'atomes. Constituant par constituant :

```
constituant = char + gaz      =>      n_E(gaz, i) = n_E(i) - n_E(char, i)
```

puis on **somme les deux gaz** — c'est là qu'entre le rapport 80/20.

### 3.1 Données d'entrée

#### Le liège : analyse élémentaire retenue

```
liège :  C 62.4   H 8.5   O 28.4   (% masse, sec sans cendres)
```

Valeurs usuellement citées pour le liège de *Quercus suber*. L'azote
(~0.6 % masse) est négligé — le retenir ajouterait N au gaz de pyrolyse, donc
une quatrième composante à `cork_pyro`.

> **Provenance à confirmer.** Ces trois nombres sont des valeurs usuelles, pas
> une donnée vérifiée sur une source primaire : les domaines d'éditeurs sont
> inaccessibles depuis l'environnement de calcul, et aucune des recherches
> menées n'a permis de remonter à un tableau d'analyse élémentaire publié pour
> *Q. suber*. À confronter à une analyse (publiée ou mesurée) avant tout usage
> engageant. Pistes : Pereira H., *Cork: Biology, Production and Uses*
> (Elsevier, 2007) ; les articles de Şen & Pereira sur la pyrolyse du liège,
> qui portent une analyse élémentaire du liège cru en plus de celle des
> biochars. Le contrôle biochimique ci-dessous, lui, est traçable ligne à
> ligne.

#### Contrôle : le liège reconstruit depuis ses unités de répétition

On additionne les constituants, chacun représenté par son unité de
répétition :

| Constituant | % masse | Unité retenue | Formule | M | C | H | O |
|---|---|---|---|---|---|---|---|
| subérine | 45 | glycérol + 3 acides 9,10-époxy-18-hydroxyoctadécanoïques, − 3 H₂O | C57H104O12 | 981.5 | 69.8 | 10.7 | 19.6 |
| lignine | 27 | alcool coniférylique (unité guaïacyle) | C10H12O3 | 180.2 | 66.7 | 6.7 | 26.6 |
| cellulose, polysaccharides | 12 | anhydroglucose | C6H10O5 | 162.1 | 44.4 | 6.2 | 49.3 |
| tanins | 6 | flavan-3-ol (catéchine) | C15H14O6 | 290.3 | 62.1 | 4.9 | 33.1 |
| céroïdes | 6 | friedéline (triterpène des cires du liège) | C30H50O | 426.7 | 84.4 | 11.8 | 3.7 |

Les parts somment à **96 %** : le complément (cendres ~1 %, humidité
résiduelle, extractibles mineurs) n'est pas de la matière C/H/O identifiée. On
renormalise sur les 96 % déclarés.

```
liège reconstruit  : C 66.16   H 8.71   O 25.13   (% masse)
```

| | C | H | O |
|---|---|---|---|
| **littérature — retenue** | **62.40** | **8.50** | **28.40** |
| reconstruit | 66.16 | 8.71 | 25.13 |
| **écart** | **+3.8** | **+0.2** | **−3.3** |

L'hydrogène tombe juste à 0.2 point près ; le carbone est surestimé d'environ
4 points, l'oxygène sous-estimé d'autant. Deux causes plausibles, non
exclusives :

- **les 4 % manquants**. S'ils sont de la matière oxygénée (humidité
  résiduelle, sels), la reconstruction rejoint la mesure à ~1 point.
- **le choix de l'unité de subérine**, qui pèse 45 % à elle seule. La variante
  sans glycérol (C18H32O3, acide époxy-hydroxy estérifié) donne C 67.6 / H 8.8
  / O 23.6, soit 1.5 point de C de plus : la subérine réelle contient plus de
  glycérol, de diacides et d'acide férulique que le modèle, tous plus
  oxygénés.

**Décision** : on retient l'**analyse de littérature**. La reconstruction
reste dans le script (`cork_elemental()`) comme contrôle indépendant — il
suffit d'écrire `CORK_MASS_PCT = cork_elemental()` pour basculer. Les deux
voies se recoupent à 4 points sur C, ce qui est un accord raisonnable vu les
choix d'unités ; l'effet sur le gaz est chiffré en §3.3.

#### La résine

| Constituant | Composition | Rendement en char |
|---|---|---|
| résine phénolique | novolac **C7H6O** (M = 106.12 g/mol) | 50 % *(énoncé)* |

Même motif que pour le TACOT, le PICA et le CPh70 — c'est la « composition
classique » des autres matériaux du dépôt.

#### Les rendements en char

| | valeur | origine |
|---|---|---|
| composite | **20 %** | mesuré (TGA du cork P50, §7) |
| résine | 50 % | énoncé |
| **liège** | **12.5 %** | **déduit** : `0.80·y + 0.20·0.50 = 0.20` |

**Contrôle par additivité des constituants.** Chaque famille a son propre
rendement (ordres de grandeur, ATG lente sous inerte) :

| | subérine | lignine | polysacch. | tanins | céroïdes | → liège |
|---|---|---|---|---|---|---|
| % masse | 45 | 27 | 12 | 6 | 6 | |
| char | 15 % | 45 % | 15 % | 40 % | 2 % | **24.2 %** |

Soit un composite à **29.3 %**, contre **20 %** mesurés sur le P50 — les deux
voies ne se rejoignent pas, d'un facteur deux sur le liège. Explications
possibles : le P50 est **plastifié au glycol**, qui part entièrement en gaz ;
son rapport liège/résine n'est pas nécessairement 80/20 ; et les rendements
par constituant ci-dessus sont indicatifs. **On garde la mesure** — c'est elle
qui fixe `k`, la grandeur qui gouverne le point de fonctionnement. L'écart est
chiffré en §3.3 et borné en §7.3.

### 3.2 Bilan sur 100 g de composite vierge

| | masse | char | gaz | moles d'atomes du gaz |
|---|---|---|---|---|
| liège | 80 g | 10 g | 70 g | C 3.353, H 6.794, O 1.430 |
| résine | 20 g | 10 g | 10 g | C 0.487, H 1.131, O 0.188 |
| **composite** | 100 g | **20 g** | **80 g** | **C 3.840, H 7.924, O 1.619** |

Fractions molaires élémentaires :

| | C | H | O |
|---|---|---|---|
| gaz du liège seul | 0.290 | 0.587 | 0.124 |
| gaz de la résine seule | 0.269 | 0.626 | 0.104 |
| **gaz du composite** | **0.287** | **0.592** | **0.121** |
| *pour mémoire, TACOT* | *0.206* | *0.679* | *0.115* |

Soit en masse : C 0.576, H 0.100, O 0.324 (le parseur convertit tout seul,
cf. `mise_en_donnees_xml.md` §3.2 — mais il faut déclarer le bon `type`).

Calcul, sensibilités et variantes : `python cork_pyrolysis_data.py`.

### 3.3 Ce que valent les deux révisions

Effet, sur le gaz, de la reconstruction du liège et de l'incertitude sur le
rendement en char :

| composition du liège | char liège | gaz C / H / O | k |
|---|---|---|---|
| **littérature 62.4/8.5/28.4 — retenue** | **12.5 %** | **0.287 / 0.592 / 0.121** | **4.00** |
| reconstruite 66.2/8.7/25.1 | 12.5 % | 0.300 / 0.594 / 0.107 | 4.00 |
| littérature | 24.2 % (additivité) | 0.243 / 0.629 / 0.128 | 2.41 |

Passer d'une composition de liège à l'autre déplace le gaz de **0.013 sur C et
0.014 sur O** (+4 % et −12 % en relatif) : réel, mais second ordre.
L'incertitude sur le rendement en char pèse **trois fois plus** sur la
composition, et surtout **divise k par 1.7** — c'est-à-dire qu'elle déplace le
point de fonctionnement bien davantage. Priorité de mesure : la TGA d'abord,
l'analyse élémentaire ensuite.

### 3.4 Ce que le rapport 80/20 change vraiment

Les deux gaz se ressemblent — celui du liège est un peu plus carboné et plus
oxygéné, celui de la résine un peu plus hydrogéné — donc **le mélange déplace
peu la composition** :

| liège (% masse) | C | H | O | k = B'g/B'c |
|---|---|---|---|---|
| 0 % (résine pure) | 0.269 | 0.626 | 0.104 | 1.000 |
| 50 % | 0.282 | 0.602 | 0.116 | 2.200 |
| **80 %** | **0.287** | **0.592** | **0.121** | **4.000** |
| 100 % (liège pur) | 0.290 | 0.587 | 0.124 | 7.000 |

Conclusion honnête : le rapport liège/résine entre bien dans la table B'
(contrairement au cas carbone/phénolique), mais son effet **direct** sur la
composition du gaz est modeste — les deux gaz se ressemblent. Son effet
**dominant** est sur la *quantité* de gaz, c.-à-d. sur `k`, qui est multiplié
par sept entre la résine pure et le liège pur.

---

## 4. Les deux fichiers XML

### 4.1 `data/mixtures/cork-air.xml` — table B' et h_w

```xml
<mixture thermo_db="NASA-9">

    <species>
       C H O N CH4 CN CO CO2 C2 C2H C2H2,acetylene C3 C4 C4H2,butadiyne C5
       HCN H2 H2O N2 CH2OH CNN CNC CNCOCN C6H6 HNC C(gr)
    </species>

    <element_compositions default="air">
        <composition name="air">N:0.79, O:0.21</composition>
        <composition name="cork_pyro">C:0.287, H:0.592, O:0.121</composition>
        <composition name="cork_char">C:1.0</composition>
    </element_compositions>

</mixture>
```

- **`C(gr)` est obligatoire** : sans phase condensée l'excès de carbone du char
  ne peut pas se condenser, Y_w,C → 1 et B'c sature sur une valeur non
  physique (≈ 200 au lieu de ≈ 0.087 à 300 K).
- 25 espèces + C(gr) : la liste de la table de référence TACOT 3.0. Elle
  couvre les mêmes éléments C/H/O/N et convient donc telle quelle. Elle
  conditionne le résultat (`mise_en_donnees_xml.md` §4.2) : pour comparer à
  une table publiée, reprendre **sa** liste.
- Rien d'autre : ni densité (466 kg/m³), ni porosité, ni rendement en char.
  Ces grandeurs vivent dans le solveur de réponse matériau.

### 4.2 `data/mixtures/cork-pyrogas.xml` — enthalpie du gaz h_g

Le gaz de pyrolyse **pur**, à l'équilibre : ni air, ni char, **ni C(gr)**.

```xml
<mixture thermo_db="NASA-9">
    <species>
       C H O CH4 CO CO2 C2 C2H C2H2,acetylene C3 C4 C4H2,butadiyne C5
       H2 H2O CH2OH C6H6
    </species>
    <element_compositions default="cork_pyro">
        <composition name="cork_pyro">C:0.287, H:0.592, O:0.121</composition>
    </element_compositions>
</mixture>
```

Deux points d'attention :

- **pas de phase condensée** : c'est la convention des tables de propriétés de
  gaz de pyrolyse (classeur TACOT 3.0, feuille *Pyrolysis model*). Ajouter
  C(gr) ferait précipiter du carbone solide et changerait h_g.

  **Ce n'est pas une hypothèse gratuite, et elle coûte cher ici.** Mesuré en
  ajoutant C(gr) au mélange, à 1 atm et 1000 K :

  | gaz | x(C(gr)) à l'équilibre |
  |---|---|
  | **liège/phénolique, C:0.287** | **0.329** |
  | *variante reconstruite, C:0.300* | *0.362* |
  | *TACOT, C:0.206* | *0.200* |

  Le carbone condense dans tous les cas — la convention est donc universelle,
  pas propre au liège — mais l'écart grandit avec la teneur en C : h_g à
  1000 K passe de **−2137 kJ/kg** (sans C(gr), retenu) à **−1434 kJ/kg**
  avec, soit **0.7 MJ/kg**. Si le modèle de réponse matériau traite la suie,
  c'est cette seconde valeur qu'il faut, pas la première.
- **pas d'azote** : le gaz n'en contient pas ; les espèces azotées sont
  retirées de la liste.

---

## 5. Char : quand `C:1.0` ne suffit plus

Le cas de base suppose que liège et résine carbonisent tous deux vers du
carbone. Un char de liège réel retient un peu d'oxygène et d'hydrogène. Le
char devient alors **multi-élément** et se pondère par les **masses de char de
chaque constituant** :

```
n_E(char) = Σ_i (m_char,i / M_i) · ν_{E,i}
```

Exemple calculé par le script (char de liège 90 % C / 2 % H / 8 % O en masse,
char de résine C pur, mêmes rendements) :

| | C | H | O |
|---|---|---|---|
| char | 0.864 | 0.108 | 0.027 |
| gaz (recalculé, les atomes retenus manquent au gaz) | 0.297 | 0.585 | 0.119 |

```xml
<composition name="cork_char">C:0.864, H:0.108, O:0.027</composition>
```

`-char-elem C` reste valable : l'élément du bilan doit être présent dans le
char et absent (ou minoritaire) au bord de couche limite. `bprime` gère le cas
multi-élément par la formule généralisée
`B'c = (Y_e + B'g·Y_g − Y_w(1+B'g))/(Y_w − Y_c)`.

> Attention à la cohérence : si le char retient des atomes, il faut les
> **retirer du gaz**. C'est automatique dans la fermeture élémentaire, pas si
> l'on recopie deux compositions de sources différentes.

---

## 6. Les commandes

```bash
export MPP_DATA_DIRECTORY=$PWD/data

# table B'c et h_w, à un B'g donné
bprime -T 300:25:5000 -P 101325 -b 0.5 \
       -m cork-air -bl air -py cork_pyro -char cork_char -char-elem C

# enthalpie et propriétés du gaz de pyrolyse
mppequil -T 200:25:4000 -P 101325 -m 0,5,9,10,18,3,32 \
         --elem-comp cork_pyro cork-pyrogas
```

ou, pour les tables complètes et les figures :

```bash
cd cork_bprime
python cork_pyrolysis_data.py   # compositions élémentaires + sensibilités
python cork_bprime.py           # table B'c / h_w, B'c(B'g), point de fonctionnement
python cork_pyrolysis_gas.py    # h_g, M, Cp, gamma, rho, mu du gaz
```

### Vérifications de recette

| Contrôle | Attendu | Obtenu |
|---|---|---|
| `checkmix cork-air` | 26 espèces, 4 éléments, C(gr) en `solid` | ✔ |
| B'c à 300 K, B'g = 0, 1 atm | ≈ 0.087 (limite C + O₂ → CO₂) | **0.0874** |
| B'c ≈ 200 | signature d'un C(gr) manquant | absent |

---

## 7. D'où vient le rendement en char, et ce qu'il pèse

### 7.1 La mesure : TGA du cork P50

La TGA publiée sur le **P50** (Amorim ; matériau de l'IXV et de QARMAN) donne
un rendement en char **du composite**, sous argon à 10 K/min :

| | valeur |
|---|---|
| début de décomposition | ~430 K (masse 98 %) |
| entièrement carbonisé | 780 K, masse **24.5 %** |
| plateau jusqu'à 1650 K | masse **20 %** ← retenu |
| ρ vierge / ρ char | 464.5 et 466.7 / 279.9 et 298.4 kg/m³ |

Source : Sakraker et al., *Performance of cork-based thermal protection
material P50 exposed to air plasma*, CEAS Space Journal 14:377-393 (2022).

L'argon est essentiel : sous air on mesurerait l'oxydation du char, pas la
pyrolyse. Le plateau est lu à 1650 K, pas à 780 K, pour être cohérent avec la
température de référence du modèle de pyrolyse. Et 10 K/min reste lent devant
une rentrée : une pyrolyse flash donnerait moins de char.

### 7.2 Du composite aux constituants

La TGA ne sépare pas le liège de la résine. Avec la résine à 50 % (énoncé) :

```
0.80 · y_liège + 0.20 · 0.50 = 0.20   =>   y_liège = 12.5 %
```

| rendement composite | liège impliqué | gaz C / H / O | k |
|---|---|---|---|
| **20.0 % (plateau, retenu)** | **12.5 %** | **0.287 / 0.592 / 0.121** | **4.00** |
| 24.5 % (plateau à 780 K) | 18.1 % | 0.266 / 0.609 / 0.124 | 3.08 |
| *29.3 % (additivité des constituants, §3.1)* | *24.2 %* | *0.243 / 0.629 / 0.128* | *2.41* |

Pour lever l'hypothèse « résine à 50 % », il faut une seconde TGA — sur la
résine seule ou sur le liège seul.

### 7.3 Sensibilité

| char liège | C | H | O | char composite | k |
|---|---|---|---|---|---|
| 5 % | 0.313 | 0.571 | 0.117 | 14.0 % | 6.143 |
| 10 % | 0.296 | 0.585 | 0.119 | 18.0 % | 4.556 |
| **12.5 %** | **0.287** | **0.592** | **0.121** | **20.0 %** | **4.000** |
| 15 % | 0.278 | 0.600 | 0.122 | 22.0 % | 3.545 |
| 20 % | 0.259 | 0.615 | 0.126 | 26.0 % | 2.846 |
| 24.2 % *(additivité)* | 0.243 | 0.629 | 0.128 | 29.4 % | 2.406 |

±5 points de rendement déplacent le carbone du gaz de ±0.018 et `k` de ∓25 %.
C'est **le** paramètre à mesurer en priorité, avant l'analyse élémentaire du
liège : `k` gouverne le point de fonctionnement (§8), la composition ne le
déplace que marginalement.

> **Piège : k ne se calcule pas sur les densités pour ce matériau.**
> L'identité `k = (ρ_v − ρ_c)/ρ_c` suppose un **volume constant**. Le P50 ne
> la respecte pas : ρ_c/ρ_v = 0.62 alors que la masse résiduelle est 0.20,
> soit V_char/V_vierge ≈ 0.32 — le char se rétracte. Prendre les densités
> donnerait k = 0.61, **six fois trop peu**. Le couplage stationnaire se
> calcule sur les **masses** : `k = (1 − y)/y`, avec y le rendement en char
> massique. Chez les carbone/phénolique (TACOT, CPh70), qui ne changent pas
> de volume, les deux formules coïncident — d'où la confusion possible.

> **Piège d'énoncé.** « 80 % de liège » : masse ou volume ? Le liège est très
> léger (≈ 120 kg/m³ contre 1200 pour la résine). 80 % en **volume** ne fait
> que **28.6 %** en masse, ce qui donne un tout autre matériau : gaz
> C:0.277, H:0.611, O:0.112, char 39.3 %, k = 1.545. Le présent document
> retient la lecture **massique**.

---

## 8. Conséquence physique : où le matériau lit la table

En ablation stationnaire le matériau impose `B'g = k·B'c` avec
`k = m_gaz/m_char = (1 − y)/y = 4.0` (80 g de gaz pour 20 g de char).
Le point de fonctionnement est le point fixe
`B'c = B'c_table(T, P, B'g = k·B'c)`.

À 1 atm :

| T [K] | B'c (B'g = 0) | **B'c point de fonct.** | B'g associé | h_w [MJ/kg] |
|---|---|---|---|---|
| 1000 | 0.1540 | **0.0598** | 0.239 | −1.11 |
| 2000 | 0.1749 | **0.0769** | 0.307 | 0.94 |
| 3000 | 0.1768 | **0.1138** | 0.455 | 4.41 |
| 3400 | 0.2097 | **0.2213** | 0.885 | 9.57 |

Le soufflage pyrolytique **divise par 2.6** l'ablation du char en régime
d'oxydation (1000 K) : le gaz, riche en H et en C, consomme l'oxygène avant
qu'il n'atteigne la paroi. Effet marginal pour un TACOT (B'g ≈ 0.04),
dominant ici (B'g ≈ 0.24 – 0.9).

Deux conséquences pratiques :

1. **Balayer B'g jusqu'à ~10, pas jusqu'à 2.** Une table s'arrêtant à B'g = 2
   ne couvre pas le point de fonctionnement au-delà de ~3500 K.
2. **Aux hautes températures le point fixe diverge** (emballement de
   sublimation) : les lignes du CSV où `Bg_ss` atteint la borne du balayage
   (10) ou où `Bc_ss` vaut 1000 ne sont pas des solutions physiques mais la
   signature de B'c → ∞. Les filtrer.

> La récession se lit sur ρ_c, pas sur B'c : `ṡ = B'c·ṁe/ρc`. La colonne
> `recession_over_mdote_m3_per_kg` du CSV donne B'c/ρc avec la **ρ_char
> mesurée** (289.1 kg/m³) et non 0.20·ρ_v = 93 kg/m³ : le char se rétracte
> (§7), et prendre l'hypothèse « sans retrait » surestimerait la récession
> d'un facteur trois.

---

## 9. Récapitulatif

> Dans le XML : **une liste d'espèces avec la phase condensée** et **trois
> compositions élémentaires**. Pour un liège/phénolique, la seule vraie
> difficulté est la deuxième : le gaz de pyrolyse est le **mélange** des gaz
> des deux constituants, obtenu par fermeture élémentaire, et il faut donc
> connaître le rendement en char **du liège** autant que celui de la résine —
> ici via le rendement du composite mesuré par TGA (20 %), dont le liège se
> déduit (12.5 %).
> Un deuxième fichier, sans phase condensée ni air, donne h_g.
