# Mise en données XML d'un bois de chêne

**Matériau traité** : bois de chêne massif (*Quercus robur* / *Q. petraea*),
**sec**, sans résine ni renfort, **rendement en char 20 %** en masse (donnée
de l'énoncé).

Objectif : produire les trois données que consomment `bprime` et `mppequil`
pour la table **B'c(T, P, B'g)**, l'enthalpie de paroi **h_w** et l'enthalpie
du gaz de pyrolyse **h_g**.

Pour la mécanique générale du fichier XML (balises, parseur, liste d'espèces,
pièges), voir `../tacot_bprime/mise_en_donnees_xml.md`. Le présent document
suit la trame de `../cork_bprime/mise_en_donnees_cork.md`.

---

## 1. Ce qui change par rapport aux autres matériaux

| | carbone/phénolique | liège/phénolique | **chêne** |
|--|--|--|--|
| constituants qui pyrolysent | la résine seule | liège **et** résine | **le bois entier** |
| gaz de pyrolyse | gaz de la résine | mélange de deux gaz | **un seul gaz : bois − char** |
| rapport de mélange dans la table B' | non | oui | **sans objet** |
| couplage stationnaire k = B'g/B'c | 0.27 (TACOT) | 4.0 | **4.0** |
| O/C du gaz | 0.56 | 0.42 | **1.10** |

C'est le cas le plus simple des ablateurs pyrolysants : une seule fermeture
élémentaire. La dernière ligne compte le plus : le gaz du bois porte **plus
d'oxygène que de carbone**, et il oxyde le char au lieu de le protéger (§6).

---

## 2. Données d'entrée

### 2.1 Analyse élémentaire retenue

```
chêne :  C 50.0   H 6.1   O 43.9   (% masse, sec sans cendres)
```

Valeurs usuellement citées pour le bois de chêne. L'azote (0.1 à 0.3 %) est
négligé.

> **Provenance à confirmer.** Ces trois nombres sont des valeurs usuelles, pas
> une donnée vérifiée sur une source primaire depuis l'environnement de
> calcul. Ils sont à confronter à une analyse (publiée ou mesurée) avant tout
> usage engageant. Le contrôle biochimique ci-dessous les retrouve à
> 0.4 point près.

### 2.2 Contrôle : le chêne reconstruit depuis ses unités de répétition

| Constituant | % masse | Unité retenue | Formule | M | C | H | O |
|---|---|---|---|---|---|---|---|
| cellulose | 43 | anhydroglucose | C6H10O5 | 162.1 | 44.4 | 6.2 | 49.3 |
| hémicelluloses | 22 | anhydroxylose (glucuronoxylane des feuillus) | C5H8O4 | 132.1 | 45.5 | 6.1 | 48.4 |
| lignine S | 16 | alcool sinapylique (syringyle) | C11H14O4 | 210.2 | 62.8 | 6.7 | 30.4 |
| lignine G | 9 | alcool coniférylique (guaïacyle) | C10H12O3 | 180.2 | 66.7 | 6.7 | 26.6 |
| tanins | 8 | vescalagine (ellagitanin du bois de cœur) | C41H26O26 | 934.6 | 52.7 | 2.8 | 44.5 |

La lignine de feuillu est répartie en 16 % S + 9 % G en masse, soit
S/G = 1.52 en moles. Les parts somment à **98 %** : le complément (cendres
~0.4 %, groupes acétyles, extractibles mineurs) n'est pas de la matière C/H/O
identifiée. On renormalise sur les 98 %.

```
chêne reconstruit  : C 50.39   H 6.04   O 43.57   (% masse)
```

| | C | H | O |
|---|---|---|---|
| **littérature — retenue** | **50.00** | **6.10** | **43.90** |
| reconstruit | 50.39 | 6.04 | 43.57 |
| **écart** | **+0.4** | **−0.1** | **−0.3** |

L'accord est bien meilleur que pour le liège (4 points sur C) :
l'anhydroglucose, qui pèse les deux tiers du bois avec le xylose, est une
unité exacte, et il n'y a pas de subérine dont l'unité serait incertaine.
L'effet sur le gaz est négligeable (x_C = 0.225 au lieu de 0.221).

### 2.3 Unité de répétition globale

Ramenée à 6 atomes de carbone :

```
C6 H8.72 O3.95     (M = 144.1 g/mol)     soit  CH1.454 O0.659
```

C'est l'anhydroglucose C6H10O5 « enrichi en carbone » par la lignine. Cette
unité n'entre pas dans le XML : c'est la même information que l'analyse
élémentaire, écrite en motif.

### 2.4 Rendement en char

| | valeur | origine |
|---|---|---|
| bois sec | **20 %** | énoncé |
| contrôle par additivité | 22.0 % | cellulose 10 %, hémicelluloses 20 %, lignine S 35 %, G 45 %, tanins 40 % (ordres de grandeur d'ATG lente) |

Les deux voies sont cohérentes. Le rendement réel dépend fortement de la
vitesse de chauffe : une pyrolyse flash descend vers 10-15 %.

---

## 3. Le gaz de pyrolyse par fermeture élémentaire

La pyrolyse ne crée ni ne détruit d'atomes :

```
n_E(gaz) = n_E(bois) − n_E(char)
```

Sur 100 g de bois sec :

| moles d'atomes | C | H | O | total |
|---|---|---|---|---|
| bois | 4.1629 | 6.0516 | 2.7439 | 12.958 |
| char (20 g de C) | 1.6651 | 0 | 0 | 1.665 |
| **gaz (80 g)** | **2.4977** | **6.0516** | **2.7439** | 11.293 |

```
x molaires   : C:0.221, H:0.536, O:0.243      <- oak_pyro
y massiques  : C:0.375, H:0.076, O:0.549
H/O = 2.2    (phénoliques 5.9-6.8, liège 4.9)
```

Le char emporte 40 % du carbone du bois et rien d'autre : tout l'hydrogène et
tout l'oxygène passent dans le gaz.

## 4. Le char

Carbone pur, `oak_char` = C:1.0, `-char-elem C`.

Variante : un charbon de bois à 90 % C / 2 % H / 8 % O en masse donne un char
multi-élément C:0.751, H:0.199, O:0.050. Le gaz devient alors
C:0.243, H:0.516, O:0.241. Il faut **retirer** du gaz les atomes gardés par le
char, sans quoi le bilan n'est plus fermé.

---

## 5. Sensibilités

### Au rendement en char

| y | C | H | O | k |
|---|---|---|---|---|
| 10 % | 0.275 | 0.499 | 0.226 | 9.00 |
| 15 % | 0.249 | 0.517 | 0.234 | 5.67 |
| **20 %** | **0.221** | **0.536** | **0.243** | **4.00** |
| 25 % | 0.191 | 0.556 | 0.252 | 3.00 |
| 30 % | 0.159 | 0.579 | 0.262 | 2.33 |

### À l'humidité du bois

L'eau part entièrement dans le gaz :

| humidité | C | H | O | char / bois humide | k |
|---|---|---|---|---|---|
| **0 % (retenu)** | **0.221** | **0.536** | **0.243** | 20.0 % | 4.00 |
| 8 % | 0.196 | 0.551 | 0.253 | 18.4 % | 4.44 |
| 12 % (bois stabilisé à l'air) | 0.184 | 0.558 | 0.258 | 17.6 % | 4.68 |

Un bois à 12 % d'humidité déplace plus la composition que ±5 points de
rendement. L'énoncé doit préciser « sec » ou le taux d'humidité.

---

## 6. Réponse matériau et table B'

```
k = B'g/B'c = m_gaz/m_char = (1 − y)/y = 80/20 = 4.0
```

Comme pour le liège, k se calcule sur les **masses**. Le bois se rétracte en
carbonisant ; avec l'hypothèse V_char/V_v = 0.55, ρ_c = 0.20 × 700 / 0.55 =
255 kg/m³. L'identité (ρ_v − ρ_c)/ρ_c donnerait 1.75, plus de deux fois trop
peu. **ρ_c n'est pas mesurée ici : c'est l'hypothèse à remplacer en premier**
pour la récession `ṡ = B'c·ṁe/ρ_c`.

Effet du soufflage à 1 atm (`oak_bprime.py`) :

| T [K] | B'g = 0 | B'g = 0.5 | B'g = 2 | B'g = 5 |
|---|---|---|---|---|
| 1000 | 0.1540 | 0.1279 | 0.0559 | 0 |
| 2000 | 0.1749 | 0.1986 | 0.2604 | 0.3801 |
| 3000 | 0.1768 | 0.2740 | 0.4632 | 0.7906 |

À basse température, le gaz dilue l'oxygène de l'air et B'c baisse. Au-dessus
de ~1500 K, l'oxygène excédentaire du gaz (O/C = 1.10, avec H2O et CO2) oxyde
le char : **B'c croît avec B'g**. Pour le liège, le soufflage protégeait le
char ; pour le bois, il l'attaque.

Point de fonctionnement stationnaire à 1 atm : B'c = 0.127 / 0.213 / 0.430 à
1000 / 2000 / 3000 K, soit de 0.8 à 2.4 fois la valeur sans pyrolyse.

Contrôle de la physique : à 300 K et B'g = 0, B'c = 0.0874. C'est la limite
C + O2 → CO2, identique aux autres chars carbonés : C(gr) est bien présent.

---

## 7. Fichiers

| Fichier | Rôle |
|---|---|
| `data/mixtures/oak-air.xml` | mélange paroi : `air`, `oak_pyro`, `oak_char`, 25 espèces + C(gr) |
| `data/mixtures/oak-pyrogas.xml` | gaz de pyrolyse seul (h_g), sans phase condensée |
| `oak_pyrolysis_data.py` | fermeture élémentaire, contrôles, sensibilités |
| `oak_bprime.py` | tables B'c / h_w, B'c(B'g), point de fonctionnement |
| `oak_pyrolysis_gas.py` | h_g, M, Cp, γ, ρ, μ du gaz de pyrolyse |
| onglet `Chêne` de `../mise_en_donnees_xlsx/mise_en_donnees_materiaux.xlsx` | la même chaîne en formules vivantes |
