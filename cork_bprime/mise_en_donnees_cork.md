# Mise en données XML d'un liège/phénolique (cork phenolic)

**Matériau traité** : 80 % liège / 20 % résine phénolique (fractions
**massiques** du solide), résine de rendement en char **50 %**, liège de
rendement en char **25 %** (valeur choisie, cf. §7).

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
| rendement en char du renfort | sans objet (0 % de gaz) | **donnée d'entrée indispensable** |
| couplage stationnaire k = B'g/B'c | 0.27 (TACOT) | **2.33** — dix fois plus de soufflage |

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

| Constituant | Composition élémentaire | Rendement en char |
|---|---|---|
| liège | C 62.4 / H 8.5 / O 28.4 (% masse, sec sans cendres) | 25 % *(choisi)* |
| résine phénolique | novolac **C7H6O** (M = 106.12 g/mol) | 50 % *(énoncé)* |

L'analyse élémentaire du liège est une donnée de **litt. matériau** (subérine,
lignine, polysaccharides) ; l'azote (~0.6 % masse) est négligé — le retenir
ajouterait N au gaz de pyrolyse, donc une quatrième composante à `cork_pyro`.

### 3.2 Bilan sur 100 g de composite vierge

| | masse | char | gaz | moles d'atomes du gaz |
|---|---|---|---|---|
| liège | 80 g | 20 g | 60 g | C 2.520, H 6.794, O 1.430 |
| résine | 20 g | 10 g | 10 g | C 0.487, H 1.131, O 0.188 |
| **composite** | 100 g | **30 g** | **70 g** | **C 3.007, H 7.924, O 1.619** |

Fractions molaires élémentaires :

| | C | H | O |
|---|---|---|---|
| gaz du liège seul | 0.235 | 0.632 | 0.133 |
| gaz de la résine seule | 0.269 | 0.626 | 0.104 |
| **gaz du composite** | **0.240** | **0.631** | **0.129** |
| *pour mémoire, TACOT* | *0.206* | *0.679* | *0.115* |

Soit en masse : C 0.516, H 0.114, O 0.370 (le parseur convertit tout seul,
cf. `mise_en_donnees_xml.md` §3.2 — mais il faut déclarer le bon `type`).

Calcul, sensibilités et variantes : `python cork_pyrolysis_data.py`.

### 3.3 Ce que le rapport 80/20 change vraiment

Les deux gaz se ressemblent (le liège est un peu plus oxygéné, la résine un
peu plus carbonée), donc **le mélange déplace peu la composition** :

| liège (% masse) | C | H | O | k = B'g/B'c |
|---|---|---|---|---|
| 0 % (résine pure) | 0.269 | 0.626 | 0.104 | 1.000 |
| 50 % | 0.249 | 0.630 | 0.122 | 1.667 |
| **80 %** | **0.240** | **0.631** | **0.129** | **2.333** |
| 100 % (liège pur) | 0.235 | 0.632 | 0.133 | 3.000 |

Conclusion honnête : le rapport liège/résine entre bien dans la table B'
(contrairement au cas carbone/phénolique), mais son effet **direct** sur la
composition du gaz est modeste. Son effet **dominant** est sur la *quantité*
de gaz, c.-à-d. sur `k` — qui triple entre la résine pure et le liège pur.

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
        <composition name="cork_pyro">C:0.240, H:0.631, O:0.129</composition>
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
- Rien d'autre : ni densité (470 kg/m³), ni porosité, ni rendement en char.
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
        <composition name="cork_pyro">C:0.240, H:0.631, O:0.129</composition>
    </element_compositions>
</mixture>
```

Deux points d'attention :

- **pas de phase condensée** : c'est la convention des tables de propriétés de
  gaz de pyrolyse (classeur TACOT 3.0, feuille *Pyrolysis model*). Ajouter
  C(gr) ferait précipiter du carbone solide et changerait h_g. Ce gaz étant
  plus carboné que celui d'un carbone/phénolique (C:0.240 vs C:0.206), la
  question du dépôt de suie mérite d'être posée si l'on vise un bilan
  d'énergie fin.
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
| char | 0.824 | 0.140 | 0.035 |
| gaz (recalculé, les atomes retenus manquent au gaz) | 0.260 | 0.616 | 0.124 |

```xml
<composition name="cork_char">C:0.824, H:0.140, O:0.035</composition>
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

## 7. Le rendement en char du liège : le paramètre choisi

**25 %** retenu ici. C'est un choix, à assumer comme tel : le liège
(subérine + lignine + polysaccharides) donne en ATG un résidu carboné de
l'ordre de 20 à 30 % à 1000 K, sensible à la vitesse de chauffe et à la
pression. Sensibilité (script §2) :

| char liège | C | H | O | char composite | k |
|---|---|---|---|---|---|
| 15 % | 0.278 | 0.600 | 0.122 | 22.0 % | 3.545 |
| 20 % | 0.259 | 0.615 | 0.126 | 26.0 % | 2.846 |
| **25 %** | **0.240** | **0.631** | **0.129** | **30.0 %** | **2.333** |
| 30 % | 0.219 | 0.649 | 0.132 | 34.0 % | 1.941 |
| 40 % | 0.174 | 0.686 | 0.140 | 42.0 % | 1.381 |

±5 points de rendement déplacent le carbone du gaz de ±0.02 (≈ 9 %) et `k` de
∓17 %. C'est **le** paramètre à mesurer en priorité (ATG), avant l'analyse
élémentaire du liège.

> **Piège d'énoncé.** « 80 % de liège » : masse ou volume ? Le liège est très
> léger (≈ 120 kg/m³ contre 1200 pour la résine). 80 % en **volume** ne fait
> que **28.6 %** en masse, ce qui donne un tout autre matériau : gaz
> C:0.256, H:0.628, O:0.115, char 42.9 %, **k = 1.333**. Le présent document
> retient la lecture **massique**.

---

## 8. Conséquence physique : où le matériau lit la table

En ablation stationnaire le matériau impose `B'g = k·B'c` avec
`k = m_gaz/m_char = (ρ_v − ρ_c)/ρ_c = 2.333` (70 g de gaz pour 30 g de char).
Le point de fonctionnement est le point fixe
`B'c = B'c_table(T, P, B'g = k·B'c)`.

À 1 atm :

| T [K] | B'c (B'g = 0) | **B'c point de fonct.** | B'g associé | h_w [MJ/kg] |
|---|---|---|---|---|
| 1000 | 0.1540 | **0.0896** | 0.209 | −1.11 |
| 2000 | 0.1749 | **0.1152** | 0.269 | 0.94 |
| 3000 | 0.1768 | **0.1705** | 0.398 | 4.41 |
| 3400 | 0.2097 | **0.3309** | 0.772 | 9.55 |

Le soufflage pyrolytique **divise par près de deux** l'ablation du char en
régime d'oxydation (1000 K) : le gaz, riche en H et en C, consomme l'oxygène
avant qu'il n'atteigne la paroi. Effet marginal pour un TACOT (B'g ≈ 0.04),
dominant ici (B'g ≈ 0.2 – 0.8).

Deux conséquences pratiques :

1. **Balayer B'g jusqu'à ~10, pas jusqu'à 2.** Une table s'arrêtant à B'g = 2
   ne couvre pas le point de fonctionnement au-delà de ~3500 K.
2. **Aux hautes températures le point fixe diverge** (emballement de
   sublimation) : les lignes du CSV où `Bg_ss` atteint la borne du balayage
   (10) ou où `Bc_ss` vaut 1000 ne sont pas des solutions physiques mais la
   signature de B'c → ∞. Les filtrer.

> La récession se lit sur ρ_c, pas sur B'c : `ṡ = B'c·ṁe/ρc`. Avec
> ρ_c ≈ 141 kg/m³ (30 % de 470, sans retrait), le liège/phénolique recule
> ~9 fois plus vite qu'un CPh70 à B'c égal. Colonne
> `recession_over_mdote_m3_per_kg` du CSV.

---

## 9. Récapitulatif

> Dans le XML : **une liste d'espèces avec la phase condensée** et **trois
> compositions élémentaires**. Pour un liège/phénolique, la seule vraie
> difficulté est la deuxième : le gaz de pyrolyse est le **mélange** des gaz
> des deux constituants, obtenu par fermeture élémentaire, et il faut donc
> connaître le rendement en char **du liège** autant que celui de la résine.
> Un deuxième fichier, sans phase condensée ni air, donne h_g.
