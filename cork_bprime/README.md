# Table B', h_w et h_g d'un liège/phénolique (cork phenolic)

**Matériau** : 80 % liège / 20 % résine phénolique (fractions **massiques**),
résine novolac C7H6O de rendement en char 50 %, et **rendement en char du
composite 20 %** (TGA du cork P50, argon, 10 K/min) — d'où un rendement du
liège de 12.5 %. Analyse élémentaire du liège : C 62.4 / H 8.5 / O 28.4
(littérature, *Quercus suber* — **provenance à confirmer**, cf.
`mise_en_donnees_cork.md` §3.1) ; la reconstruction depuis les constituants
(45 % subérine, 27 % lignine, 12 % polysaccharides, 6 % tanins, 6 % céroïdes)
donne C 66.2 / H 8.7 / O 25.1 et sert de contrôle.

Particularité par rapport à un TACOT / PICA / Zuram / CPh70 : **le renfort
pyrolyse aussi**. Le gaz de pyrolyse est le mélange du gaz du liège et de
celui de la résine, et le rapport liège/résine entre donc dans la table B'.

| | |
|--|--|
| mélange paroi | `data/mixtures/cork-air.xml` (25 espèces + C(gr)) |
| mélange gaz de pyrolyse | `data/mixtures/cork-pyrogas.xml` (17 espèces C/H/O, sans phase condensée) |
| gaz de pyrolyse (`cork_pyro`) | C:0.287, H:0.592, O:0.121 |
| arête de couche limite (`air`) | N:0.79, O:0.21 |
| char (`cork_char`) | C:1.0 |
| bilan sur 100 g vierge | 20 g de char + 80 g de gaz |
| couplage stationnaire | **k = B'g/B'c = (1−y)/y = 4.0** (TACOT : 0.273) |
| ρ vierge / ρ char (mesurées, P50) | 465.6 / 289.1 kg/m³ |

**L'explication complète de la mise en données est dans
[`mise_en_donnees_cork.md`](mise_en_donnees_cork.md).**

```bash
export MPP_DATA_DIRECTORY=$PWD/../data
python cork_pyrolysis_data.py   # compositions élémentaires + sensibilités
python cork_bprime.py           # table B'c / h_w, B'c(B'g), point de fonctionnement
python cork_pyrolysis_gas.py    # h_g, M, Cp, gamma, rho, mu du gaz de pyrolyse
```

---

## 1. `cork_pyrolysis_data.py` — les compositions élémentaires

Fermeture élémentaire constituant par constituant
(`n_E(gaz) = n_E(constituant) − n_E(char)`), puis somme des deux gaz au
prorata des masses dégazées.

| | masse | char | gaz | gaz : C / H / O (x molaires) |
|---|---|---|---|---|
| liège | 80 g | 10 g | 70 g | 0.290 / 0.587 / 0.124 |
| résine | 20 g | 10 g | 10 g | 0.269 / 0.626 / 0.104 |
| **composite** | 100 g | 20 g | 80 g | **0.287 / 0.592 / 0.121** |

Le script imprime aussi les sensibilités (rendement en char du liège, rapport
liège/résine), la variante « char de liège non purement carboné » et la
variante « 80 % en volume au lieu de 80 % en masse ».

## 2. `cork_bprime.py` — table B'c et h_w

25 isobares de 10⁻³ à 10³ atm, 189 températures de 300 à 5000 K, pour sept
valeurs de B'g : 0, 0.1, 0.2, 0.5, 1.0, 2.0 et **5.0** — ce matériau
fonctionne à fort soufflage, une table s'arrêtant à B'g = 2 ne suffit pas.

À 1 atm, effet du soufflage pyrolytique sur l'ablation du char :

| T [K] | B'g = 0 | B'g = 0.5 | B'g = 2 | B'g = 5 |
|---|---|---|---|---|
| 1000 | 0.1540 | 0.0000 | 0 | 0 |
| 2000 | 0.1749 | 0.0140 | 0 | 0 |
| 3000 | 0.1768 | 0.1050 | 0 | 0 |
| 3500 | 0.2472 | 0.3255 | 0.4030 | 0.4621 |

Point de fonctionnement stationnaire (`B'c = table(T, P, B'g = 4·B'c)`),
1 atm :

| T [K] | B'c sans pyrolyse | **B'c stationnaire** | B'g | h_w [MJ/kg] |
|---|---|---|---|---|
| 1000 | 0.1540 | **0.0598** | 0.239 | −1.11 |
| 2000 | 0.1749 | **0.0769** | 0.307 | 0.94 |
| 3000 | 0.1768 | **0.1138** | 0.455 | 4.41 |
| 3400 | 0.2097 | **0.2213** | 0.885 | 9.57 |

> **Plafond numérique.** `Thermodynamics::surfaceMassBalance` ajoute une
> quantité finie de char, `max(100·B'g, 200)` ; au-delà de la sublimation
> complète B'c sature exactement sur cette valeur — signature de B'c → ∞, pas
> une solution physique. De même, les lignes du CSV de point de fonctionnement
> où `Bg_ss` atteint la borne du balayage (10) sont hors table : les filtrer.

## 3. `cork_pyrolysis_gas.py` — enthalpie du gaz de pyrolyse

Gaz de pyrolyse pur, à l'équilibre chimique, **sans phase condensée** : ni
air, ni char, ni C(gr) — l'état du gaz avant qu'il n'atteigne la paroi. h_g
ferme le bilan d'énergie de surface (avec h_w) et le terme source de pyrolyse
en profondeur (h_g − h_s). 25 isobares de 10⁻³ à 10³ atm, T de 200 à 4000 K.

À 1 atm :

| T [K] | h_g [kJ/kg] | M [kg/kmol] | Cp [kJ/kg/K] |
|---|---|---|---|
| 300 | −5256.9 | 29.932 | 1.346 |
| 1000 | −2136.7 | 21.182 | 5.284 |
| 2000 | 5571.6 | 14.366 | 3.411 |
| 3000 | 10864.2 | 13.481 | 10.339 |
| 4000 | 32672.7 | 8.892 | 25.458 |

---

## Fichiers produits

| Fichier | Contenu |
|---------|---------|
| `cork_pyrolysis_data.csv` | sensibilité des compositions au rendement en char du liège |
| `cork_bprime_Bg{0p0,0p1,0p2,0p5,1p0,2p0,5p0}.csv` | table B' complète : 25 isobares × 189 T × (B'c, h_w, 26 fractions molaires) |
| `cork_bprime_Bg*.png` | isobares B'c et h_w pour chaque B'g |
| `cork_bprime_bg_comparison.png` | influence de B'g à 1 atm |
| `cork_bprime_Bc_vs_Bg.csv` / `.png` | B'c(B'g) à T fixée, avec les points de fonctionnement |
| `cork_bprime_steady_state.csv` | point de fonctionnement stationnaire : `Bc_ss`, `Bg_ss`, `hw_ss`, `Bc_Bg0`, `B'c/ρc` |
| `cork_pyrolysis_gas.csv` / `.png` | h_g, M, Cp, γ, ρ, μ du gaz de pyrolyse |
| `cork_pyrolysis_gas_enthalpy.png` | zoom sur h_g et sa sensibilité à la pression |

Voir `../tacot_bprime/` (TACOT et mécanique générale du XML),
`../cph70_bprime/` et `../zuram_bprime/` pour les carbone/phénolique.

---

## 4. `cork_bprime_oat.py` — même matériau sous torche oxyacétylénique

Seul le bord de couche limite change (`data/mixtures/cork-oat.xml`, même
gaz de pyrolyse, même char, mêmes 26 espèces) : produits de
C2H2 + r O2 en sortie de buse, **sans air entraîné** (pas d'azote),
C:H:O = 2:2:2r. Trois réglages : r = 1.0 (neutre), **1.3** (réglage OAT
usuel), 2.5 (stœchiométrique).

Ce qui pilote la table est l'**oxygène en excès sur le carbone** de la
flamme, seul capable d'oxyder le char (→ CO) :
plateau B'c = M_C(2r−2)/(M_C2H2 + r·M_O2) = 0 / 0.107 / 0.340 (air : 0.175).

1 atm :

| T [K] | 1000 | 1500 | 2000 | 2500 | 3000 | 3400 |
|---|---|---|---|---|---|---|
| B'c, B'g=0 — air | 0.154 | 0.175 | 0.175 | 0.175 | 0.177 | 0.210 |
| B'c, B'g=0 — OAT r=1.0 | 0 | 0 | 0.0003 | 0.0044 | 0.024 | 0.092 |
| B'c, B'g=0 — OAT r=1.3 | 0.001 | 0.106 | 0.107 | 0.110 | 0.127 | 0.190 |
| B'c, B'g=0 — OAT r=2.5 | 0.210 | 0.339 | 0.340 | 0.342 | 0.353 | 0.405 |
| **B'c stationnaire (B'g=4B'c)** — air | 0.060 | 0.075 | 0.077 | 0.086 | 0.114 | 0.221 |
| — OAT r=1.0 | 0 | 0 | 0.0003 | 0.0023 | 0.012 | 0.064 |
| — OAT r=1.3 | 0.001 | 0.046 | 0.046 | 0.048 | 0.062 | 0.132 |
| — OAT r=2.5 | 0.083 | 0.146 | 0.146 | 0.150 | 0.172 | 0.281 |

- Le seuil de sublimation (~3500–4000 K à 1 atm) ne bouge pas : il dépend
  du char, pas du bord.
- h_w est **plus bas** qu'à l'air à basse T (−10 MJ/kg à 300 K contre
  −2.6) : les éléments C/H de la flamme se recombinent en CH4/H2O/CO2 à la
  paroi. h_w et h_e (enthalpie de récupération de la flamme) doivent être
  pris dans la **même référence** (NASA-9, enthalpies de formation
  incluses) — ne pas mélanger avec une h_e « air » ou sensible seule.

### 4.1 Flamme + air entraîné

Loin de la buse, le jet entraîne de l'air ambiant. Cas traités : fraction
**massique** f = 0.25 et 0.50 d'air (N2/O2 = 79/21) mélangée aux produits
de flamme, pour r = 1.0 et 1.3 (compositions `oat_r*_air25/50` de
`cork-oat.xml`). L'azote entre alors au bord de couche limite.

Le palier d'oxydation est la moyenne massique des deux :
B'c = (1 − f)·B'c(flamme) + f·0.175, retrouvé exactement par Mutation++.

1 atm :

| T [K] | 1000 | 1500 | 2000 | 2500 | 3000 | 3400 |
|---|---|---|---|---|---|---|
| B'c, B'g=0 — r=1.0 + 25 % air | 0 | 0.044 | 0.046 | 0.056 | 0.084 | 0.151 |
| — r=1.0 + 50 % air | 0.033 | 0.087 | 0.090 | 0.100 | 0.125 | 0.183 |
| — r=1.3 + 25 % air | 0.042 | 0.124 | 0.126 | 0.135 | 0.160 | 0.222 |
| — r=1.3 + 50 % air | 0.081 | 0.141 | 0.143 | 0.152 | 0.175 | 0.229 |
| **B'c stationnaire** — r=1.0 + 25 % air | 0 | 0.019 | 0.020 | 0.025 | 0.044 | 0.118 |
| — r=1.0 + 50 % air | 0.013 | 0.038 | 0.039 | 0.046 | 0.069 | 0.156 |
| — r=1.3 + 25 % air | 0.016 | 0.053 | 0.054 | 0.061 | 0.083 | 0.172 |
| — r=1.3 + 50 % air | 0.032 | 0.060 | 0.062 | 0.070 | 0.095 | 0.193 |

L'entraînement rapproche la table de celle de l'air sans l'atteindre tant
que r < 2.3 (au-delà, la flamme seule oxyde déjà plus que l'air). Pour une
flamme neutre il change la nature du régime : sans air le char ne s'oxyde
pas avant ~2500 K, avec 25 % d'air il a déjà un palier à 0.044. La
distance buse–éprouvette pèse donc autant que le réglage O2/C2H2.

- Fichiers : `cork_oat_bprime_bc_table.csv` / `_hw_table.csv` (format long
  avec colonnes `O2_C2H2` et `f_air`), `cork_oat_bprime_steady_state.csv`,
  `cork_oat_bprime_vs_air.png`, `cork_oat_bprime_steady_state.png`,
  `cork_oat_air_bprime.png` (effet de l'air entraîné).
