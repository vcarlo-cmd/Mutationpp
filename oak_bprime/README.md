# Table B', h_w et h_g d'un bois de chêne

**Matériau** : bois de chêne massif (*Quercus robur* / *Q. petraea*), **sec**,
sans résine, **rendement en char 20 %** (donnée de l'énoncé). Analyse
élémentaire : C 50.0 / H 6.1 / O 43.9 % masse (littérature, sec sans cendres —
**provenance à confirmer**, cf. `mise_en_donnees_chene.md` §2) ; la
reconstruction depuis les constituants (43 % cellulose, 22 % hémicelluloses,
25 % lignine S+G, 8 % ellagitanins) donne C 50.4 / H 6.0 / O 43.6 et sert de
contrôle. Unité de répétition globale : **C₆H₈.₇O₃.₉₅** (CH₁.₄₅O₀.₆₆).

Particularité : c'est le plus simple des ablateurs pyrolysants du dépôt — **un
seul constituant**, qui pyrolyse entièrement. Pas de renfort inerte (TACOT),
pas de mélange de deux gaz (liège/phénolique).

| | |
|--|--|
| mélange paroi | `data/mixtures/oak-air.xml` (25 espèces + C(gr)) |
| mélange gaz de pyrolyse | `data/mixtures/oak-pyrogas.xml` (17 espèces C/H/O, sans phase condensée) |
| gaz de pyrolyse (`oak_pyro`) | C:0.221, H:0.536, O:0.243 |
| arête de couche limite (`air`) | N:0.79, O:0.21 |
| char (`oak_char`) | C:1.0 |
| bilan sur 100 g de bois sec | 20 g de char + 80 g de gaz |
| couplage stationnaire | **k = B'g/B'c = (1−y)/y = 4.0** (TACOT : 0.273) |
| ρ vierge / ρ char | 700 (ordre de grandeur) / 255 kg/m³ (**hypothèse** V_char/V_v = 0.55, à mesurer) |

**L'explication complète de la mise en données est dans
[`mise_en_donnees_chene.md`](mise_en_donnees_chene.md).** L'onglet `Chêne` du
classeur `../mise_en_donnees_xlsx/mise_en_donnees_materiaux.xlsx` refait la
chaîne en formules vivantes.

```bash
export MPP_DATA_DIRECTORY=$PWD/../data
python oak_pyrolysis_data.py   # compositions élémentaires + sensibilités
python oak_bprime.py           # table B'c / h_w, B'c(B'g), point de fonctionnement
python oak_pyrolysis_gas.py    # h_g, M, Cp, gamma, rho, mu du gaz de pyrolyse
```

---

## 1. `oak_pyrolysis_data.py` — les compositions élémentaires

Fermeture élémentaire `n_E(gaz) = n_E(bois) − n_E(char)` sur 100 g de bois sec :

| | C | H | O |
|---|---|---|---|
| bois (mol d'atomes) | 4.163 | 6.052 | 2.744 |
| char, C pur (20 g) | 1.665 | 0 | 0 |
| **gaz** (80 g) | 2.498 | 6.052 | 2.744 |
| **gaz, x molaires** | **0.221** | **0.536** | **0.243** |

Le script imprime aussi les sensibilités (rendement en char, humidité du bois)
et la variante « char non purement carboné ».

## 2. `oak_bprime.py` — table B'c et h_w

25 isobares de 10⁻³ à 10³ atm, 189 températures de 300 à 5000 K, pour sept
valeurs de B'g : 0, 0.1, 0.2, 0.5, 1.0, 2.0 et 5.0.

À 1 atm, effet du soufflage pyrolytique sur l'ablation du char :

| T [K] | B'g = 0 | B'g = 0.5 | B'g = 2 | B'g = 5 |
|---|---|---|---|---|
| 300 | 0.0874 | 0 | 0 | 0 |
| 1000 | 0.1540 | 0.1279 | 0.0559 | 0 |
| 2000 | 0.1749 | 0.1986 | 0.2604 | 0.3801 |
| 3000 | 0.1768 | 0.2740 | 0.4632 | 0.7906 |
| 3500 | 0.2472 | 0.4725 | 1.0184 | 2.0259 |

**Au-dessus de ~1500 K, B'c croît avec B'g** : le gaz du bois porte plus
d'oxygène que de carbone (O/C = 1.10), il oxyde le char au lieu de le
protéger. C'est l'inverse du liège/phénolique (O/C = 0.42), pour lequel le
soufflage fait tomber B'c à zéro jusqu'à 3000 K.

Point de fonctionnement stationnaire (`B'c = table(T, P, B'g = 4·B'c)`),
1 atm :

| T [K] | B'c sans pyrolyse | **B'c stationnaire** | B'g | h_w [MJ/kg] |
|---|---|---|---|---|
| 1000 | 0.1540 | **0.1274** | 0.510 | −1.77 |
| 2000 | 0.1749 | **0.2135** | 0.854 | 0.49 |
| 3000 | 0.1768 | **0.4305** | 1.722 | 4.30 |

À 3400 K, `Bg_ss` atteint la borne du balayage (10) : ce point est hors
table. Filtrer ces lignes du CSV, comme pour le liège.

> **Plafond numérique.** `Thermodynamics::surfaceMassBalance` ajoute une
> quantité finie de char, `max(100·B'g, 200)` ; au-delà de la sublimation
> complète, B'c sature exactement sur cette valeur. C'est la signature de B'c → ∞, pas
> une solution physique.

## 3. `oak_pyrolysis_gas.py` — enthalpie du gaz de pyrolyse

Gaz de pyrolyse pur, à l'équilibre chimique, **sans phase condensée**. 25
isobares de 10⁻³ à 10³ atm, T de 200 à 4000 K. À 1 atm :

| T [K] | h_g [kJ/kg] | M [kg/kmol] | Cp [kJ/kg/K] |
|---|---|---|---|
| 300 | −8197.4 | 27.720 | 1.301 |
| 1000 | −3406.4 | 15.640 | 11.440 |
| 2000 | −358.8 | 14.475 | 2.618 |
| 3000 | 3942.1 | 13.726 | 7.927 |
| 4000 | 18415.6 | 10.201 | 14.865 |

h_g à 300 K est bien plus négatif que pour le liège (−5257 kJ/kg) : le gaz
est dominé par H₂O et CO₂, très stables.

---

## Fichiers produits

| Fichier | Contenu |
|---------|---------|
| `oak_pyrolysis_data.csv` | sensibilité de la composition au rendement en char |
| `oak_bprime_Bg{0p0,0p1,0p2,0p5,1p0,2p0,5p0}.csv` | table B' complète : 25 isobares × 189 T × (B'c, h_w, 26 fractions molaires) |
| `oak_bprime_Bg*.png` | isobares B'c et h_w pour chaque B'g |
| `oak_bprime_bc_table.csv` | B'c au format long, SI : `Tw_K, P_bar, Bg, Bc` |
| `oak_bprime_hw_table.csv` | h_w à B'g = 0, SI : `Tw_K, P_bar, hw_Jkg` |
| `oak_bprime_bg_comparison.png` | influence de B'g à 1 atm |
| `oak_bprime_Bc_vs_Bg.csv` / `.png` | B'c(B'g) à T fixée, avec les points de fonctionnement |
| `oak_bprime_steady_state.csv` | point de fonctionnement stationnaire : `Bc_ss`, `Bg_ss`, `hw_ss`, `Bc_Bg0`, `B'c/ρc` |
| `oak_pyrolysis_gas.csv` / `.png` | h_g, M, Cp, γ, ρ, μ du gaz de pyrolyse |
| `oak_pyrolysis_gas_enthalpy.png` | zoom sur h_g et sa sensibilité à la pression |

Voir `../cork_bprime/` (liège/phénolique, même trame) et `../tacot_bprime/`
(mécanique générale du XML).
