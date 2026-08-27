# Table B', h_w et h_g d'un liège/phénolique (cork phenolic)

**Matériau** : 80 % liège / 20 % résine phénolique (fractions **massiques**),
résine novolac C7H6O de rendement en char 50 %, et **rendement en char du
composite 20 %** (TGA du cork P50, argon, 10 K/min) — d'où un rendement du
liège de 12.5 %. Le liège est reconstruit depuis ses constituants (45 %
subérine, 27 % lignine, 12 % polysaccharides, 6 % tanins, 6 % céroïdes) et
leurs unités de répétition : C 66.2 / H 8.7 / O 25.1 % masse.

Particularité par rapport à un TACOT / PICA / Zuram / CPh70 : **le renfort
pyrolyse aussi**. Le gaz de pyrolyse est le mélange du gaz du liège et de
celui de la résine, et le rapport liège/résine entre donc dans la table B'.

| | |
|--|--|
| mélange paroi | `data/mixtures/cork-air.xml` (25 espèces + C(gr)) |
| mélange gaz de pyrolyse | `data/mixtures/cork-pyrogas.xml` (17 espèces C/H/O, sans phase condensée) |
| gaz de pyrolyse (`cork_pyro`) | C:0.300, H:0.594, O:0.107 |
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
| liège | 80 g | 10 g | 70 g | 0.304 / 0.589 / 0.107 |
| résine | 20 g | 10 g | 10 g | 0.269 / 0.626 / 0.104 |
| **composite** | 100 g | 20 g | 80 g | **0.300 / 0.594 / 0.107** |

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
| 1000 | 0.1540 | 0 | 0 | 0 |
| 2000 | 0.1749 | 0 | 0 | 0 |
| 3000 | 0.1768 | 0.0767 | 0 | 0 |
| 3500 | 0.2472 | 0.2977 | 0.2903 | 0.1791 |

Point de fonctionnement stationnaire (`B'c = table(T, P, B'g = 4·B'c)`),
1 atm :

| T [K] | B'c sans pyrolyse | **B'c stationnaire** | B'g | h_w [MJ/kg] |
|---|---|---|---|---|
| 1000 | 0.1540 | **0.0551** | 0.221 | −1.06 |
| 2000 | 0.1749 | **0.0698** | 0.279 | 0.96 |
| 3000 | 0.1768 | **0.1010** | 0.404 | 4.35 |
| 3400 | 0.2097 | **0.1858** | 0.743 | 9.30 |

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
| 300 | −4713.3 | 30.945 | 1.341 |
| 1000 | −1748.3 | 22.390 | 5.121 |
| 2000 | 6435.6 | 14.667 | 3.478 |
| 3000 | 11801.4 | 13.744 | 10.588 |
| 4000 | 34610.1 | 8.858 | 26.609 |

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
