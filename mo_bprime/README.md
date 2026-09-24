# Table B' du molybdène (Mo, TZM)

**Molybdène** nu dans l'air, sans pyrolyse (B'g = 0). La table s'applique au
**TZM** (Mo + 0.5 % Ti, 0.08 % Zr, ~0.02 % C). Ces additions améliorent la
tenue mécanique, pas l'oxydation, qui reste celle du molybdène. Même grille et
même format de sortie que les autres matériaux : 25 isobares de 10⁻³ à
10³ atm, Tw de 300 à 5000 K par pas de 25 K.

| | |
|--|--|
| mélange | `data/mixtures/mo-air.xml` : 18 espèces (air, Mo et ses oxydes gazeux jusqu'à Mo5O15, 5 phases condensées) |
| char | `mo` : Mo:1.0 |
| binaire | `bprime` généralisé : `-char mo -char-elem N` |
| couplage | k = B'g/B'c = 0 |

```bash
cd mo_bprime
export MPP_DATA_DIRECTORY=../data     # si ce n'est pas déjà fait
python mo_bprime.py                    # ~5 s
```

| fichier | contenu |
|--|--|
| `mo_bprime_table.csv`, `_bc_table.csv`, `_hw_table.csv` | **équilibre complet** (borne basse), même format que carbone, silice, SiC |
| `mo_bprime_nu_table.csv`, `_nu_bc_table.csv`, `_nu_hw_table.csv` | **surface nue**, sans oxydes condensés (borne haute, limite de diffusion) |
| `mo_bprime_table.png` | B'c des deux tables (trait plein / tirets) et h_w |

Le solveur de Mutation++ converge sur les 4725 points, sans le traitement
point par point nécessaire au SiC. Le contrôle de stabilité des phases du
script SiC, appliqué à la table d'équilibre, ne signale aucun point instable.

---

## 1. Pourquoi deux tables

À l'équilibre avec un excès de métal, l'oxygène qui atteint la paroi est
fixé en **MoO2 solide**. La table d'équilibre donne donc B'c = 0 jusqu'à la
disparition de ce dioxyde :

| P | 0.001 atm | 0.01 | 0.1 | 1 atm | 10 | 100 | 1000 atm |
|--|--|--|--|--|--|--|--|
| oxydes condensés instables dès | 1875 K | 2025 K | 2175 K | **2350 K** | 2600 K | 2925 K | 4775 K |

Or le molybdène et le TZM s'oxydent de façon **catastrophique** dans l'air
dès ~750–800 °C. La face externe de l'oxyde voit l'air et forme **MoO3**, qui
fond à 795 °C et se sublime aussi vite qu'il se forme. Le MoO2 ne reste qu'en
sous-couche mince, sans rôle protecteur. L'hypothèse « paroi à l'équilibre
avec un réservoir infini de métal » est trop réductrice pour le molybdène
dans cette plage.

La seconde table retire donc les oxydes condensés MoO2(cr), MoO3(cr) et
MoO3(L). Tout l'oxygène qui atteint la paroi repart en oxydes gazeux. On
obtient le **régime limité par la diffusion** :

```
Mo + 3/2 O2 → MoO3(g)      B'c = y_O · M_Mo / (3 M_O) = 0.466
```

La table donne 0.4657 sur tout le plateau, en accord avec la valeur
analytique.

## 2. Quelle table utiliser

| T_w | régime réel | table à retenir |
|--|--|--|
| < ~1000 K | oxydation lente, MoO3 non volatil, limitée par la réaction | aucune : B'c ≈ 0, modèle cinétique linéaire si besoin |
| ~1000 K → disparition de MoO2 | MoO3 volatil, vitesse croissante vers la limite de diffusion | **surface nue** (borne haute, conservative) ; le vrai B'c est entre 0 et 0.466 selon la cinétique |
| au-delà | surface nue réellement, oxydes instables | les deux tables coïncident |
| ≳ 2900 K | Mo liquide (fusion à 2896 K), puis sublimation | idem ; B'c = 200 = plafond du solveur (gazéification totale) |

Dans la plage intermédiaire, la vitesse vraie dépend de la cinétique de
surface et de l'écoulement. Si l'écoulement est rapide et la température
élevée, on s'approche de la limite de diffusion, et la table « surface nue »
est une bonne estimation. Les constantes de la littérature (Gulbransen et
Wysong notamment) ou des essais permettent de raffiner cette zone.

## 3. Choix du suivi sur N

Comme pour le SiC (`../sic_bprime/README.md`), le bilan est suivi sur
l'azote : B'c = y_e,N / y_w,N − 1, soit la masse nette gazéifiée. Pour un
char monoélément, c'est identique au suivi sur Mo tant que seul le métal
est condensé. Dans la table d'équilibre, là où MoO2 solide fixe l'oxygène, le
bilan net est une prise de masse, écrêtée à 0.

## 4. Limites

- Le TZM **s'emploie revêtu** dans l'air (siliciures de type MoSi2). La table
  décrit le métal mis à nu par une défaillance du revêtement, pas le
  revêtement lui-même.
- Pas de données thermodynamiques pour Ti, Zr, C du TZM : ils sont ignorés.
- MoO2(cr) est décrit par NASA-9 jusqu'à 6000 K, sans phase liquide. Sa
  stabilité à très haute pression (> 100 atm) est donc à prendre avec
  prudence.
