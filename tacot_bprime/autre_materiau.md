# Changer la composition fibres/résine : qu'est-ce qui change pour la table B' ?

Cas traité : **CPh70** — 70 % fibres / 30 % résine (en volume solide),
porosité 0.01, fibres et résine identiques à celles du TACOT.

---

## 1. Réponse courte

**La table B' ne change pas.** Elle est rigoureusement identique à celle du
TACOT — vérifié numériquement : écart maximal `0.000e+00` sur B'c(T, B'g) aux
quatre pressions.

Ce qui change est la **réponse matériau** : quel B'g le matériau atteint
réellement.

---

## 2. Pourquoi

La table B' est une **condition limite de thermochimie de surface**. Le solveur
`Thermodynamics::surfaceMassBalance`, appelé par `bprime`, ne consomme que
(cf. `src/apps/bprime.cpp:304-330`) :

| Entrée | Origine | Change avec fibres/résine/porosité ? |
|--------|---------|--------------------------------------|
| `Yke` — bord de couche limite | air | non |
| `Ykg` — gaz de pyrolyse | **résine seule** | **non** (voir ci-dessous) |
| `Ychar` — char | fibres + résine carbonisée | **non** (voir ci-dessous) |
| `T`, `P`, `B'g` | conditions de vol | non |

Ni la densité, ni la porosité, ni la fraction volumique de fibres n'apparaissent
dans ce bilan — ce sont des grandeurs **volumiques**, alors que le bilan de
surface est écrit sur des **fractions massiques élémentaires**, donc normalisé.

Les deux points clés :

- **Gaz de pyrolyse** : il est produit par la résine seule. À résine identique,
  sa composition élémentaire reste `C:0.206, H:0.679, O:0.115`. Changer le
  rapport fibres/résine change la **quantité** de gaz produite, pas sa
  **composition**. Or B'g est justement la quantité, et c'est un paramètre
  d'entrée de la table, balayé de 0 à 10.

- **Char** : les fibres de carbone et la résine carbonisée sont **toutes deux du
  carbone pur**. Donc `Y_c = C:1.0` quelles que soient les proportions. Si les
  deux constituants avaient des compositions élémentaires différentes, il
  faudrait au contraire pondérer — voir §5.

---

## 3. Ce qui change réellement : le bilan de phase solide

Avec ρ_fibre = 1600 kg/m³, ρ_résine = 1200 kg/m³, rendement en char de la
résine 50 % (valeurs TACOT — origine de chacune dans `resine_tacot.md`) :

| | TACOT | CPh70 |
|--|-------|------|
| fibres / résine / pores (volume) | 0.10 / 0.10 / 0.80 | 0.693 / 0.297 / 0.01 |
| **ρ vierge** | 280.0 kg/m³ | **1465.2 kg/m³** (×5.2) |
| **ρ char** | 220.0 kg/m³ | **1287.0 kg/m³** |
| fibres (% masse, vierge) | 57.1 % | 75.7 % |
| résine (% masse, vierge) | 42.9 % | 24.3 % |
| gaz de pyrolyse libéré | 60.0 kg/m³ (21.4 % masse) | 178.2 kg/m³ (**12.2 %** masse) |
| **couplage k = B'g/B'c** | **0.2727** | **0.1385** (×0.51) |

Le CPh70 est cinq fois plus dense, mais **riche en fibres** : il produit
relativement **deux fois moins** de gaz de pyrolyse par unité de char consommé.

---

## 4. Où la composition entre : la ligne de fonctionnement

En ablation stationnaire (matériau entièrement carbonisé reculant à vitesse
constante), les flux de gaz de pyrolyse et de char sont liés par le matériau :

$$
\frac{B'_g}{B'_c} = \frac{\dot m_g}{\dot m_c}
= \frac{\rho_v - \rho_c}{\rho_c} = k
$$

Le point de fonctionnement à (T, P) est donc la solution du **point fixe**

$$
B'_c = B'_c^{\text{table}}\big(T,\, P,\, B'_g = k\,B'_c\big)
$$

La table est commune ; **c'est la droite B'g = k·B'c qui la traverse
différemment**. Résultat à 1 atm :

| T [K] | TACOT B'c | TACOT B'g | CPh70 B'c | CPh70 B'g | écart |
|-------|-----------|-----------|----------|----------|-------|
| 1000 | 0.14221 | 0.0388 | 0.14774 | 0.0205 | +3.9 % |
| 2000 | 0.16684 | 0.0455 | 0.17105 | 0.0237 | +2.5 % |
| 3000 | 0.18962 | 0.0517 | 0.18564 | 0.0257 | −2.1 % |
| 3400 | 0.24090 | 0.0657 | 0.22659 | 0.0314 | −5.9 % |
| 3600 | 0.42054 | 0.1147 | 0.37210 | 0.0515 | −11.5 % |
| 3800 | 3.34781 | 0.9130 | 1.56790 | 0.2171 | **−53.2 %** |

En régime d'oxydation (T < 3000 K) l'écart reste sous 4 %. Il explose au genou
de sublimation : moins de soufflage pyrolytique ⇒ le char survit à plus haute
température ⇒ le genou du CPh70 est décalé vers la droite (visible sur
`material_response_operating_lines.png`).

> Attention : à densité très différente, un même B'c ne donne **pas** la même
> vitesse de recul. La récession est `ṡ = ṁ_c/ρ_c = B'_c·ṁ_e/ρ_c`. À B'c égal,
> le CPh70 recule **5.85 fois moins vite** que le TACOT (1287/220).

---

## 5. Quand la mise en données changerait vraiment

| Changement | Table B' affectée ? | Quoi modifier |
|------------|--------------------|---------------|
| rapport fibres/résine | **non** | rien |
| porosité | **non** | rien |
| densités intrinsèques | **non** | rien |
| **autre résine** (autre polymère) | **oui** | recalculer `*_pyro` depuis la composition du gaz |
| **fibres non carbonées** (silice, SiC…) | **oui** | `*_char` devient multi-élément, et `-char-elem` change ; ajouter les espèces correspondantes |
| char non purement carboné (résidus H/O) | oui | pondérer `*_char` |
| liste d'espèces | oui | cf. `tacot_bprime_validation.md` §4 |

Pour un char multi-élément, `bprime` gère déjà le cas via `-char` / `-char-elem`
(formule généralisée `B'c = (Y_e + B'g·Y_g − Y_w(1+B'g))/(Y_w − Y_c)`) — c'est
ce que fait `silice_bprime` avec `Si:1.0, O:2.0`.

---

## 6. Fichiers

| Fichier | Contenu |
|---------|---------|
| `data/mixtures/cph70-air_25.xml` | mise en données CPh70 (25 espèces + C(gr)) |
| `tacot_bprime/material_response.py` | bilan solide, couplage k, point fixe, tracé |
| `tacot_bprime/material_response_operating_lines.png` | lignes de fonctionnement TACOT vs CPh70 |

```bash
cd tacot_bprime && python material_response.py
```

Le script `solid_balance(f_solid, porosity, rho_f, rho_m, char_yield)` est
générique : il suffit de changer les arguments dans le dictionnaire `MATERIALS`
pour traiter une autre composition.
