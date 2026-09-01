# Table B' du MX-4926

**MX-4926** — carbone/phénolique chargé, matériau de tuyère de moteur à
propergol solide. Renfort tissé en fibres de carbone ex-rayon (satin de 5 ou
de 8), matrice résol SC-1008, charge de noir de carbone.

| | |
|--|--|
| composition retenue | 52.2 % fibres / 34.2 % résine / 13.6 % charge (masse) |
| ρ vierge | 1451 kg/m³ *(porosité 0.02)* |
| ρ char | 1228 kg/m³ |
| gaz de pyrolyse libéré | 223 kg/m³ (15.4 % de la masse vierge) |
| rendement en char du composite | Y_comp = 0.8463 |
| couplage stationnaire | **k = B'g/B'c = 0.1817** *(enveloppe 0.1621 – 0.1998)* |
| mélange | `data/mixtures/mx4926-air.xml` (25 espèces + phénol + C(gr)) |
| gaz de pyrolyse | C:0.2526, H:0.6407, O:0.1068 — **hérité du SC-1008** |
| char | C:1.0 |

```bash
cd mx4926_bprime
python composition_mx4926_verification.py    # les 6 contrôles, ~2 s
python mx4926_bprime.py                      # tables et figures, ~1 min
```

> Traçabilité complète — fermeture du domaine de composition, héritage du gaz,
> rendement en char, densités : [`composition_mx4926.md`](composition_mx4926.md).
>
> La résine elle-même est documentée dans
> [`../sc1008_bprime/resine_sc1008.md`](../sc1008_bprime/resine_sc1008.md).

---

## « La table B' est-elle différente de celle du SC-1008 ? »

Il faut distinguer deux objets — c'est la même situation que le couple
CPh70 / TACOT.

**1. La table comme fonction tabulée B'c(T, P, B'g) — identique.**
Pas « proche » : **identique caractère pour caractère**, vérifié à B'g =
0, 0.1, 0.5 et 2.0. Le bilan de masse de surface ne voit que trois
compositions élémentaires, et les trois sont celles du SC-1008 :

- le **gaz** est produit par la résine seule — le renfort ex-rayon est déjà
  carbonisé, le noir de carbone est inerte : ni l'un ni l'autre ne dégaze.
  La matrice étant le SC-1008, le gaz est celui du SC-1008 ;
- le **char** vaut `C:1.0` : les trois constituants carbonisent vers le
  carbone pur.

Ni la densité, ni la porosité, ni le taux de charge, ni le rapport
fibres/résine n'entrent dans ce bilan.

**2. Le B'c réellement atteint — différent.**
En ablation stationnaire le matériau impose

$$B'_g = k\,B'_c, \qquad k = \frac{\rho_v - \rho_c}{\rho_c} = \frac{1 - Y_{comp}}{Y_{comp}}$$

Le point de fonctionnement est l'intersection de la courbe B'c(B'g) — commune
— avec cette droite, propre à chaque matériau. C'est ce que montre
`mx4926_bprime_Bc_vs_Bg.png` : les courbes colorées sont communes, les ★
(MX-4926), ◆ (PICA) et ■ (TACOT) sont les intersections.

| Matériau | k | |
|---|---|---|
| TACOT (10/10/80 vol.) | 0.2727 | poreux |
| PICA (SC-1008) | 0.2346 | même résine, sans charge |
| **MX-4926** | **0.1817** | dense et chargé |
| CPh70 (69/30/1 vol.) | 0.1385 | dense, riche en fibres |

Le noir de carbone agit comme des fibres du point de vue du bilan : il ajoute
de la masse qui ne dégaze pas, donc il **abaisse** k.

### Points de fonctionnement

| P [atm] | T [K] | MX-4926 B'c | MX-4926 B'g | PICA B'c | TACOT B'c | MX / TACOT |
|---|---|---|---|---|---|---|
| 0.01 | 2000 | 0.16645 | 0.0302 | 0.16399 | 0.16223 | +2.60 % |
| 0.01 | 3000 | 0.24789 | 0.0450 | 0.24811 | 0.24827 | −0.15 % |
| 1 | 2000 | 0.16647 | 0.0302 | 0.16401 | 0.16226 | +2.60 % |
| 1 | 3000 | 0.18187 | 0.0330 | 0.18214 | 0.18198 | −0.06 % |
| 1 | 3400 | 0.22321 | 0.0406 | 0.22607 | 0.22818 | −2.18 % |
| 1 | 3700 | 0.62848 | 0.1142 | 0.67160 | 0.70522 | −10.88 % |

Sous 3000 K l'écart au TACOT reste sous ~3 % ; il se creuse au genou de
sublimation, où la courbe B'c(B'g) devient raide. Moins de soufflage
pyrolytique ⇒ le char du MX-4926 survit à plus haute température.

> **Les points à 0.01 atm au-delà de ~3300 K sont absents du tableau à
> dessein.** Le solveur y renvoie 200 ou 1000, qui ne sont pas des valeurs de
> B'c : ce sont les `char_amount = max(100·B'g, 200)` injectés par
> `Thermodynamics::surfaceMassBalance` pour simuler une surface ablative
> infinie. À basse pression et haute température le gaz pariétal devient du
> carbone pur, le dénominateur (y_w,C − y_c,C) s'annule et B'c n'est plus borné
> que par cette quantité finie. Physiquement : sublimation franche, le transfert
> de masse ne limite plus rien. Ce plafond est commun aux trois matériaux.

### Sensibilité au domaine de composition

Le rendement en char du composite ne dépend **que** de la fraction de résine :

```
    Y_comp = w_fibre + w_charge + Y·w_résine  =  1 − (1 − Y)·w_résine
```

par la fermeture Σw = 1. Le partage fibres/charge disparaît — les deux sont du
carbone à perte de masse nulle. L'incertitude sur le renfort (§1 de
`composition_mx4926.md` : la plage annoncée 41–56 % se referme sur 47–56 %) est
donc **sans effet** sur la réponse matériau.

Reste la plage 31–37 % de résine, soit k ∈ [0.1621, 0.1998]. Effet sur le point
de fonctionnement :

| P [atm] | T [K] | k = 0.1621 | k = 0.1817 | k = 0.1998 | étendue |
|---|---|---|---|---|---|
| 0.01 | 2000 | 0.16738 | 0.16645 | 0.16560 | 1.07 % |
| 0.01 | 3000 | 0.24781 | 0.24789 | 0.24797 | 0.06 % |
| 1 | 2000 | 0.16740 | 0.16647 | 0.16562 | 1.07 % |
| 1 | 3000 | 0.18171 | 0.18187 | 0.18203 | 0.18 % |
| 1 | 3400 | 0.22192 | 0.22321 | 0.22418 | 1.02 % |
| 1 | 3700 | 0.61376 | 0.62848 | 0.64270 | 4.61 % |

**Toute la plage de composition annoncée vaut ~1 % sur B'c** hors genou de
sublimation. Ce n'est pas le paramètre à resserrer en priorité.

---

## Contenu du répertoire

| Fichier | Rôle |
|---|---|
| `composition_mx4926.md` | traçabilité : fermeture du domaine, héritage du gaz, k, densités |
| `composition_mx4926_verification.py` | rejoue les 6 contrôles du document |
| `mx4926_bprime.py` | génère les tables et les figures |
| `mx4926_bprime_Bg*.csv` | table B' complète, une par B'g (25 isobares × 300–5000 K) |
| `mx4926_bprime_Bg*.png` | isobares B'c(T) et h_w(T), une figure par B'g |
| `mx4926_bprime_bg_comparison.png` | influence de B'g à 1 atm |
| `mx4926_bprime_Bc_vs_Bg.csv/.png` | B'c(B'g) à T fixée + points de fonctionnement |
| `mx4926_bprime_steady_state.csv` | point de fonctionnement stationnaire, avec l'enveloppe k_min/k_max |

## Contrôles

```
checkmix mx4926-air                     27 espèces, 4 éléments, C(gr) solide   OK
bprime   mx4926-air                     table complète 300–5000 K              OK
identité bit à bit avec sc1008-air      B'g = 0, 0.1, 0.5, 2.0                 OK
B'c(300 K, B'g = 0, 1 atm) = 0.087427   limite C + O2 -> CO2, attendu 0.0874   OK
```
