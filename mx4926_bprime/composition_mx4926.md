# Le MX-4926 : composition, fermeture du domaine, réponse matériau

Document de traçabilité, pendant de `../sc1008_bprime/resine_sc1008.md`, mais
un cran plus haut : la résine y était l'objet, ici c'est le **composite**.

Toutes les vérifications numériques de ce document sont rejouables :

```bash
python composition_mx4926_verification.py
```

---

## 0. Réponse courte

| Question | Réponse |
|---|---|
| **Matériau** | **MX-4926** — carbone/phénolique chargé, tuyère de moteur à propergol solide |
| **Renfort** | fibres de carbone **ex-rayon** (ex-cellulose), tissu **satin de 5 ou de 8** |
| **Matrice** | **résol** phénol-formaldéhyde **SC-1008** — la même que le PICA |
| **Charge** | **noir de carbone** |
| **Composition** | 41–56 / 31–37 / 11–16 % masse *(plages annoncées, non indépendantes — §1)* |
| **Gaz de pyrolyse** | **hérité du SC-1008** : `C:0.2526, H:0.6407, O:0.1068` |
| **Char** | **`C:1.0`** — les trois constituants carbonisent vers le carbone pur |
| **Table B'** | **rigoureusement celle du SC-1008**, identité vérifiée bit à bit |
| **Ce qui est propre au MX-4926** | le couplage **k = B'g/B'c = 0.1817** |

Ce document ne recalcule donc **pas** une composition de gaz : il justifie
qu'il n'y a rien à recalculer, et établit ce qui distingue réellement le
MX-4926 — sa **réponse matériau**.

---

## 1. La fermeture du domaine de composition

Les trois plages annoncées ne sont **pas indépendantes** :

```
    somme des bornes basses :  41 + 31 + 11  =   83 %
    somme des bornes hautes :  56 + 37 + 16  =  109 %
```

Toute composition réelle doit sommer à 100 %. Le domaine admissible est donc
l'intersection de la boîte des trois plages avec l'hyperplan Σw = 1 — un
**polygone à 5 sommets** dans le plan (w_fibre, w_résine) :

| | fibre | résine | charge |
|---|---|---|---|
| S1 | 53.00 | 31.00 | 16.00 |
| S2 | 56.00 | 31.00 | 13.00 |
| S3 | 56.00 | 33.00 | 11.00 |
| S4 | 52.00 | 37.00 | 11.00 |
| S5 | 47.00 | 37.00 | 16.00 |

**Plages effectives** après fermeture :

| Constituant | annoncé | effectif | |
|---|---|---|---|
| renfort | 41 – 56 % | **47 – 56 %** | **borne basse déplacée** |
| matrice | 31 – 37 % | 31 – 37 % | intacte |
| charge | 11 – 16 % | 11 – 16 % | intacte |

> Les **41 %** de renfort annoncés ne sont atteignables avec **aucun** couple
> (résine, charge) admissible : il faudrait 59 % de résine + charge, alors que le
> maximum admissible est 37 + 16 = 53 %. C'est le seul défaut de fermeture du jeu de données.

**Centroïde exact** du polygone (centroïde d'aire, pas une moyenne de grille) :

```
    fibre 52.2024 %      résine 34.1667 %      charge 13.6310 %
```

C'est la composition de référence retenue partout ailleurs, arrondie à
**52.2 / 34.2 / 13.6**.

---

## 2. Pourquoi la table B' est celle du SC-1008

Le solveur de bilan de masse de surface (`Thermodynamics::surfaceMassBalance`)
ne consomme que **trois compositions élémentaires** — bord de couche limite,
gaz de pyrolyse, char — plus T, P et B'g. Cf.
`../tacot_bprime/mise_en_donnees_xml.md`.

### 2.1 Gaz de pyrolyse — hérité

Il est produit par la **résine seule** :

- le **renfort ex-rayon** est déjà carbonisé (une fibre de carbone est le
  produit d'une pyrolyse) : il ne dégaze pas ;
- le **noir de carbone** est du carbone quasi pur, thermiquement inerte
  jusqu'au domaine de sublimation : il ne dégaze pas non plus.

La matrice étant le SC-1008, la composition élémentaire du gaz est celle
établie dans `../sc1008_bprime/resine_sc1008.md` §4.4 :

```
    C:0.2526, H:0.6407, O:0.1068          (H/O = 6.00)
```

Le taux de charge et le rapport fibres/résine changent la **quantité** de gaz
produite — donc B'g, qui est un paramètre d'**entrée** de la table, balayé de
0 à 10 — pas sa composition.

### 2.2 Char — `C:1.0`

La pondération générale du char sur ses constituants

```
    n_E(char) = Σ_i (m_char,i / M_i) · ν_E,i
```

se réduit à `C:1.0` dès que **tous** les constituants carbonisent vers le même
élément. C'est le cas ici pour les trois : fibres ex-rayon, noir de carbone et
résol carbonisé sont du carbone. Les proportions n'entrent donc pas dans Y_c.

> C'est exactement ce qui rendrait le calcul faux avec des **fibres de silice**
> ou un char retenant H/O : le rapport entrerait alors dans Y_c et la table
> changerait. Cf. `../tacot_bprime/cph70_vs_tacot.md`.

### 2.3 Vérification

Même liste d'espèces, même gaz, même char ⇒ la table doit être identique
**caractère pour caractère**. Elle l'est, à quatre valeurs de B'g :

```
  [OK] B'g = 0.0  : tables identiques (149 lignes)
  [OK] B'g = 0.1  : tables identiques (149 lignes)
  [OK] B'g = 0.5  : tables identiques (149 lignes)
  [OK] B'g = 2.0  : tables identiques (149 lignes)
```

Même situation que `cph70-air_25.xml` face à `tacot-air_25.xml`. Le fichier
`data/mixtures/mx4926-air.xml` n'est justifié que par la **traçabilité** :
nommer le matériau.

Contrôle physique associé, indépendant : à 300 K et B'g = 0 dans l'air, un char
carboné doit donner la limite d'oxydation C + O₂ → CO₂.

```
    B'c(300 K, B'g = 0, 1 atm) = 0.087427          attendu 0.0874
```

Une valeur ≈ 200 signalerait l'absence de `C(gr)` dans la liste d'espèces.

---

## 3. Ce qui EST propre au MX-4926 : le rendement en char

### 3.1 Il ne dépend que de la fraction de résine

Le renfort et la charge ne perdent pas de masse ; seule la résine pyrolyse,
avec Y = 0.55 (`../sc1008_bprime/resine_sc1008.md` §5.2). Donc

```
    Y_comp = w_fibre + w_charge + Y · w_résine
```

et, en utilisant la fermeture w_fibre + w_charge = 1 − w_résine :

```
    Y_comp = 1 − (1 − Y) · w_résine
```

> **Le partage fibres/charge disparaît.** Deux constituants distincts, mais
> tous deux du carbone à perte de masse nulle : seule leur **somme** compte, et
> la fermeture la fixe. Le rendement en char du composite — donc k, donc le
> point de fonctionnement — ne dépend que de **w_résine**.

C'est ce qui rend l'incertitude sur le renfort (§1) sans conséquence : la seule
plage qui pilote la réponse matériau est 31–37 % de résine, et elle survit
intacte à la fermeture.

### 3.2 Le couplage k

En ablation stationnaire (matériau entièrement carbonisé reculant à vitesse
constante), les flux de gaz de pyrolyse et de char sont dans le rapport

```
    k = B'g / B'c = (ρ_v − ρ_c) / ρ_c = (1 − Y_comp) / Y_comp
```

**k est un rapport de masses** : ni les densités, ni la porosité n'y entrent.

| w_résine | Y_comp | k |
|---|---|---|
| 0.31 (borne basse) | 0.8605 | **0.1621** |
| **0.3417 (centroïde)** | **0.8463** | **0.1817** |
| 0.37 (borne haute) | 0.8335 | **0.1998** |

Situation dans le dépôt :

| Matériau | k | |
|---|---|---|
| TACOT (10/10/80 vol.) | 0.2727 | poreux, beaucoup de résine par unité de masse |
| PICA (SC-1008) | 0.2346 | `resine_sc1008.md` §7 : 0.19/0.81 |
| **MX-4926 (52/34/14 wt)** | **0.1817** | dense et chargé |
| CPh70 (69/30/1 vol.) | 0.1385 | dense et très riche en fibres |

Le MX-4926 se place entre le PICA — même résine, mais sans charge et à taux de
résine plus élevé — et le CPh70. Le noir de carbone agit exactement comme des
fibres du point de vue du bilan : il **abaisse** k en ajoutant de la masse qui
ne dégaze pas.

### 3.3 Densités

Elles n'entrent ni dans le XML, ni dans k — seulement dans la récession
`ṡ = B'c·ṁe/ρ_c`. Avec ρ = 1600 (fibres ex-cellulose, valeur du classeur
TACOT 3.0), 1250 (résol cuit) et 1800 kg/m³ (noir de carbone) :

| porosité | ρ_v | ρ_c | gaz libéré |
|---|---|---|---|
| 0.00 | 1480.8 | 1253.1 | 227.7 kg/m³ |
| **0.02** | **1451.2** | **1228.0** | **223.1 kg/m³** |
| 0.04 | 1421.5 | 1203.0 | 218.6 kg/m³ |

La porosité **0.02** n'est pas ajustée librement : c'est celle qui restitue la
densité publiée du MX-4926 (≈ 1.45 g/cm³). Le calcul inverse donne 0.0208.

Le gaz libéré vaut **15.4 % de la masse vierge**, contre 21.4 % pour le TACOT
et 12.2 % pour le CPh70 — cohérent avec la position de k.

---

## 4. Ce qui entre dans le XML

```xml
<element_compositions default="air">
    <composition name="air">N:0.79, O:0.21</composition>
    <composition name="mx4926_pyro">C:0.2526, H:0.6407, O:0.1068</composition>
    <composition name="mx4926_char">C:1.0</composition>
</element_compositions>
```

```bash
bprime -T 300:25:4000 -P 101325 -b 0.5 \
       -m mx4926-air -bl air -py mx4926_pyro -char mx4926_char -char-elem C
```

Liste d'espèces : celle de `sc1008-air.xml`, reprise à l'identique pour que les
deux tables soient comparables terme à terme (25 espèces de la table B' de
référence TACOT 3.0 + phénol + `C(gr)`).

---

## 5. Récapitulatif des sources

| Donnée | Valeur | Origine |
|---|---|---|
| Plages de composition | 41–56 / 31–37 / 11–16 % masse | donnée d'entrée |
| Nature du renfort | carbone ex-rayon, satin de 5 ou 8 | donnée d'entrée |
| Nature de la matrice | résol SC-1008 | donnée d'entrée |
| Nature de la charge | noir de carbone | donnée d'entrée |
| Gaz de pyrolyse | C:0.2526 H:0.6407 O:0.1068 | `resine_sc1008.md` §4.4 (fermeture, Y = 0.55) |
| Rendement en char de la résine | 0.55 | `resine_sc1008.md` §5.2 |
| Char | C:1.0 | §2.2 de ce document |
| Domaine admissible, centroïde | 52.20 / 34.17 / 13.63 % | §1, dérivé |
| **Y_comp, k** | **0.8463, 0.1817** | §3, dérivé |
| ρ fibres / résine / charge | 1600 / 1250 / 1800 kg/m³ | classeur TACOT 3.0 `Description!A11` ; valeurs usuelles |
| Densité publiée du MX-4926 | ≈ 1450 kg/m³ | fiche matériau |

---

## 6. Fichiers liés

| Fichier | Rôle |
|---|---|
| `composition_mx4926_verification.py` | rejoue les 6 vérifications de ce document |
| `mx4926_bprime.py` | tables B' et figures, points de fonctionnement |
| `README.md` | résultats et mode d'emploi |
| `../sc1008_bprime/resine_sc1008.md` | la résine — d'où sortent les trois nombres du gaz |
| `../tacot_bprime/mise_en_donnees_xml.md` | ce qu'il faut renseigner dans un XML de mélange |
| `../tacot_bprime/cph70_vs_tacot.md` | le même raisonnement d'héritage, pour le couple TACOT/CPh70 |
| `../tacot_bprime/material_response.py` | couplage k et lignes de fonctionnement |
| `../../data/mixtures/mx4926-air.xml` | le mélange |
| `../../data/mixtures/sc1008-pyrogas.xml` | gaz de pyrolyse seul, pour h_g(T, P) — valable tel quel |
