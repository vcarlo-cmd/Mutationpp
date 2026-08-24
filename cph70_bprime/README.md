# Table B' du CPh70

**CPh70** — composite carbone/phénolique dense
70 % fibres de carbone / 30 % résine phénolique (volume solide), porosité 0.01,
fibres et résine identiques à celles du TACOT.

| | |
|--|--|
| ρ vierge | 1465.2 kg/m³ |
| ρ char | 1287.0 kg/m³ |
| gaz de pyrolyse libéré | 178.2 kg/m³ (12.2 % de la masse vierge) |
| couplage stationnaire | **k = B'g/B'c = 0.1385** |
| mélange | `data/mixtures/cph70-air_25.xml` (25 espèces + C(gr)) |
| gaz de pyrolyse | C:0.206, H:0.679, O:0.115 |
| char | C:1.0 |

```bash
cd cph70_bprime && python cph70_bprime.py
```

---


> Bloc de **pyrolyse** (cinétique), dans ses trois écritures — native,
> `F_i`, `Δρ_i` : [`mise_en_donnees_pyrolyse.md`](mise_en_donnees_pyrolyse.md).
> Attention au facteur 9 sur la phase B.

## « La table B'c(B'g) est-elle différente de celle du TACOT ? »

Il faut distinguer deux objets.

**1. La table comme fonction tabulée B'c(T, P, B'g) — identique.**
Pour un même triplet (T, P, B'g), CPh70 et TACOT donnent exactement le même B'c
(écart max mesuré : `0.000e+00`). Le bilan de masse de surface ne voit que les
compositions élémentaires de l'air, du gaz de pyrolyse et du char — toutes trois
inchangées.

**2. Le B'c réellement atteint par le matériau — différent.**
Et c'est là que l'intuition est juste : puisque B'c dépend de B'g, et que les
deux matériaux ne fonctionnent pas au même B'g, ils ne lisent pas la table au
même endroit.

En ablation stationnaire le matériau impose

$$B'_g = k\,B'_c, \qquad k = \frac{\rho_v - \rho_c}{\rho_c}$$

soit **k = 0.1385 pour le CPh70** contre **0.2727 pour le TACOT**. Le point de
fonctionnement est l'intersection de la courbe B'c(B'g) — commune — avec cette
droite, propre à chaque matériau.

C'est exactement ce que montre `cph70_bprime_Bc_vs_Bg.png` : les courbes colorées
sont communes, les ★ (CPh70) et ■ (TACOT) sont les intersections.

### Points de fonctionnement

| P [atm] | T [K] | CPh70 B'c | CPh70 B'g | TACOT B'c | TACOT B'g | écart B'c |
|---------|-------|----------|----------|-----------|-----------|-----------|
| 0.01 | 2000 | 0.17102 | 0.0237 | 0.16682 | 0.0455 | +2.5 % |
| 0.01 | 3000 | 0.25315 | 0.0351 | 0.25920 | 0.0707 | −2.3 % |
| 1 | 2000 | 0.17105 | 0.0237 | 0.16684 | 0.0455 | +2.5 % |
| 1 | 3000 | 0.18564 | 0.0257 | 0.18961 | 0.0517 | −2.1 % |
| 1 | 3400 | 0.22660 | 0.0314 | 0.24089 | 0.0657 | −5.9 % |
| 1 | 3700 | 0.62677 | 0.0868 | 0.78927 | 0.2152 | **−20.6 %** |

Sous 3000 K l'écart reste sous ~2.5 % ; il se creuse au genou de sublimation,
où la courbe B'c(B'g) devient raide.

> **Récession.** À B'c égal, la vitesse de recul n'est pas la même :
> `ṡ = B'c·ṁe/ρc`. Le CPh70 étant 5.85 fois plus dense en char, il recule
> d'autant moins vite. La colonne `recession_over_mdote_m3_per_kg` de
> `cph70_bprime_steady_state.csv` donne directement `B'c/ρc`.

---

## Fichiers produits

| Fichier | Contenu |
|---------|---------|
| `cph70_bprime_Bg{0p0,0p1,0p2,0p5,1p0,2p0}.csv` | table B' complète : 25 isobares (10⁻³ → 10³ atm) × 189 T (300–5000 K, pas 25 K) × 26 fractions molaires |
| `cph70_bprime_Bg*.png` | isobares B'c et h_w pour chaque B'g |
| `cph70_bprime_bg_comparison.png` | influence de B'g à 1 atm |
| `cph70_bprime_Bc_vs_Bg.csv` / `.png` | **B'c en fonction de B'g** à T fixée, avec les points de fonctionnement |
| `cph70_bprime_steady_state.csv` | point de fonctionnement stationnaire : `Bc_ss`, `Bg_ss`, `hw_ss`, plus `Bc_Bg0` (référence sans pyrolyse) et `B'c/ρc` |

### Limite de validité du point fixe

Au-delà de la sublimation complète (≈ 3900 K à 1 atm, ≈ 3300 K à 0.01 atm), le
point fixe **diverge** : B'g sature sur le maximum de la grille (10) et B'c sur
le plafond interne de `surfaceMassBalance`. Les lignes correspondantes de
`cph70_bprime_steady_state.csv` (B'g_ss = 10) sont des **artefacts numériques**,
pas une solution physique — il n'existe alors plus d'état stationnaire au sens
de ce bilan. Les filtrer avec `Bg_ss < 10`.

---

Voir `tacot_bprime/autre_materiau.md` pour le raisonnement complet et
`tacot_bprime/material_response.py` pour le bilan de phase solide générique.
