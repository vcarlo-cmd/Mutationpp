#!/usr/bin/env python
"""
Génère et trace la table B'c de l'AVCOAT (5026-39HC/G) dans l'air.
Utilise le binaire C++ MutationPP `bprime` comme moteur de calcul.

Matériau : résine époxy-novolac phénolique + microsphères creuses SiO₂
Système élémentaire : C-H-O-N-Si

Composition du gaz de pyrolyse (Goldstein 1992, base phénolique comme PICA) :
    CH4:0.13, CO:0.23, CO2:0.11, H2:0.12, H2O:0.39, C2H4:0.02 (mol. frac.)
    → éléments (moles) : C:0.51, H:1.62, O:0.84

Composition du char (Milos & Chen 2014, fractions massiques) :
    C:0.568, O:0.147, Si:0.285
    → proportions molaires : C:4.73, O:0.92, Si:1.01

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique, 25 isobares) et une plage de températures de 300 à 5000 K,
pour quatre valeurs de B'g : 0, 0.5, 1.0, 2.0.

Usage :
    python avcoat_bprime.py

Prérequis :
    - Le binaire `bprime` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/avcoat-air.xml doit exister
    - matplotlib et numpy installés
"""

import subprocess
import shutil
import sys
import os
import csv
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paramètres du calcul
# ---------------------------------------------------------------------------
BPRIME_CMD = "bprime"
T_RANGE    = "300:25:5000"   # Tw de 300 à 5000 K, pas 25 K
MIXTURE    = "avcoat-air"
BL_COMP    = "air"
PYRO_COMP  = "avcoat_pyro"
CHAR_COMP  = "avcoat_char"
CHAR_ELEM  = "C"

# Valeurs de B'g à calculer
BG_VALUES = [0.0, 0.5, 1.0, 2.0]

ONEATM = 101325.0  # Pa

# 25 pressions log-espacées : 0.001 → 1000 atm
PRESSURES_ATM = np.logspace(-3, 3, 25)

# 7 pressions puissances de 10 pour les graphiques isobares
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)


# ---------------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------------

def find_bprime():
    """Localise le binaire bprime (PATH ou répertoires build courants)."""
    cmd = shutil.which(BPRIME_CMD)
    if cmd:
        return cmd
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "../build/src/apps/bprime"),
        os.path.join(script_dir, "../../build/src/apps/bprime"),
        "build/src/apps/bprime",
        "../build/src/apps/bprime",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def _bprime_cmd(bprime_path, T_range, pressure_pa, bg):
    """Construit la commande bprime pour une plage de températures donnée."""
    return [
        bprime_path,
        "-T", T_range,
        "-P", str(pressure_pa),
        "-b", str(bg),
        "-m", MIXTURE,
        "-bl", BL_COMP,
        "-py", PYRO_COMP,
        "-char", CHAR_COMP,
        "-char-elem", CHAR_ELEM,
    ]


def _run_cmd(cmd, timeout_s=30):
    """Exécute une commande et retourne (stdout, ok). Retourne ('', False) si timeout."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        if result.returncode != 0:
            return "", False
        return result.stdout, True
    except subprocess.TimeoutExpired:
        return "", False


def run_bprime(bprime_path, pressure_pa, bg):
    """
    Exécute bprime pour une pression (Pa) et un B'g donnés.

    Stratégie robuste :
    1. Essai rapide sur la plage complète (timeout 45 s).
    2. Si accrochage : repasse point par point — les points accrochés sont
       interpolés linéairement entre voisins valides.
    """
    # ── Essai rapide ──────────────────────────────────────────────────────
    stdout, ok = _run_cmd(_bprime_cmd(bprime_path, T_RANGE, pressure_pa, bg),
                          timeout_s=45)
    if ok:
        return stdout

    # ── Reprise point par point ────────────────────────────────────────────
    # Reconstruit la grille de températures à partir de T_RANGE "T1:dT:T2"
    parts = T_RANGE.split(":")
    T1, dT, T2 = float(parts[0]), float(parts[1]), float(parts[2])
    temps = []
    T = T1
    while T <= T2 + 1e-6:
        temps.append(T)
        T += dT

    header = None
    rows = {}   # T → liste de valeurs

    for T in temps:
        stdout_pt, ok_pt = _run_cmd(
            _bprime_cmd(bprime_path, str(int(T)), pressure_pa, bg),
            timeout_s=5
        )
        if ok_pt:
            lines = [l.strip() for l in stdout_pt.strip().splitlines() if l.strip()]
            if header is None:
                header = lines[0]
            if len(lines) >= 2:
                rows[T] = lines[1]   # ligne de données

    if not rows:
        print(f"\nAvertissement : aucun point convergé à P={pressure_pa:.0f} Pa,"
              f" B'g={bg}. Pression ignorée.")
        return ""

    # Interpolation linéaire des points manquants
    valid_temps = sorted(rows.keys())
    all_rows = {}
    for T in temps:
        if T in rows:
            all_rows[T] = rows[T]
        else:
            # Cherche les voisins valides les plus proches
            lo = max((t for t in valid_temps if t < T), default=None)
            hi = min((t for t in valid_temps if t > T), default=None)
            if lo is None or hi is None:
                continue   # hors plage, ignorer
            # Interpolation numérique colonne par colonne
            vals_lo = [float(v) for v in rows[lo].split()]
            vals_hi = [float(v) for v in rows[hi].split()]
            alpha = (T - lo) / (hi - lo)
            interp = [(1 - alpha) * a + alpha * b
                      for a, b in zip(vals_lo, vals_hi)]
            interp[0] = T   # forcer Tw exact
            all_rows[T] = " ".join(f"{v:.6e}" for v in interp)

    # Reconstruction du flux de sortie au format attendu par parse_output
    lines = [header] + [all_rows[T] for T in sorted(all_rows.keys())]
    n_interp = len(temps) - len(rows)
    if n_interp > 0:
        print(f"[{n_interp} pts interpolés]", end=" ")
    return "\n".join(lines)


def parse_output(output):
    """
    Parse la sortie de bprime.
    Retourne (header: list[str], data: np.ndarray).
    """
    lines = [l.strip() for l in output.strip().splitlines() if l.strip()]
    header = [h.strip('"') for h in lines[0].split()]
    data = []
    for line in lines[1:]:
        try:
            data.append([float(v) for v in line.split()])
        except ValueError:
            continue
    return header, np.array(data)


def bg_label(bg):
    """Étiquette propre pour B'g."""
    return rf"$B'_g = {bg}$"


def bg_filename(bg):
    """Suffixe sûr pour les noms de fichiers."""
    return str(bg).replace(".", "p")


# ---------------------------------------------------------------------------
# Visualisation — table isobares (un graphique par B'g)
# ---------------------------------------------------------------------------

def plot_isobar_table(all_data, pressures_atm, bg, out_png):
    """
    Trace la table B'c classique (isobares) pour un B'g donné.
    Gauche : B'c vs Tw (log), Droite : h_w vs Tw.
    """
    data_map = {P: d for P, d in zip(pressures_atm, all_data)}
    n = len(PLOT_PRESSURES_ATM)
    colors = plt.get_cmap("plasma", n + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        rf"Table B' — AVCOAT dans l'air  ({bg_label(bg)},  "
        r"$P \in [10^{-3},\,10^3]$ atm)",
        fontsize=13,
    )

    for idx, P_atm in enumerate(PLOT_PRESSURES_ATM):
        closest = min(pressures_atm, key=lambda p: abs(np.log10(p) - np.log10(P_atm)))
        header, data = data_map[closest]
        Tw = data[:, 0]
        Bc = data[:, 1]
        hw = data[:, 2]

        exp = int(round(np.log10(P_atm)))
        lbl = rf"$10^{{{exp}}}$ atm" if exp != 0 else "1 atm"
        color = colors(idx)
        ax1.plot(Tw, Bc, color=color, lw=2, label=lbl)
        ax2.plot(Tw, hw, color=color, lw=2, label=lbl)

    ax1.set_yscale("log")
    ax1.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(rf"Taux d'ablation char $B'_c$  ({bg_label(bg)})")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    ax2.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(rf"Enthalpie de paroi $h_w$  ({bg_label(bg)})")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"  Figure isobares sauvegardée : {out_png}")
    plt.close()


# ---------------------------------------------------------------------------
# Visualisation — comparaison des B'g à 1 atm
# ---------------------------------------------------------------------------

def plot_bg_comparison(results_by_bg):
    """
    Trace B'c et h_w vs T_w à 1 atm pour toutes les valeurs de B'g.
    """
    colors = plt.get_cmap("tab10", len(BG_VALUES))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        r"Influence de $B'_g$ sur la table B' — AVCOAT dans l'air (P = 1 atm)",
        fontsize=13,
    )

    for idx, bg in enumerate(BG_VALUES):
        pressures_atm, all_data = results_by_bg[bg]
        _, data = min(
            zip(pressures_atm, all_data),
            key=lambda x: abs(np.log10(x[0])),
        )
        header, data = data
        Tw = data[:, 0]
        Bc = data[:, 1]
        hw = data[:, 2]

        color = colors(idx)
        ax1.plot(Tw, Bc, color=color, lw=2, label=bg_label(bg))
        ax2.plot(Tw, hw, color=color, lw=2, label=bg_label(bg))

    ax1.set_yscale("log")
    ax1.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(r"Taux d'ablation $B'_c$ vs $B'_g$")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=10, loc="upper left")

    ax2.set_xlabel(r"Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"Enthalpie de paroi $h_w$ vs $B'_g$")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=10, loc="upper left")

    plt.tight_layout()
    out_png = "avcoat_bprime_bg_comparison.png"
    plt.savefig(out_png, dpi=150)
    print(f"\nFigure de comparaison B'g sauvegardée : {out_png}")
    plt.close()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    bprime_path = find_bprime()
    if bprime_path is None:
        print(
            f"Binaire '{BPRIME_CMD}' introuvable.\n"
            "Compilez MutationPP :\n"
            "  cmake -B build -DCMAKE_BUILD_TYPE=Release .\n"
            "  cmake --build build --target bprime\n"
            "puis ajoutez build/src/apps/ au PATH."
        )
        sys.exit(1)
    print(f"Binaire trouvé : {bprime_path}\n")

    results_by_bg = {}

    for bg in BG_VALUES:
        print(f"=== B'g = {bg} ===")
        all_data = []

        for P_atm in PRESSURES_ATM:
            P_pa = P_atm * ONEATM
            print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
                  end=" ", flush=True)
            output = run_bprime(bprime_path, P_pa, bg)
            header, data = parse_output(output)
            all_data.append((header, data))
            print(f"{len(data)} points")

        # Sauvegarde CSV
        out_csv = f"avcoat_bprime_Bg{bg_filename(bg)}.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Bg", "P_atm"] + header)
            for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
                for row in data:
                    writer.writerow(
                        [f"{bg}", f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row]
                    )
        print(f"  Table CSV sauvegardée : {out_csv}")

        # Graphique isobares pour ce B'g
        out_png = f"avcoat_bprime_Bg{bg_filename(bg)}.png"
        plot_isobar_table(all_data, PRESSURES_ATM, bg, out_png)

        results_by_bg[bg] = (PRESSURES_ATM, all_data)
        print()

    # Graphique de comparaison entre B'g à 1 atm
    plot_bg_comparison(results_by_bg)

    print("\nTerminé.")
