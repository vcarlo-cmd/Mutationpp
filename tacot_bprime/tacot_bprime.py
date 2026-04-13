#!/usr/bin/env python
"""
Génère et trace la table B'c du TACOT (Theoretical Ablative Composite for Open
Testing) dans l'air.
Utilise le binaire C++ MutationPP `bprime` comme moteur de calcul.

Composition du gaz de pyrolyse (phénol, C6H5OH) :
    C:6, H:6, O:1 (fractions molaires élémentaires)

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique, 25 isobares) et une plage de températures de 300 à
5000 K, pour quatre valeurs de B'g : 0, 0.5, 1.0, 2.0.

Usage :
    python tacot_bprime.py

Prérequis :
    - Le binaire `bprime` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/tacot-air_35.xml doit exister
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
MIXTURE    = "tacot-air_35"
BL_COMP    = "air"
PYRO_COMP  = "tacot_pyro"
CHAR_COMP  = "tacot_char"
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


def make_env():
    """Retourne l'environnement d'exécution avec MPP_DATA_DIRECTORY défini."""
    env = os.environ.copy()
    if "MPP_DATA_DIRECTORY" not in env:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, "../data")
        env["MPP_DATA_DIRECTORY"] = os.path.abspath(data_dir)
    return env


def run_bprime(bprime_path, pressure_pa, bg):
    """Exécute bprime pour une pression (Pa) et un B'g donnés."""
    cmd = [
        bprime_path,
        "-T", T_RANGE,
        "-P", str(pressure_pa),
        "-b", str(bg),
        "-m", MIXTURE,
        "-bl", BL_COMP,
        "-py", PYRO_COMP,
        "-char", CHAR_COMP,
        "-char-elem", CHAR_ELEM,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if result.returncode != 0:
        print(f"\nERREUR bprime à P={pressure_pa:.2f} Pa, B'g={bg} :")
        print(result.stderr)
        sys.exit(1)
    return result.stdout


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
        rf"Table B' — TACOT dans l'air  ({bg_label(bg)},  "
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
        r"Influence de $B'_g$ sur la table B' — TACOT dans l'air (P = 1 atm)",
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
    out_png = "tacot_bprime_bg_comparison.png"
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
        out_csv = f"tacot_bprime_Bg{bg_filename(bg)}.csv"
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
        out_png = f"tacot_bprime_Bg{bg_filename(bg)}.png"
        plot_isobar_table(all_data, PRESSURES_ATM, bg, out_png)

        results_by_bg[bg] = (PRESSURES_ATM, all_data)
        print()

    # Graphique de comparaison entre B'g à 1 atm
    plot_bg_comparison(results_by_bg)

    print("\nTerminé.")
