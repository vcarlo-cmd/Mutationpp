#!/usr/bin/env python
"""
Génère et trace la table B' de la silice (SiO2) dans l'air.
Utilise le binaire C++ MutationPP `bprime_silica` comme moteur de calcul.

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique) et une plage de températures de 300 à 5000 K.

Usage :
    python silice_bprime.py

Prérequis :
    - Le binaire `bprime_silica` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/silica-air.xml doit exister
    - matplotlib et numpy installés
"""

import subprocess
import shutil
import sys
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Paramètres du calcul
# ---------------------------------------------------------------------------
BPRIME_CMD = "bprime_silica"    # nom du binaire (ou chemin absolu)
T_RANGE    = "300:25:5000"      # Tw de 300 à 5000 K, pas de 25 K
BG         = "0.0"              # débit pyrolyse adimensionné (0 = silice pure)
MIXTURE    = "silica-air"
BL_COMP    = "air"

ONEATM = 101325.0   # Pa

# Pressions : 0.001 -> 1000 atm, espacement logarithmique (25 isobares)
PRESSURES_ATM = np.logspace(-3, 3, 25)


# ---------------------------------------------------------------------------
# Fonctions utilitaires
# ---------------------------------------------------------------------------

def find_bprime_silica():
    """Localise le binaire bprime_silica (PATH ou répertoires build courants)."""
    cmd = shutil.which(BPRIME_CMD)
    if cmd:
        return cmd
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "../build/src/apps/bprime_silica"),
        os.path.join(script_dir, "../../build/src/apps/bprime_silica"),
        "build/src/apps/bprime_silica",
        "../build/src/apps/bprime_silica",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def run_bprime_silica(bprime_path, pressure_pa):
    """Exécute bprime_silica à une pression donnée (Pa) et retourne la sortie texte."""
    cmd = [
        bprime_path,
        "-T", T_RANGE,
        "-P", str(pressure_pa),
        "-b", BG,
        "-m", MIXTURE,
        "-bl", BL_COMP,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\nERREUR bprime_silica à P = {pressure_pa:.2f} Pa :")
        print(result.stderr)
        sys.exit(1)
    return result.stdout


def parse_output(output):
    """
    Parse la sortie de bprime_silica.
    La première ligne est l'en-tête avec les noms entre guillemets.
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


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

# Pressions à tracer : uniquement les puissances de 10
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)   # 0.001, 0.01, 0.1, 1, 10, 100, 1000 atm


def plot_bprime_table(all_data, pressures_atm):
    """
    Trace la table B' silice :
      - Gauche  : B'c vs Tw en échelle log10, une courbe par puissance de 10
      - Droite  : h_w vs Tw, une courbe par puissance de 10
    Seules les 7 pressions puissances de 10 sont tracées.
    """
    data_map = {P: d for P, d in zip(pressures_atm, all_data)}

    n = len(PLOT_PRESSURES_ATM)
    colors = plt.get_cmap("plasma", n + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        r"Table B' — Silice (SiO$_2$) dans l'air  "
        r"($B'_g = 0$,  $P \in [10^{-3},\,10^3]$ atm)",
        fontsize=13
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
    ax1.set_xlabel("Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(r"Taux d'ablation silice $B'_c$")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    ax2.set_xlabel("Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"Enthalpie de paroi $h_w$")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    plt.tight_layout()
    out_png = "silice_bprime_table.png"
    plt.savefig(out_png, dpi=150)
    print(f"Figure sauvegardée : {out_png}")
    plt.show()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # 1. Localiser le binaire
    bprime_path = find_bprime_silica()
    if bprime_path is None:
        print(
            f"Binaire '{BPRIME_CMD}' introuvable.\n"
            "Compilez MutationPP :\n"
            "  cmake -B build -DCMAKE_BUILD_TYPE=Release .\n"
            "  cmake --build build --target bprime_silica\n"
            "puis ajoutez build/src/apps/ au PATH."
        )
        sys.exit(1)
    print(f"Binaire trouvé : {bprime_path}\n")

    # 2. Calcul pour chaque pression
    all_data = []
    for P_atm in PRESSURES_ATM:
        P_pa = P_atm * ONEATM
        print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
              end=" ", flush=True)
        output = run_bprime_silica(bprime_path, P_pa)
        header, data = parse_output(output)
        all_data.append((header, data))
        print(f"{len(data)} points")

    # 3. Sauvegarde CSV globale
    out_csv = "silice_bprime_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["P_atm"] + header)
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            for row in data:
                writer.writerow(
                    [f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row]
                )
    print(f"\nTable complète sauvegardée : {out_csv}")

    # 4. Visualisation
    plot_bprime_table(all_data, PRESSURES_ATM)
