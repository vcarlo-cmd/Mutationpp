#!/usr/bin/env python
"""
Génère et trace la table B' du carbone (graphite) dans l'air.
Utilise le binaire C++ MutationPP `bprime` comme moteur de calcul.

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique) et une plage de températures de 300 à 5000 K.

Usage :
    python carbon_bprime.py

Prérequis :
    - Le binaire `bprime` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/carbon-air.xml doit exister
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
BPRIME_CMD = "bprime"       # nom du binaire (ou chemin absolu)
T_RANGE    = "300:100:5000" # Tw de 300 à 5000 K, pas de 100 K
BG         = "0.0"          # débit pyrolyse adimensionné (0 = graphite pur)
MIXTURE    = "carbon-air"
BL_COMP    = "air"
PYRO_COMP  = "pyro"

ONEATM = 101325.0   # Pa

# Pressions : 0.001 -> 1000 atm, espacement logarithmique (25 isobares)
PRESSURES_ATM = np.logspace(-3, 3, 25)


# ---------------------------------------------------------------------------
# Fonctions utilitaires
# ---------------------------------------------------------------------------

def find_bprime():
    """Localise le binaire bprime (PATH ou répertoires build courants)."""
    cmd = shutil.which(BPRIME_CMD)
    if cmd:
        return cmd
    # Recherche dans les répertoires build typiques, relatifs au script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "../../build/src/apps/bprime"),
        os.path.join(script_dir, "../../../build/src/apps/bprime"),
        "build/src/apps/bprime",
        "../build/src/apps/bprime",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def run_bprime(bprime_path, pressure_pa):
    """Exécute bprime à une pression donnée (Pa) et retourne la sortie texte."""
    cmd = [
        bprime_path,
        "-T", T_RANGE,
        "-P", str(pressure_pa),
        "-b", BG,
        "-m", MIXTURE,
        "-bl", BL_COMP,
        "-py", PYRO_COMP,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\nERREUR bprime à P = {pressure_pa:.2f} Pa :")
        print(result.stderr)
        sys.exit(1)
    return result.stdout


def parse_output(output):
    """
    Parse la sortie de bprime.
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

def plot_bprime_table(all_data, pressures_atm):
    """
    Trace la table B' classique :
      - Gauche  : B'c vs Tw, une courbe par pression
      - Droite  : h_w vs Tw, une courbe par pression
    """
    n = len(pressures_atm)
    cmap = cm.colormaps.get_cmap("plasma").resampled(n)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "Table B' — Graphite dans l'air  "
        r"($B'_g = 0$,  $P \in [10^{-3},\,10^3]$ atm)",
        fontsize=13
    )

    for idx, (P_atm, (header, data)) in enumerate(zip(pressures_atm, all_data)):
        Tw    = data[:, 0]
        Bc    = data[:, 1]
        hw    = data[:, 2]   # MJ/kg (déjà converti par bprime)
        color = cmap(idx)
        lbl   = f"{P_atm:.3g} atm"

        ax1.plot(Tw, Bc, color=color, lw=1.5, label=lbl)
        ax2.plot(Tw, hw, color=color, lw=1.5, label=lbl)

    for ax, ylabel, title in [
        (ax1, r"$B'_c$ (adimensionné)",   r"Taux d'ablation char $B'_c$"),
        (ax2, r"$h_w$ [MJ/kg]",           r"Enthalpie de paroi $h_w$"),
    ]:
        ax.set_xlabel("Température de paroi $T_w$ [K]")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(fontsize=6, ncol=2, loc="upper left",
                  title="Pression", title_fontsize=7)

    plt.tight_layout()
    out_png = "carbon_bprime_table.png"
    plt.savefig(out_png, dpi=150)
    print(f"Figure sauvegardée : {out_png}")
    plt.show()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # 1. Localiser le binaire
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

    # 2. Calcul pour chaque pression
    all_data = []
    for P_atm in PRESSURES_ATM:
        P_pa = P_atm * ONEATM
        print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
              end=" ", flush=True)
        output = run_bprime(bprime_path, P_pa)
        header, data = parse_output(output)
        all_data.append((header, data))
        print(f"{len(data)} points")

    # 3. Sauvegarde CSV globale
    out_csv = "carbon_bprime_table.csv"
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
