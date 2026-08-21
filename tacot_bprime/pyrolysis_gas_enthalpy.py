#!/usr/bin/env python
"""
Génère et trace les propriétés du gaz de pyrolyse du TACOT à l'équilibre
chimique : enthalpie h_g, masse molaire M, Cp, gamma, masse volumique rho
et viscosité mu, en fonction de la température, pour plusieurs pressions.
Utilise l'outil MutationPP `mppequil` comme moteur de calcul.

h_g est la grandeur qui manque à la table B' pour fermer le bilan
d'énergie de surface (avec h_w) et le terme source de pyrolyse en
profondeur (h_g - h_s).

Hypothèse : équilibre chimique du gaz, sans phase condensée (le gaz de
pyrolyse ne contient pas de carbone solide à l'équilibre dans la gamme de
température considérée).

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique) et une plage de températures de 200 à 4000 K.

Usage :
    python pyrolysis_gas_enthalpy.py

Prérequis :
    - Le binaire `mppequil` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/tacot-pyrogas.xml doit exister
    - matplotlib et numpy installés

Pour comparer ces résultats à la table de référence du classeur TACOT 3.0,
voir pyrolysis_gas_enthalpy_validation.py.
"""

import subprocess
import shutil
import sys
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paramètres du calcul
# ---------------------------------------------------------------------------
MPPEQUIL_CMD = "mppequil"
T_RANGE      = "200:25:4000"   # T de 200 à 4000 K, pas de 25 K
MIXTURE      = "tacot-pyrogas"
ELEM_COMP    = "tacot_pyro"

ONEATM = 101325.0   # Pa

# Pressions : 0.001 -> 1000 atm, espacement logarithmique (25 isobares)
PRESSURES_ATM = np.logspace(-3, 3, 25)

# Grandeurs demandées à mppequil (voir mppequil --help pour la liste complète)
#   0 : Th [K]          5 : Mw [kg/mol]      9 : Cp_eq [J/kg-K]
#   10: H  [J/kg]        18: gam_eq [-]       3 : rho [kg/m3]
#   32: mu [Pa-s]
MPP_INDICES = "0,5,9,10,18,3,32"

# colonnes de sortie, dans l'ordre de MPP_INDICES : (nom, facteur, unité)
QUANTITIES = [
    ("M",     1e3,  "kg/kmol"),
    ("Cp",    1e-3, "kJ/kg-K"),
    ("h",     1e-3, "kJ/kg"),
    ("gamma", 1.0,  "-"),
    ("rho",   1.0,  "kg/m3"),
    ("mu",    1e4,  "millipoise"),
]


# ---------------------------------------------------------------------------
# Fonctions utilitaires
# ---------------------------------------------------------------------------

def find_mppequil():
    """Localise le binaire mppequil (PATH ou répertoires build courants)."""
    cmd = shutil.which(MPPEQUIL_CMD)
    if cmd:
        return cmd
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "../build/src/apps/mppequil"),
        os.path.join(script_dir, "../../build/src/apps/mppequil"),
        "build/src/apps/mppequil",
        "../build/src/apps/mppequil",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def make_env():
    """Environnement d'exécution avec MPP_DATA_DIRECTORY défini."""
    env = os.environ.copy()
    if "MPP_DATA_DIRECTORY" not in env:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        env["MPP_DATA_DIRECTORY"] = os.path.abspath(
            os.path.join(script_dir, "../data"))
    return env


def run_mppequil(mppequil_path, pressure_pa):
    """Exécute mppequil à une pression donnée (Pa) et retourne la sortie texte."""
    cmd = [
        mppequil_path,
        "-T", T_RANGE,
        "-P", str(pressure_pa),
        "-m", MPP_INDICES,
        "--elem-comp", ELEM_COMP,
        MIXTURE,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, env=make_env())
    if result.returncode != 0:
        print(f"\nERREUR mppequil à P = {pressure_pa:.2f} Pa :")
        print(result.stderr)
        sys.exit(1)
    return result.stdout


def parse_output(output):
    """
    Parse la sortie de mppequil.
    Retourne un ndarray de colonnes [Th, M, Cp, h, gamma, rho, mu] (unités SI).
    """
    lines = [l.strip() for l in output.strip().splitlines() if l.strip()]
    data = []
    for line in lines[1:]:
        try:
            data.append([float(v) for v in line.split()])
        except ValueError:
            continue
    return np.array(data)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

# Pressions à tracer : uniquement les puissances de 10
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)   # 0.001, 0.01, 0.1, 1, 10, 100, 1000 atm


def plot_gas_properties(all_data, pressures_atm):
    """
    Trace les six propriétés du gaz de pyrolyse en fonction de T, une
    courbe par puissance de 10 de pression.
    """
    data_map = {P: d for P, d in zip(pressures_atm, all_data)}
    n = len(PLOT_PRESSURES_ATM)
    colors = plt.get_cmap("plasma", n + 1)

    fig, axes = plt.subplots(2, 3, figsize=(19, 10))
    fig.suptitle(
        "Propriétés du gaz de pyrolyse du TACOT à l'équilibre  "
        r"($P \in [10^{-3},\,10^3]$ atm)",
        fontsize=14,
    )
    ylabels = {"M": "M [kg/kmol]", "Cp": r"$C_p$ [kJ/kg/K]",
               "h": r"$h_g$ [kJ/kg]", "gamma": r"$\gamma$",
               "rho": r"$\rho$ [kg/m³]", "mu": r"$\mu$ [millipoise]"}

    for ax, (name, fac, unit) in zip(axes.ravel(), QUANTITIES):
        col = 1 + [q[0] for q in QUANTITIES].index(name)
        for idx, P_atm in enumerate(PLOT_PRESSURES_ATM):
            closest = min(pressures_atm,
                          key=lambda p: abs(np.log10(p) - np.log10(P_atm)))
            data = data_map[closest]
            Tw = data[:, 0]
            val = data[:, col] * fac
            exp = int(round(np.log10(P_atm)))
            lbl = rf"$10^{{{exp}}}$ atm" if exp != 0 else "1 atm"
            ax.plot(Tw, val, color=colors(idx), lw=1.8, label=lbl)

        ax.set_xlabel(r"$T$ [K]")
        ax.set_ylabel(ylabels[name])
        ax.grid(True, ls="--", alpha=0.35)
        if name in ("rho", "mu"):
            ax.set_yscale("log")
        ax.legend(fontsize=7, ncol=2, title="Pression", title_fontsize=7)

    plt.tight_layout()
    out_png = "pyrolysis_gas_enthalpy.png"
    plt.savefig(out_png, dpi=150)
    print(f"Figure sauvegardée : {out_png}")
    plt.close()


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # 1. Localiser le binaire
    mppequil_path = find_mppequil()
    if mppequil_path is None:
        print(
            f"Binaire '{MPPEQUIL_CMD}' introuvable.\n"
            "Compilez MutationPP :\n"
            "  cmake -B build -DCMAKE_BUILD_TYPE=Release .\n"
            "  cmake --build build --target mppequil\n"
            "puis ajoutez build/src/apps/ au PATH."
        )
        sys.exit(1)
    print(f"Binaire trouvé : {mppequil_path}\n")

    # 2. Calcul pour chaque pression
    all_data = []
    for P_atm in PRESSURES_ATM:
        P_pa = P_atm * ONEATM
        print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
              end=" ", flush=True)
        data = parse_output(run_mppequil(mppequil_path, P_pa))
        all_data.append(data)
        print(f"{len(data)} points")

    # 3. Sauvegarde CSV globale
    out_csv = "pyrolysis_gas_enthalpy.csv"
    header = ["T_K"] + [f"{name}[{unit}]" for name, _, unit in QUANTITIES]
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["P_atm"] + header)
        for P_atm, data in zip(PRESSURES_ATM, all_data):
            for row in data:
                T = row[0]
                vals = [row[1 + i] * fac for i, (_, fac, _) in enumerate(QUANTITIES)]
                writer.writerow(
                    [f"{P_atm:.6g}", f"{T:.6g}"] + [f"{v:.6e}" for v in vals]
                )
    print(f"\nTable complète sauvegardée : {out_csv}")

    # 4. Visualisation
    plot_gas_properties(all_data, PRESSURES_ATM)
