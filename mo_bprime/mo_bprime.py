#!/usr/bin/env python
"""
Génère et trace la table B' du molybdène (Mo, applicable au TZM) dans l'air.
Utilise le binaire C++ MutationPP `bprime` (bilan de masse généralisé,
options -char/-char-elem) comme moteur de calcul.

La table est calculée pour une plage de pressions de 0.001 à 1000 atm
(espacement logarithmique) et une plage de températures de 300 à 5000 K.

Deux tables sont produites :

  1. Équilibre complet (mo_bprime_*) : mélange data/mixtures/mo-air.xml,
     toutes les phases condensées admises, comme pour les autres matériaux.
     Avec un excès de métal, l'équilibre fixe l'oxygène en MoO2 solide :
     B'c = 0 tant que MoO2 est stable (~2000 K à 1 atm). C'est une BORNE
     BASSE : en réalité la face externe de l'oxyde voit l'air, forme MoO3
     qui s'évapore dès ~1000 K, et le TZM s'oxyde vite.

  2. Surface nue (mo_bprime_nu_*) : même mélange sans les oxydes condensés
     MoO2(cr), MoO3(cr), MoO3(L). Tout l'oxygène qui atteint la paroi part
     en oxydes gazeux : régime limité par la diffusion, B'c = y_O M_Mo /
     (3 M_O) = 0.466 (Mo + 3/2 O2 -> MoO3). C'est une BORNE HAUTE, à
     retenir au-dessus de ~1000-1100 K (MoO3 volatil) ; en dessous,
     l'oxydation est lente et limitée par la réaction.

Au-delà de la température où MoO2 devient instable, les deux tables
coïncident. Le bilan est suivi sur N (masse nette gazéifiée, cf.
sic_bprime/README.md) ; pour un char monoélément comme Mo, le résultat
est identique au suivi sur Mo tant que seul le métal est condensé.

Usage :
    python mo_bprime.py

Prérequis :
    - Le binaire `bprime` doit être dans le PATH ou dans build/src/apps/
    - Le fichier data/mixtures/mo-air.xml doit exister
    - MPP_DATA_DIRECTORY pointe sur data/ (sinon ../data est utilisé)
    - matplotlib et numpy installés
"""

import subprocess
import shutil
import sys
import os
import csv
import re
import tempfile
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paramètres du calcul
# ---------------------------------------------------------------------------
BPRIME_CMD = "bprime"           # nom du binaire (ou chemin absolu)
T_RANGE    = "300:25:5000"      # Tw de 300 à 5000 K, pas de 25 K
BG         = "0.0"              # débit pyrolyse adimensionné (0 = Mo pur)
MIXTURE    = "mo-air"
MIXTURE_NU = "mo-air-nu"        # variante sans oxydes condensés (écrite à la volée)
OXIDES_CONDENSES = ["MoO2(cr)", "MoO3(cr)", "MoO3(L)"]
BL_COMP    = "air"
CHAR_COMP  = "mo"               # composition élémentaire du char (Mo:1)
CHAR_ELEM  = "N"                # élément de suivi pour B'c (masse nette gazéifiée)

ONEATM = 101325.0   # Pa

# Pressions : 0.001 -> 1000 atm, espacement logarithmique (25 isobares)
PRESSURES_ATM = np.logspace(-3, 3, 25)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.environ.get("MPP_DATA_DIRECTORY",
                          os.path.join(SCRIPT_DIR, "..", "data"))


# ---------------------------------------------------------------------------
# Fonctions utilitaires
# ---------------------------------------------------------------------------

def find_bprime():
    """Localise le binaire bprime (PATH ou répertoires build courants)."""
    cmd = shutil.which(BPRIME_CMD)
    if cmd:
        return cmd
    candidates = [
        os.path.join(SCRIPT_DIR, "../build/src/apps/bprime"),
        os.path.join(SCRIPT_DIR, "../../build/src/apps/bprime"),
        "build/src/apps/bprime",
        "../build/src/apps/bprime",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def write_bare_surface_mixture(tmpdir):
    """Écrit mo-air-nu.xml : mo-air.xml privé de ses oxydes condensés."""
    txt = open(os.path.join(DATA_DIR, "mixtures", MIXTURE + ".xml")).read()
    txt = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
    for sp in OXIDES_CONDENSES:
        txt = re.sub(r"(?<=\s)" + re.escape(sp) + r"(?=\s)", "", txt)
    with open(os.path.join(tmpdir, MIXTURE_NU + ".xml"), "w") as f:
        f.write(txt)


def run_bprime_mo(bprime_path, pressure_pa, mixture, cwd=None):
    """Exécute bprime (char Mo) à une pression donnée (Pa) et retourne la sortie texte."""
    cmd = [
        bprime_path,
        "-T", T_RANGE,
        "-P", str(pressure_pa),
        "-b", BG,
        "-m", mixture,
        "-bl", BL_COMP,
        "-char", CHAR_COMP,
        "-char-elem", CHAR_ELEM,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:
        print(f"\nERREUR bprime ({mixture}) à P = {pressure_pa:.2f} Pa :")
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


def write_tables(prefix, all_data):
    """Table complète + tables B'c et h_w au format long (unités SI)."""
    header = all_data[0][0]
    out_csv = f"{prefix}_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["P_atm"] + header)
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            for row in data:
                writer.writerow(
                    [f"{P_atm:.6g}"] + [f"{v:.6e}" for v in row]
                )
    print(f"Table complète sauvegardée : {out_csv}")

    # Table B'c au format long, unites SI : Tw_K, P_bar, Bc
    # (nTw x nP ; B'g = 0 fixe, pas de colonne Bg)
    out_csv = f"{prefix}_bc_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Tw_K", "P_bar", "Bc"])
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            P_bar = P_atm * ONEATM / 1.0e5
            for row in data:
                writer.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}",
                                  f"{row[1]:.6e}"])
    print(f"Table B'c (nTw x nP) sauvegardée : {out_csv}")

    # Table h_w au format long, unites SI : Tw_K, P_bar, hw_Jkg (nTw x nP)
    out_csv = f"{prefix}_hw_table.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Tw_K", "P_bar", "hw_Jkg"])
        for P_atm, (_, data) in zip(PRESSURES_ATM, all_data):
            P_bar = P_atm * ONEATM / 1.0e5
            for row in data:
                writer.writerow([f"{row[0]:.6g}", f"{P_bar:.6g}",
                                  f"{row[2] * 1.0e6:.6e}"])
    print(f"Table h_w (nTw x nP) sauvegardée : {out_csv}")


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

# Pressions à tracer : uniquement les puissances de 10
PLOT_PRESSURES_ATM = np.logspace(-3, 3, 7)   # 0.001, 0.01, 0.1, 1, 10, 100, 1000 atm


def plot_bprime_table(all_eq, all_nu, pressures_atm):
    """
    Trace la table B' Mo :
      - Gauche  : B'c vs Tw en échelle log10 ; trait plein = équilibre
                  complet (borne basse), tirets = surface nue (borne haute)
      - Droite  : h_w vs Tw (équilibre complet)
    Seules les 7 pressions puissances de 10 sont tracées.
    """
    eq_map = {P: d for P, d in zip(pressures_atm, all_eq)}
    nu_map = {P: d for P, d in zip(pressures_atm, all_nu)}

    n = len(PLOT_PRESSURES_ATM)
    colors = plt.get_cmap("plasma", n + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        r"Table B' — Molybdène (Mo, TZM) dans l'air  "
        r"($B'_g = 0$,  $P \in [10^{-3},\,10^3]$ atm)",
        fontsize=13
    )

    for idx, P_atm in enumerate(PLOT_PRESSURES_ATM):
        closest = min(pressures_atm, key=lambda p: abs(np.log10(p) - np.log10(P_atm)))
        _, deq = eq_map[closest]
        _, dnu = nu_map[closest]

        exp = int(round(np.log10(P_atm)))
        lbl = rf"$10^{{{exp}}}$ atm" if exp != 0 else "1 atm"

        color = colors(idx)
        ax1.plot(deq[:, 0], deq[:, 1], color=color, lw=2, label=lbl)
        ax1.plot(dnu[:, 0], dnu[:, 1], color=color, lw=1.2, ls="--")
        ax2.plot(deq[:, 0], deq[:, 2], color=color, lw=2, label=lbl)

    ax1.set_yscale("log")
    ax1.set_ylim(1e-2, 1e3)
    ax1.set_xlabel("Température de paroi $T_w$ [K]")
    ax1.set_ylabel(r"$B'_c$ (échelle log$_{10}$)")
    ax1.set_title(r"$B'_c$ — plein : équilibre complet, tirets : surface nue")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    ax2.set_xlabel("Température de paroi $T_w$ [K]")
    ax2.set_ylabel(r"$h_w$ [MJ/kg]")
    ax2.set_title(r"Enthalpie de paroi $h_w$ (équilibre complet)")
    ax2.grid(True, ls="--", alpha=0.4)
    ax2.legend(fontsize=9, loc="upper left", title="Pression", title_fontsize=9)

    plt.tight_layout()
    out_png = "mo_bprime_table.png"
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

    # 2. Calcul pour chaque pression, équilibre complet puis surface nue
    all_eq, all_nu = [], []
    with tempfile.TemporaryDirectory() as tmpdir:
        write_bare_surface_mixture(tmpdir)
        for P_atm in PRESSURES_ATM:
            P_pa = P_atm * ONEATM
            print(f"  P = {P_atm:8.4g} atm  ({P_pa:12.2f} Pa) ...",
                  end=" ", flush=True)
            all_eq.append(parse_output(run_bprime_mo(bprime_path, P_pa, MIXTURE)))
            all_nu.append(parse_output(run_bprime_mo(bprime_path, P_pa, MIXTURE_NU,
                                                     cwd=tmpdir)))
            print(f"{len(all_eq[-1][1])} points")

    # 3. Sauvegarde des tables
    print("\nÉquilibre complet (borne basse) :")
    write_tables("mo_bprime", all_eq)
    print("\nSurface nue, limite de diffusion (borne haute) :")
    write_tables("mo_bprime_nu", all_nu)
    print()

    # 4. Visualisation
    plot_bprime_table(all_eq, all_nu, PRESSURES_ATM)
