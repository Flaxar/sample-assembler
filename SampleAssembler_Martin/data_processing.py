import matplotlib.pyplot as plt

# Skript hledá soubor log.lammps ve stejné složce
filename = "log.lammps"

timesteps = []
pressures = []
temperatures = []
energies = []

print("Načítám data z log.lammps...")

with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
    in_table = False
    idx_step = -1
    idx_press = -1
    idx_temp = -1
    idx_etotal = -1

    for line in f:
        parts = line.split()
        if not parts:
            continue

        # Detekce začátku termodynamické tabulky
        if "Step" in parts and "Temp" in parts:
            in_table = True
            idx_step = parts.index("Step")
            idx_temp = parts.index("Temp")

            if "Press" in parts:
                idx_press = parts.index("Press")
            elif "Pxx" in parts:
                idx_press = parts.index("Pxx")

            if "etotal" in parts:
                idx_etotal = parts.index("etotal")
            elif "TotEng" in parts:
                idx_etotal = parts.index("TotEng")
            continue

        # Detekce konce tabulky
        if in_table and line.startswith("Loop time"):
            in_table = False

        # Zpracování samotných čísel
        if in_table:
            try:
                step = int(parts[idx_step])
                temp = float(parts[idx_temp])

                # Ošetření pro případ, že by nějaký sloupec v logu chyběl
                press = float(parts[idx_press]) if idx_press != -1 else 0.0
                etotal = float(parts[idx_etotal]) if idx_etotal != -1 else 0.0

                timesteps.append(step)
                temperatures.append(temp)
                pressures.append(press)
                energies.append(etotal)
            except ValueError:
                pass  # Přeskočí řádky, které nejsou čísla (např. varování v logu)

if not timesteps:
    print("Nepodařilo se najít žádná data. Zkontroluj log.lammps.")
else:
    print(f"Úspěšně načteno {len(timesteps)} kroků.")

    # Vytvoření 2x2 mřížky pro grafy
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Termalizace a relaxace vzorku', fontsize=16, fontweight='bold', y=0.98)

    # 1. Tlak vs Krok (Tvůj původní graf vlevo nahoře)
    axs[0, 0].plot(timesteps, pressures, color='#d62728', linewidth=1.5)
    axs[0, 0].set_title('Celkový tlak vs. Krok simulace')
    axs[0, 0].set_xlabel('Krok simulace (Timestep)')
    axs[0, 0].set_ylabel('Tlak [bar]')
    axs[0, 0].axhline(0, color='black', linewidth=1, alpha=0.8)
    axs[0, 0].grid(True, linestyle='--', alpha=0.7)

    # 2. Teplota vs Krok (Vpravo nahoře)
    axs[0, 1].plot(timesteps, temperatures, color='#ff7f0e', linewidth=1.5)
    axs[0, 1].set_title('Vývoj teploty')
    axs[0, 1].set_xlabel('Krok simulace (Timestep)')
    axs[0, 1].set_ylabel('Teplota [K]')
    axs[0, 1].grid(True, linestyle='--', alpha=0.7)

    # 3. Totální energie vs Krok (Vlevo dole)
    axs[1, 0].plot(timesteps, energies, color='#1f77b4', linewidth=1.5)
    axs[1, 0].set_title('Celková energie vs. Krok simulace')
    axs[1, 0].set_xlabel('Krok simulace (Timestep)')
    axs[1, 0].set_ylabel('Celková energie [eV]')
    axs[1, 0].grid(True, linestyle='--', alpha=0.7)

    # 4. Totální energie vs Teplota (Vpravo dole)
    # Zde používám malé tečky i čáru, hezky to ukazuje směr smyčky ohřevu a chlazení
    axs[1, 1].plot(temperatures, energies, color='#2ca02c', linewidth=1, marker='.', markersize=4)
    axs[1, 1].set_title('Celková energie vs. Teplota (Fázový prostor)')
    axs[1, 1].set_xlabel('Teplota [K]')
    axs[1, 1].set_ylabel('Celková energie [eV]')
    axs[1, 1].grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig("termalization_relaxation_plot.png", dpi=600)
    plt.show()