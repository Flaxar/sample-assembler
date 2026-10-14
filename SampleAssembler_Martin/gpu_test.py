import torch
from ase.io import read, write
from ase.optimize import LBFGS
from mace.calculators import mace_mp

# 1. Bezpečná detekce dostupného hardwaru (podle doporučení Copilota)
if torch.cuda.is_available():
    device_name = "cuda"
    gpu_name = torch.cuda.get_device_name(0)
    print(f"[INFO] Paráda! PyTorch vidí GPU: {gpu_name}")
    print("[INFO] Relaxace MACE poběží na akcelerátoru.")
else:
    device_name = "cpu"
    print("[VAROVÁNÍ] PyTorch nevidí žádné kompatibilní GPU (ROCm/CUDA).")
    print("[VAROVÁNÍ] Ujistěte se, že máte správně nainstalovaný ROCm PyTorch.")
    print("[INFO] Záchranná brzda: Spouštím MACE na pomalém CPU...")

# 2. Načtení geometrie
# atoms = read("sestaveny_vzorek.xyz")

# 3. Inicializace univerzálního potenciálu s detekovaným zařízením
# Nyní se device rovná "cuda" jen tehdy, pokud je to bezpečné
# macenet = mace_mp(model="medium", default_dtype="float32", device=device_name)
# atoms.calc = macenet

# 4. Spuštění statické relaxace
print("[INFO] Začínám minimalizaci energie...")
# optimizer = LBFGS(atoms, trajectory='relaxace.traj')
# optimizer.run(fmax=0.05)

# 5. Uložení výsledku
# write("zrelaxovany_vzorek_MLIP.xyz", atoms)
print("[INFO] Hotovo, vzorek je uložen.")