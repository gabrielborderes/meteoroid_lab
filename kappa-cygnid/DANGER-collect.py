import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob


save_file = "out/ephemerides"
output_file = "out/ephemerides_combined.h5"

file_list = sorted(glob.glob(f"{save_file}_*.h5"))  # Ordena para garantir sequência


all_index, all_time = [], []
all_x, all_y, all_z = [], [], []
all_vx, all_vy, all_vz = [], [], []
all_a, all_e, all_I = [], [], []
all_Omega, all_omega, all_f, all_M = [], [], [], []
all_met_code, all_met_group = [], []


for file_name in file_list:
    with h5py.File(file_name, "r") as hf:
        all_index.append(hf["index"][:])
        all_time.append(hf["time"][:])
        all_x.append(hf["x"][:])
        all_y.append(hf["y"][:])
        all_z.append(hf["z"][:])
        all_vx.append(hf["vx"][:])
        all_vy.append(hf["vy"][:])
        all_vz.append(hf["vz"][:])
        all_a.append(hf["a"][:])
        all_e.append(hf["e"][:])
        all_I.append(hf["i"][:])
        all_Omega.append(hf["Omega"][:])
        all_omega.append(hf["omega"][:])
        all_f.append(hf["f"][:])
        all_M.append(hf["M"][:])
        all_met_code.append(hf["met_code"][:])
        all_met_group.append(hf["met_group"][:])

print(f"{len(file_list)} arquivos carregados e concatenados.")

# Convertendo listas para arrays contínuos
index = np.concatenate(all_index)
time = np.concatenate(all_time)
x = np.concatenate(all_x)
y = np.concatenate(all_y)
z = np.concatenate(all_z)
vx = np.concatenate(all_vx)
vy = np.concatenate(all_vy)
vz = np.concatenate(all_vz)
a = np.concatenate(all_a)
e = np.concatenate(all_e)
I = np.concatenate(all_I)
Omega = np.concatenate(all_Omega)
omega = np.concatenate(all_omega)
f = np.concatenate(all_f)
M = np.concatenate(all_M)
met_code = np.concatenate(all_met_code)
met_group = np.concatenate(all_met_group)


with h5py.File(output_file, "w") as hf:
    hf.create_dataset("index", data=index)
    hf.create_dataset("time", data=time)
    hf.create_dataset("x", data=x)
    hf.create_dataset("y", data=y)
    hf.create_dataset("z", data=z)
    hf.create_dataset("vx", data=vx)
    hf.create_dataset("vy", data=vy)
    hf.create_dataset("vz", data=vz)
    hf.create_dataset("a", data=a)
    hf.create_dataset("e", data=e)
    hf.create_dataset("i", data=I)
    hf.create_dataset("Omega", data=Omega)
    hf.create_dataset("omega", data=omega)
    hf.create_dataset("f", data=f)
    hf.create_dataset("M", data=M)
    hf.create_dataset("met_code", data=met_code)
    hf.create_dataset("met_group", data=met_group)

print(f"Dados combinados salvos em {output_file}")

plt.figure(figsize=(8, 6))
plt.scatter(x, y, s=1, alpha=0.5, color="blue")  # s=1 define tamanho dos pontos, alpha=0.5 deixa mais transparente
plt.xlabel("X ")
plt.ylabel("Y ")
plt.title("Órbitas das Partículas (x vs y)")
plt.grid()
plt.savefig("out/orbit_plot.png", dpi=300)


plt.figure(figsize=(8, 6))
plt.scatter(time, x, s=1, alpha=0.5, color="blue")  # s=1 define tamanho dos pontos, alpha=0.5 deixa mais transparente
plt.xlabel("T ")
plt.ylabel("X ")
plt.title("Órbitas das Partículas (x vs y)")
plt.grid()
plt.savefig("out/time_plot.png", dpi=300)


print("Figura salva em out/orbit_plot.png")
