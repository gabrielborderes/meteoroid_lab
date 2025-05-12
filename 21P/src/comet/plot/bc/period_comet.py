import pathlib
from tqdm import tqdm
import numpy as np
import h5py
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import configparser

from astropy.constants import GM_sun, au
import re
import glob
from astropy import units




def extract_number(arquivo):
    match = re.search(r'_(\d+)\.h5$', arquivo)
    return int(match.group(1)) if match else float('inf')

year = (units.year).to(units.second)



config = configparser.ConfigParser()
config.read("config/config_ress.ini")

out_folder = config["system"]["out_folder"]
save_file = config["system"]["save_file"]
save_time = config["system"]["save_time"]

out = pathlib.Path(out_folder).resolve()

input_files = out / save_file

major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")


files_list_not_ord = glob.glob(f"{input_files}_*.h5")
file_list = sorted(files_list_not_ord, key=extract_number)

ni_met = len(major_bodies) + len(minor_bodies)
time = None
met_group = None

file_name = out / save_time
with h5py.File(file_name, "r") as hf:
        time = hf["time"][:]
        index = hf["index"][:]
        met_group = hf["met_group"][:]

index  = [indx.decode() for indx in index]
met_group = [met.decode() for met in met_group]
n_G1 = met_group.count("G1")
n_G1A = met_group.count("G1A")

first = True
n_file = 0
for file_name in tqdm(file_list):
    with h5py.File(file_name, "r") as hf:
        if first:
            first = False
            ecc = hf["e"][:]
            xc = hf["x"][:]
            yc = hf["y"][:]
            zc = hf["z"][:]
            per_tmp = hf["Per"][:]
            per = np.zeros((ecc.shape[1]*len(file_list)))
            ec = np.zeros((ecc.shape[1]*len(file_list)))
            d = np.zeros((ecc.shape[1]*len(file_list)))
        else:
            ecc = hf["e"][:]
            xc = hf["x"][:]
            yc = hf["y"][:]
            zc = hf["z"][:]
            per_tmp = hf["Per"][:]

    n_c = index.index('1900 Y1')
    for t_tmp in range(ecc.shape[1]):
        per[t_tmp + ecc.shape[1]*n_file] = per_tmp[n_c][t_tmp]
        ec[t_tmp + ecc.shape[1]*n_file] = ecc[n_c][t_tmp]
        d[t_tmp + ecc.shape[1]*n_file] = np.sqrt(xc[n_c][t_tmp]**2+yc[n_c][t_tmp]**2+zc[n_c][t_tmp]**2)

    n_file = n_file + 1

ti = -145
tf = -135

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time/year, per/year, color="blue")
ax.axhline(y=7.11552, color='r', linestyle='--')
ax.set_xlim([ti,tf])
ax.set_ylabel("T (year)")
ax.set_xlabel("time (year)")
plt.tight_layout()
figure_name = f"per_fast1.png"
plt.savefig(figure_name)



fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time/year, d/au.value, color="blue")
ax.set_xlim([ti,tf])
ax.set_ylabel("d(au)")
ax.set_xlabel("time (year)")
plt.tight_layout()
figure_name = f"d_fast1.png"
plt.savefig(figure_name)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time/year, ec, color="blue")

ax.set_ylabel("e")
ax.set_xlabel("time (year)")
ax.set_xlim([ti,tf])
plt.tight_layout()
figure_name = f"ec_fast1.png"
plt.savefig(figure_name)


# Criando o gráfico com um tamanho específico
fig, ax1 = plt.subplots(figsize=(10, 6))

# Primeiro gráfico - Curva per
ax1.plot(time/year, per/year, color="blue", label="per/year")
ax1.axhline(y=7.11552, color='r', linestyle='--')
ax1.set_xlim([ti, tf])
ax1.set_xlabel("time (year)")
ax1.set_ylabel("T (year)", color='blue')

# Criando o segundo eixo y (lado direito)
ax2 = ax1.twinx()
ax2.plot(time/year, d/au.value, color="green", label="d (au)")
ax2.set_ylabel("d (au)", color='green')

# Ajustando layout e salvando a figura
plt.tight_layout()
figure_name = "combined_plot1.png"
plt.savefig(figure_name)







