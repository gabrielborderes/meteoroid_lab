from tqdm import tqdm
import numpy as np
import pandas as pd
import h5py
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import configparser

import pathlib

from astropy.constants import GM_sun, au
import re
import glob
from astropy import units

year = (units.year).to(units.second)


config = configparser.ConfigParser()
config.read("param.config")



data_folder = pathlib.Path("plot/orbital_elem/").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()



input_files = "out/all_meteors/ephemerides"

major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")

#active_cl = config["clones"].getboolean("active")
#n_clones = config["clones"].getint("n_clones")



def extract_number(arquivo):
    match = re.search(r'_(\d+)\.h5$', arquivo)
    return int(match.group(1)) if match else float('inf')


files_list_not_ord = glob.glob(f"{input_files}_*.h5")

file_list = sorted(files_list_not_ord, key=extract_number)





ni_met = len(major_bodies) + len(minor_bodies)
time = None
met_group = None




first = True
n_file = 0
for file_name in tqdm(file_list):
    with h5py.File(file_name, "r") as hf:
        if first:
            time = hf["time"][:]
            index = hf["index"][:]
            met_group = hf["met_group"][:]
            # all_met_code.append(hf["met_code"][:])
            first = False
            etmp = hf["e"][:]
            e = np.zeros((etmp.shape[0],etmp.shape[1]*len(file_list)))
        else:
            etmp = hf["e"][:]


    for t_tmp in range(etmp.shape[1]):
        for bd in range(1,len(index)):
            if etmp[bd][t_tmp] < 1.:
                e[bd][t_tmp + etmp.shape[1]*n_file] = etmp[bd][t_tmp]
            else:
                e[bd][t_tmp + etmp.shape[1]*n_file] = None

    n_file = n_file + 1




a_values = np.arange(42)
b_values = np.arange(6)

norm_G1 = mcolors.Normalize(vmin=a_values.min(), vmax=a_values.max())
norm_G1A = mcolors.Normalize(vmin=b_values.min(), vmax=b_values.max())

cmap_G1 = cm.Blues
cmap_G1A = cm.Reds

tic_names_G1 = []
tic_names_G1A = []

fig, ax = plt.subplots(figsize=(10, 6))

n_G1 = 0
n_G1A = 0
index


index  = [indx.decode() for indx in index]
met_group = [met.decode() for met in met_group]


for bd in range(ni_met,len(index)):
    if met_group [bd-ni_met] == "G1":
        color = cmap_G1(norm_G1(a_values[n_G1]))
        tic_names_G1.append(index[bd])
        n_G1 = n_G1 + 1
    else:
        color = cmap_G1A(norm_G1A(b_values[n_G1A]))
        tic_names_G1A.append(index[bd])
        n_G1A = n_G1A + 1

    ax.plot(time/year, e[bd], color=color)




sm_G1 = cm.ScalarMappable(cmap=cmap_G1, norm=norm_G1)
sm_G1A = cm.ScalarMappable(cmap=cmap_G1A, norm=norm_G1A)
sm_G1.set_array([])
sm_G1A.set_array([])


cbar_G1 = plt.colorbar(sm_G1, ax=ax, orientation="vertical", pad=0.1)
cbar_G1A = plt.colorbar(sm_G1A, ax=ax, orientation="vertical", pad=0.02)


cbar_G1.set_ticks(a_values)
cbar_G1.set_ticklabels(tic_names_G1)

cbar_G1A.set_ticks(b_values[:len(tic_names_G1A)])
cbar_G1A.set_ticklabels(tic_names_G1A)

#plt.axis("equal")

ax.set_xlabel("e")
ax.set_ylabel("time (year)")

#ax.legend()
ax.grid(True)

plt.tight_layout()
figure_name = data_folder / "all_e.png"
plt.savefig(figure_name)

ax.clear()

n_G1 = 0
n_G1A = 0





for bd in range(ni_met,len(index)):
    fig, ax = plt.subplots(figsize=(7, 6))

    fig.suptitle(f"{index[bd]} - {met_group [bd-ni_met]}", fontsize=16)

    if met_group [bd-ni_met] == "G1":
        ax.plot(time/year, e[bd], color="blue")
    else:
        ax.plot(time/year, e[bd], color="red")



    #plt.axis("equal")
    ax.set_xlabel("e")
    ax.set_ylabel("time (year)")
    #ax.grid(True)

    plt.tight_layout()
    figure_name = data_folder / f"e_{index[bd]}.png"
    plt.savefig(figure_name)
    ax.clear()









