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

f_plot = pathlib.Path("plot").resolve()
f_orbit = f_plot / "orbit"
f_orbit.mkdir(parents=True, exist_ok=True)

config = configparser.ConfigParser()
config.read("param.config")


input_files = config["system"]["save_file"]

major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")

#active_cl = config["clones"].getboolean("active")
#n_clones = config["clones"].getint("n_clones")






files_list_not_ord = glob.glob(f"{input_files}_*.h5")

file_list = sorted(files_list_not_ord, key=extract_number)


data_sets = {}

state_keys = [
    "x",
    "y",
]


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
            index  = [indx.decode() for indx in index]
            met_group = [met.decode() for met in met_group]
            n_G1 = met_group.count("G1")
            n_G1A = met_group.count("G1A")
            first = False
            xtmp = hf["x"][:]
            ytmp = hf["y"][:]
            ecc = hf["e"][:]
            x = np.zeros((xtmp.shape[0],xtmp.shape[1]*len(file_list)))
            y = np.zeros((xtmp.shape[0],xtmp.shape[1]*len(file_list)))
        else:
            xtmp = hf["x"][:]
            ytmp = hf["y"][:]
            ecc = hf["e"][:]


    for t_tmp in range(xtmp.shape[1]):
        for bd in range(0,len(index)):
            if ecc[bd][t_tmp] < 1.:
                x[bd][t_tmp + xtmp.shape[1]*n_file] = xtmp[bd][t_tmp]
                y[bd][t_tmp + xtmp.shape[1]*n_file] = ytmp[bd][t_tmp]
            else:
                x[bd][t_tmp + xtmp.shape[1]*n_file] = None
                y[bd][t_tmp + xtmp.shape[1]*n_file] = None

    n_file = n_file + 1


a_values = np.arange(n_G1)
b_values = np.arange(n_G1A)

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

lim_time = len(x[0]) // 10


for bd in range(ni_met,len(index)):
    if met_group [bd-ni_met] == "G1":
        color = cmap_G1(norm_G1(a_values[n_G1]))
        tic_names_G1.append(index[bd])
        n_G1 = n_G1 + 1
    else:
        color = cmap_G1A(norm_G1A(b_values[n_G1A]))
        tic_names_G1A.append(index[bd])
        n_G1A = n_G1A + 1

    ax.plot(x[bd][:lim_time]/au.value, y[bd][:lim_time]/au.value, color=color)
    ax.plot(
        x[index.index("Jupiter")][:lim_time]/au.value,
        y[index.index("Jupiter")][:lim_time]/au.value,
        color="black"
    )
    ax.plot(
        x[index.index("Sun")][:lim_time]/au.value,
        y[index.index("Sun")][:lim_time]/au.value,
        color="black"
    )



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

plt.axis("equal")

ax.set_xlabel("x (au)")
ax.set_ylabel("y (au)")

#ax.legend()
ax.grid(True)

plt.tight_layout()

figure_name = f_orbit / "all_orbit.png"
plt.savefig(figure_name)


ax.clear()

n_G1 = 0
n_G1A = 0





for bd in range(ni_met,len(index)):
    fig, ax = plt.subplots(figsize=(7, 6))

    fig.suptitle(f"{index[bd]} - {met_group [bd-ni_met]}", fontsize=16)

    if met_group [bd-ni_met] == "G1":
        color = cmap_G1(norm_G1(a_values[n_G1]))
        n_G1 = n_G1 + 1
    else:
        color = cmap_G1A(norm_G1A(b_values[n_G1A]))
        n_G1A = n_G1A + 1

    ax.plot(x[bd]/au.value, y[bd]/au.value, color=color)
    ax.plot(
        x[index.index("Jupiter")]/au.value,
        y[index.index("Jupiter")]/au.value,
        color="black"
    )
    ax.plot(
        x[index.index("Sun")]/au.value,
        y[index.index("Sun")]/au.value,
        color="black"
    )

    plt.axis("equal")
    ax.set_xlabel("x (au)")
    ax.set_ylabel("y (au)")
    ax.grid(True)

    plt.tight_layout()
    figure_name = f_orbit / f"orb_{index[bd]}.png"
    plt.savefig(figure_name)
    ax.clear()









