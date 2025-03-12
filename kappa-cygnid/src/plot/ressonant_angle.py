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
f_ress = f_plot / "ressonance"
f_ress.mkdir(parents=True, exist_ok=True)

config = configparser.ConfigParser()
config.read("param_ress.config")

input_files = config["system"]["save_file"]

major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")

active_cl = config["clones"].getboolean("active")
n_clones = config["clones"].getint("n_clones")

files_list_not_ord = glob.glob(f"{input_files}_*.h5")
file_list = sorted(files_list_not_ord, key=extract_number)

ni_met = len(major_bodies) + len(minor_bodies)
time = None
met_group = None


p = int(3)
q = int(2)

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
            Mtmp = hf["M"][:]
            Otmp = hf["Omega"][:]
            otmp = hf["omega"][:]
            ecc = hf["e"][:]
            phi = np.zeros((Mtmp.shape[0],Mtmp.shape[1]*len(file_list)))

        else:
            Mtmp = hf["M"][:]
            Otmp = hf["Omega"][:]
            otmp = hf["omega"][:]
            ecc = hf["e"][:]

    n_j = index.index('Jupiter')
    for t_tmp in range(Mtmp.shape[1]):
        l_J = Otmp[n_j][t_tmp] + otmp[n_j][t_tmp] + Mtmp[n_j][t_tmp]

        for bd in range(1,len(index)):
            l_bd = Otmp[bd][t_tmp] + otmp[bd][t_tmp] + Mtmp[bd][t_tmp]
            varpi = Otmp[bd][t_tmp] + otmp[bd][t_tmp]
            if ecc[bd][t_tmp] < 1.:
                phi0 = ( p + q ) * l_J - p * l_bd - q * varpi
                phi[bd][t_tmp + Mtmp.shape[1]*n_file] = np.arctan2(np.sin(phi0),np.cos(phi0))
            else:
                phi[bd][t_tmp + Mtmp.shape[1]*n_file] = None

    n_file = n_file + 1




if active_cl:

    fig, ax = plt.subplots(figsize=(10, 6))

    index_met = index[ni_met:]
    for bd, item in enumerate(index_met):
        item = item.strip()
        if "_cl_" in item:
            ax.plot(time/year, np.degrees(phi[bd+ni_met]), color=color_cl, alpha=0.2)
        else:
            ng = int(bd/(n_clones+1))
            if met_group[ng] == "G1":
                color_cl = "blue"
            else:
                color_cl = "red"

            #fig.suptitle(f"{index_met[bd]} - {met_group[ng]}", fontsize=16)
            ax.plot(time/year, np.degrees(phi[bd+ni_met]), color=color_cl)


    plt.xlabel("Time (years)")
    plt.ylabel(r"$\varphi$ (degrees)")

    ax.set_xlim([-20,0])
    ax.set_ylim([-180,180])
    ax.grid(True)
    #ax.set_yscale('log')
    plt.tight_layout()
    figure_name = f_ress / f"ress_20.png"
    plt.savefig(figure_name)
    ax.clear()

else:
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



    for bd in range(ni_met,len(index)):
        if met_group [bd-ni_met] == "G1":
            color = cmap_G1(norm_G1(a_values[n_G1]))
            tic_names_G1.append(index[bd])
            n_G1 = n_G1 + 1
        else:
            color = cmap_G1A(norm_G1A(b_values[n_G1A]))
            tic_names_G1A.append(index[bd])
            n_G1A = n_G1A + 1

        ax.plot(time/year, np.degrees(phi[bd]), color=color)


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

    plt.xlabel("Time (years)")
    plt.ylabel(r"$\varphi$ (degrees)")

    #ax.legend()
    ax.grid(True)
    #ax.set_yscale('log')

    plt.tight_layout()
    figure_name = f_ress / "all_ress.png"
    plt.savefig(figure_name)

    ax.clear()

    n_G1 = 0
    n_G1A = 0


    for bd in range(ni_met,len(index)):
        fig, ax = plt.subplots(figsize=(7, 6))

        fig.suptitle(f"{index[bd]} - {met_group [bd-ni_met]}", fontsize=16)

        if met_group [bd-ni_met] == "G1":
            ax.plot(time/year, np.degrees(phi[bd]), color="blue")
        else:
            ax.plot(time/year, np.degrees(phi[bd]), color="red")



        #plt.axis("equal")
        plt.xlabel("Time (years)")
        plt.ylabel(r"$\varphi$ (degrees)")
        ax.grid(True)
        #ax.set_yscale('log')

        plt.tight_layout()
        figure_name = f_ress/ f"ress_{index[bd]}.png"
        plt.savefig(figure_name)
        ax.clear()
























