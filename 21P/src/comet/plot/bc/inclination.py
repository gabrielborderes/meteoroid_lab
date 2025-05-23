import pathlib

from tqdm import tqdm
import numpy as np
import pandas as pd
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

year = (units.year).to(units.second)


config = configparser.ConfigParser()
#config.read("param.config")
config.read("config/config.ini")


f_plot = pathlib.Path("plot").resolve()
f_orbele = f_plot / "orbital_elements/inc"
f_orbele.mkdir(parents=True, exist_ok=True)



input_files = config["system"]["save_file"]
input_time = config["system"]["save_time"]
folder = config["system"]["out_folder"]


major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")

active_cl = config["clones"].getboolean("active")
n_clones = config["clones"].getint("n_clones")



def extract_number(arquivo):
    match = re.search(r'_(\d+)\.h5$', arquivo)
    return int(match.group(1)) if match else float('inf')




files_list_not_ord = glob.glob(f"{folder}/{input_files}_*.h5")
file_list = sorted(files_list_not_ord, key=extract_number)



with h5py.File(f"{folder}/{input_time}", "r") as hf:
    time = hf["time"][:]
    index = hf["index"][:]
    met_group = hf["met_group"][:]


met_group = [met.decode() for met in met_group]
n_G1 = met_group.count("G1")
n_G1A = met_group.count("G1A")
index  = [indx.decode() for indx in index]


ni_met = len(major_bodies) + len(minor_bodies)

first = True
n_file = 0

for file_name in tqdm(file_list):
    with h5py.File(file_name, "r") as hf:
        if first:
            first = False
            itmp = hf["i"][:]
            ecc = hf["e"][:]
            inc = np.zeros((itmp.shape[0],itmp.shape[1]*len(file_list)))
        else:
            itmp = hf["i"][:]
            ecc = hf["e"][:]


    for t_tmp in range(itmp.shape[1]):
        for bd in range(1,len(index)):
            if ecc[bd][t_tmp] < 1.:
                inc[bd][t_tmp + itmp.shape[1]*n_file] = itmp[bd][t_tmp]
            else:
                inc[bd][t_tmp + itmp.shape[1]*n_file] = None

    n_file = n_file + 1

if active_cl:

    fig, ax = plt.subplots(figsize=(10, 6))

    index_met = index[ni_met:]
    for bd, item in enumerate(index_met):
        item = item.strip()

        if "_cl_" in item:
            ax.plot(time/year, np.degrees(inc[bd+ni_met]), color=color_cl, alpha=0.2)
        else:
            ng = int(bd/(n_clones+1))
            if met_group[ng] == "G1":
                color_cl = "blue"
            else:
                color_cl = "red"
            #fig.suptitle(f"{index_met[bd]} - {met_group[bd]}", fontsize=16)
            ax.plot(time/year, np.degrees(inc[bd+ni_met]), color=color_cl)

    ax.set_ylabel("i (degrees)")
    ax.set_xlabel("time (year)")

    plt.tight_layout()
    figure_name = f_orbele / f"e.png"
    plt.savefig(figure_name)
    ax.clear()


else:

    # a_values = np.arange(n_G1)
    # b_values = np.arange(n_G1A)
    #
    # norm_G1 = mcolors.Normalize(vmin=a_values.min(), vmax=a_values.max())
    # norm_G1A = mcolors.Normalize(vmin=b_values.min(), vmax=b_values.max())
    #
    # cmap_G1 = cm.Blues
    # cmap_G1A = cm.Reds
    #
    # tic_names_G1 = []
    # tic_names_G1A = []
    #
    # fig, ax = plt.subplots(figsize=(10, 6))
    #
    # n_G1 = 0
    # n_G1A = 0
    # index
    #
    #
    #
    # for bd in range(ni_met,len(index)):
    #     if met_group [bd-ni_met] == "G1":
    #         color = cmap_G1(norm_G1(a_values[n_G1]))
    #         tic_names_G1.append(index[bd])
    #         n_G1 = n_G1 + 1
    #     else:
    #         color = cmap_G1A(norm_G1A(b_values[n_G1A]))
    #         tic_names_G1A.append(index[bd])
    #         n_G1A = n_G1A + 1
    #
    #     ax.plot(time/year, np.degrees(inc[bd]), color=color)
    #
    #
    #
    #
    # sm_G1 = cm.ScalarMappable(cmap=cmap_G1, norm=norm_G1)
    # sm_G1A = cm.ScalarMappable(cmap=cmap_G1A, norm=norm_G1A)
    # sm_G1.set_array([])
    # sm_G1A.set_array([])
    #
    #
    # cbar_G1 = plt.colorbar(sm_G1, ax=ax, orientation="vertical", pad=0.1)
    # cbar_G1A = plt.colorbar(sm_G1A, ax=ax, orientation="vertical", pad=0.02)
    #
    #
    # cbar_G1.set_ticks(a_values)
    # cbar_G1.set_ticklabels(tic_names_G1)
    #
    # cbar_G1A.set_ticks(b_values[:len(tic_names_G1A)])
    # cbar_G1A.set_ticklabels(tic_names_G1A)
    #
    # #plt.axis("equal")
    #
    # ax.set_ylabel("i (degrees)")
    # ax.set_xlabel("time (year)")
    #
    # #ax.legend()
    # ax.grid(True)
    #
    # plt.tight_layout()
    #
    # figure_name = f_orbele / "all_inc.png"
    # plt.savefig(figure_name)
    #
    # ax.clear()
    #
    # n_G1 = 0
    # n_G1A = 0




    #for bd in range(ni_met,len(index)):
    for bd in range(1,ni_met):


        print(f"{index[bd]}",time[np.argmax(inc[bd])],inc[bd][np.argmax(inc[bd])])
        fig, ax = plt.subplots(figsize=(7, 6))

        #fig.suptitle(f"{index[bd]} - {met_group [bd-ni_met]}", fontsize=16)
        fig.suptitle(f"{index[bd]}", fontsize=16)

        # if met_group [bd-ni_met] == "G1":
        #     ax.plot(time/year, np.degrees(inc[bd]), color="blue")
        # else:
        #     ax.plot(time/year, np.degrees(inc[bd]), color="red")

        ax.plot(time/year, np.degrees(inc[bd]), color="blue")

        #plt.axis("equal")
        ax.set_ylabel("i (degrees)")
        ax.set_xlabel("time (year)")
        #ax.grid(True)

        plt.tight_layout()
        figure_name = f_orbele / f"inc_{index[bd]}.png"
        plt.savefig(figure_name)
        ax.clear()









