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

year = (units.year).to(units.second)


config = configparser.ConfigParser()
config.read("param.config")


input_files = config["system"]["save_file"]

major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")

#active_cl = config["clones"].getboolean("active")
#n_clones = config["clones"].getint("n_clones")



def extract_number(arquivo):
    match = re.search(r'_(\d+)\.h5$', arquivo)
    return int(match.group(1)) if match else float('inf')


files_list_not_ord = glob.glob(f"{input_files}_*.h5")

file_list = sorted(files_list_not_ord, key=extract_number)

all_H_K = []


ni_met = len(major_bodies) + len(minor_bodies)
first = True
for file_name in tqdm(file_list):
    with h5py.File(file_name, "r") as hf:
        if first:
            time = hf["time"][:]
            index = hf["index"][:]
            met_group = hf["met_group"][:]
            # all_met_code.append(hf["met_code"][:])
            first = False


        e = hf["e"][:]
        inc = hf["i"][:]


    H_K = np.zeros((len(index) - ni_met, e.shape[1]))
    for t_tmp in range(e.shape[1]):
        for bd in range(ni_met,len(index)):
            H_K[bd-ni_met][t_tmp] = np.cos(inc[bd][t_tmp])*np.sqrt(1. - e[bd][t_tmp]**2)

    all_H_K.append(H_K)





















