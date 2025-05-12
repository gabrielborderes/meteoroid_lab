import seaborn as sns
import pathlib
from tqdm import tqdm
import numpy as np
import h5py
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import configparser
from matplotlib.colors import LogNorm

from astropy.constants import GM_sun, au
import re
import glob
from astropy import units as units

m_2_au = (units.m).to(units.AU)


def extract_number(arquivo):
    match = re.search(r'_(\d+)\.h5$', arquivo)
    return int(match.group(1)) if match else float('inf')

year = (units.year).to(units.second)



config = configparser.ConfigParser()
config.read("config/config.ini")

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



file_name = out / save_time
with h5py.File(file_name, "r") as hf:
        t = hf["time"][:]
        index = hf["index"][:]


index  = [indx.decode() for indx in index]


t = t[:int(len(t)/5)]

inc = np.zeros((len(index),len(t)))

t_s = 0
t_e = 0
for file_name in tqdm(file_list):
    with h5py.File(file_name, "r") as hf:
        i_tmp = hf["i"][:]

    t_e = int(t_e +  i_tmp.shape[1])
    if t_e > len(t):
        inc[1:,t_s:len(t)] = i_tmp[1:,:len(t)-t_s]
        break

    inc[1:,t_s:t_e] = i_tmp[1:,:]
    t_s = t_e

bins = np.linspace(0, 50, 30)


heatmap_data = np.zeros((len(bins) - 1, len(t)))


inc = np.degrees(inc)

i_len = np.zeros((len(bins) - 1))
for i in range(len(bins) - 1):
    i_len[i] = bins[i + 1]/2. + bins[i]/2.


for j in tqdm(range(inc.shape[1])):

    for i in range(ni_met,inc.shape[0]):
        i_eval = inc[i, j]

        for k in range(len(bins) - 1):
            if bins[k] <= i_eval < bins[k + 1]:
                heatmap_data[k, j] += 1

#heatmap_data[heatmap_data==0.0] = np.nan

# fig, ax1 = plt.subplots(figsize=(10, 6))
# sns.heatmap(heatmap_data, cmap='jet_r', cbar_kws={'label': 'Particles'})
#
# ax1.set_title('')
# ax1.set_xlabel('Time (year)')
# ax1.set_ylabel('a (au)')
#
#
# xticks = np.linspace(0, len(t) - 1, 20, dtype=int)
# yticks = np.linspace(0, len(a_len) - 1, 10, dtype=int)
# ytick_labels = np.round(a_len[yticks], 1)
# ax1.set_yticks(yticks, ytick_labels)
# ax1.invert_yaxis()
# ax1.invert_xaxis()
# ax1.set_xticks(xticks, np.round(t[xticks]/year))
#
#
#
#
# ax2 = ax1.twinx()
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_ylim(ax1.get_ylim())
# ax3 = ax2.twiny()
# ax3.set_ylim(ax2.get_ylim())
#
#
#
#
#
#
# ax3.plot(t/year,a[13,:], color='black')




fig, ax1 = plt.subplots(figsize=(10, 6))

extent = [t[0]/year, t[-1]/year, i_len[-1], i_len[0]]  # inverte y para bater com heatmap
# img = ax1.imshow(heatmap_data, aspect='auto', extent=extent, cmap='jet_r', origin='upper')
img = ax1.imshow(heatmap_data, aspect='auto', extent=extent,
                 cmap='jet_r', origin='upper',
                 norm=LogNorm(vmin=heatmap_data[heatmap_data > 0].min(), vmax=heatmap_data.max()))

ax1.set_xlabel('Time (year)')
ax1.set_ylabel('inc (degrees)')


#ax1.plot(t / year, inc[14, :], color='black')
ax1.scatter(t / year, inc[14, :], s=0.2,color='black')
ax1.axhline(y=35.5, color='r', linestyle='--')
ax1.axhline(y=41.7, color='r', linestyle='--')
plt.colorbar(img, ax=ax1, label='Particles')

ax1.invert_yaxis()
ax1.invert_xaxis()
ax1.set_ylim([0,50])

#ax2.set_yticks([])
plt.savefig(f"data/i_log-zoom.png", bbox_inches="tight")







