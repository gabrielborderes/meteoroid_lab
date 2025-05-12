# import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
import h5py
import seaborn as sns
import pathlib
import rebound
import numpy as np
import scipy.constants as constants
import matplotlib.font_manager as font_manager
from tqdm import tqdm
from astropy import units
import os
from astropy.constants import GM_sun, au
from matplotlib.colors import LogNorm
#try:
#    from mpi4py import MPI

#    comm = MPI.COMM_WORLD
#except ImportError:
#    class COMM_WORLD:
#        rank = 0
#        size = 1

#        def barrier(self):
#            pass


from config import (
    solar_system_objects,
    Nog,
    init_sim_file,
    min_size_log,
    max_size_log,
    factor_flush,
)


nb = len(solar_system_objects)
year = (units.year).to(units.second)

data_folder = pathlib.Path("data").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

particle_file = data_folder / "particles_21P.h5"



# sim = rebound.Simulation()
# sim.units = ("m", "s", "kg")
# date = "2000-01-01 00:00"
# for i in range(len(solar_system_objects)):
#     sim.add(solar_system_objects[i], date=date, hash=solar_system_objects[i])
#
# sim.move_to_com()
# # Configure simulation
# sim.integrator = "ias15"
# # Let's save it for next time
# # Note: sim.save_to_file() only saves the particle data, not the integrator settings, etc.
# # sim.collision = "direct"
# # sim.collision = "none"
# ps = sim.particles
#
#
# extra_sim = np.arange(0,59900.*24.*3600.,24.*3600.)
#
# x_bodies = np.zeros((nb, len(extra_sim)))
# y_bodies = np.zeros((nb, len(extra_sim)))
#
#
#
#
# for ti, t in tqdm(enumerate(extra_sim), total=len(extra_sim)):
#     sim.integrate(t)
#
#     for bind in range(nb):
#         p = ps[f"{solar_system_objects[bind]}"]
#         x_bodies[bind][ti] = p.x
#         y_bodies[bind][ti] = p.y





my_cmap = plt.get_cmap("viridis")


plot_folder = data_folder / f"map"
#if comm.rank == 0:
if not plot_folder.is_dir():
    plot_folder.mkdir()



with h5py.File(particle_file, 'r') as f:
    shell_part = f["shell_part"][:]
    N_activity = f["N_activity"][:]


N_part = shell_part.shape[0]
nb = len(solar_system_objects)
first = True
K_per = (2.*np.pi)/np.sqrt(GM_sun.value)

n_file = 0
for fl in tqdm(range(factor_flush)):
#for fl in tqdm(range(2)):

    print(f"Part_{fl} running")
    ephem_file = data_folder / f"part_{fl}/ephemerides_21P_all.h5"

    if first:
        first = False
        with h5py.File(ephem_file, "r") as hf:
            t = np.concatenate([
                hf["t"][()],
            ], axis=0)
            x = np.concatenate([
                hf["body_x"][()],
                hf["x"][()],
            ], axis=0)
            y = np.concatenate([
                hf["body_y"][()],
                hf["y"][()],
            ], axis=0)
            z = np.concatenate([
                hf["body_z"][()],
                hf["z"][()],
            ], axis=0)
            Per_tmp = np.concatenate([
                hf["body_Per"][()],
                hf["Per"][()],
            ], axis=0)
            # a = np.concatenate([
            #     hf["body_Per"][()],
            #     hf["a"][()],
            # ], axis=0)
            # e = np.concatenate([
            #     hf["body_Per"][()],
            #     hf["e"][()],
            # ], axis=0)

        dist_tmp = np.sqrt(x**2+y**2+z**2)

        Per = np.zeros((Per_tmp.shape[0],len(t)))

        dist = np.zeros((Per_tmp.shape[0],len(t)))

            # for i in range(a.shape[0]):
            #     for j in range(a.shape[1]):
            #         if i > 0:
            #             if e[i,j] < 1.0:
            #                 Per_n[i,j] = K_per*(a[i,j])**(3./2.)



        Per[:,:Per_tmp.shape[1]] = Per_tmp
        dist[:,:Per_tmp.shape[1]] = dist_tmp
        n_sam_fl = Per_tmp.shape[1]

    else:
        with h5py.File(ephem_file, "r") as hf:
            x = np.concatenate([
                hf["body_x"][()],
                hf["x"][()],
            ], axis=0)
            y = np.concatenate([
                hf["body_y"][()],
                hf["y"][()],
            ], axis=0)
            z = np.concatenate([
                hf["body_z"][()],
                hf["z"][()],
            ], axis=0)
            Per_tmp = np.concatenate([
                hf["body_Per"][()],
                hf["Per"][()],
            ], axis=0)
            # a = np.concatenate([
            #     hf["body_Per"][()],
            #     hf["a"][()],
            # ], axis=0)
            # e = np.concatenate([
            #     hf["body_Per"][()],
            #     hf["e"][()],
            # ], axis=0)

        dist_tmp = np.sqrt(x**2+y**2+z**2)

        Per[:,fl*n_sam_fl:(fl+1)*Per_tmp.shape[1]] = Per_tmp
            #Per_n[1:,fl*n_sam_fl:(fl+1)*a.shape[1]] = K_per*(a[1:,:])**(3./2.)

            # for i in range(a.shape[0]):
            #     for j in range(a.shape[1]):
            #         if i > 0:
            #             if e[i,j] < 1.0:
            #                 Per_n[i,fl*n_sam_fl + j] = K_per*(a[i,j])**(3./2.)
            #             else:
            #                 Per_n[i,fl*n_sam_fl + j] = 0.0

       #
        dist[:,fl*n_sam_fl:(fl+1)*Per_tmp.shape[1]] = dist_tmp


    x[x==0.0] = np.nan
    y[y==0.0] = np.nan
    #Per[Per==0.0] = np.nan
    dist[dist==0.0] = np.nan


    size = shell_part[:, 0]

    sizefont = "18"

    n_frames = x.shape[1]
    iter_inds = x.shape[1]
    #iter_inds = np.arange(comm.rank, n_frames, comm.size, dtype=np.int64)
    #N_process_part = len(iter_inds)

    cbar = None
    #pbar = tqdm(total=N_process_part, position=comm.rank)

    # axis_granuralrity = 0.5
    # particle_span_x = axis_granuralrity
    # particle_span_y = axis_granuralrity

# Definir o intervalo de períodos entre 5.9 e 7.7 e dividir em 10 bins
bins = np.linspace(5.9, 7.7, 101)


heatmap_data = np.zeros((len(bins) - 1, len(t)))


per_len = np.zeros((len(bins) - 1))
for i in range(len(bins) - 1):
    per_len[i] = bins[i + 1]/2. + bins[i]/2.


for b in tqdm(range(Per.shape[1])):
#for b in tqdm(range(2*n_sam_fl)):


    for a in range(Per.shape[0]):
        period = Per[a, b]/year

        for i in range(len(bins) - 1):
            if bins[i] <= period < bins[i + 1]:
                heatmap_data[i, b] += 1


#heatmap_data[heatmap_data==0.0] = np.nan

fig, ax1 = plt.subplots(figsize=(10, 6))
sns.heatmap(heatmap_data, cmap='viridis', norm=LogNorm(vmin=1, vmax=np.max(heatmap_data)), cbar_kws={'label': 'Particles'})
#sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': 'Particles'})


ax1.set_title('')
ax1.set_xlabel('Time (year)')
ax1.set_ylabel('T (year)')
#plt.ylim(-0.5, len(bins) - 1 + 0.5)

xticks = np.linspace(0, len(t) - 1, 20, dtype=int)
yticks = np.linspace(0, len(per_len) - 1, 10, dtype=int)
ytick_labels = np.round(per_len[yticks], 1)
ax1.set_yticks(yticks, ytick_labels)
ax1.invert_yaxis()
ax1.set_xticks(xticks, np.round(t[xticks]/year))
ax2 = ax1.twinx()

line_value = 7.11552
# Encontrando o índice mais próximo ao valor de 7.11552
closest_idx = np.argmin(np.abs(per_len - line_value))
ax2.axhline(y=closest_idx, color='r', linestyle='--')
#ax2.set_ylim([5.9,7.7])
ax2.set_ylim(ax1.get_ylim())
ax2.set_yticks([])
plt.savefig(plot_folder / f"map.png", bbox_inches="tight")



