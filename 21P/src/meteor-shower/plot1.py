# import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
import h5py
import pathlib
import rebound
import numpy as np
import scipy.constants as constants
import matplotlib.font_manager as font_manager
from tqdm import tqdm
from astropy import units
import os
from astropy.constants import GM_sun, au
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

rh = np.zeros((nb))

sim = rebound.Simulation()
sim.units = ("m", "s", "kg")
date = "2000-01-01 00:00"
for i in range(len(solar_system_objects)):
    sim.add(solar_system_objects[i], date=date, hash=solar_system_objects[i])

sim.move_to_com()
# Configure simulation
sim.integrator = "ias15"
# Let's save it for next time
# Note: sim.save_to_file() only saves the particle data, not the integrator settings, etc.
# sim.collision = "direct"
# sim.collision = "none"
ps = sim.particles

for n in range(1,nb):
    bd = ps[f"{solar_system_objects[n]}"]
    rh[n] = bd.a*(1-bd.e)*(bd.m/(3.*(bd.m+ps[f"{solar_system_objects[0]}"].m)))**(1./3.)

print(rh)


nb = len(solar_system_objects)
year = (units.year).to(units.second)

data_folder = pathlib.Path("data").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

particle_file = data_folder / "particles_21P.h5"


my_cmap = plt.get_cmap("viridis")


plot_folder = data_folder / f"21P_plot_dist"

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

for fl in range(factor_flush):
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

        dist = np.zeros((nb,len(t)))
        n_sam_fl = x.shape[1]

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

    for n in tqdm(range(nb)):
        for nt in range(x.shape[1]):
            min_d = np.sqrt((x[n][nt]-x[nb][nt])**2 + (y[n][nt]-y[nb][nt])**2 + (z[n][nt]-z[nb][nt])**2)
            for n_p in range(nb+1,x.shape[0]):
                val = np.sqrt((x[n][nt]-x[n_p][nt])**2 + (y[n][nt]-y[n_p][nt])**2 + (z[n][nt]-z[n_p][nt])**2)
                if min_d > val:
                    min_d = val

            dist[n][nt + n_sam_fl*fl] = min_d





for n in range(1,nb):
    fig, ax = plt.subplots(figsize=(7, 6))
    fig.suptitle(f"{solar_system_objects[n]}", fontsize=16)
    ax.plot(t/year, dist[n], color="blue")
    ax.axhline(y=rh[n], color='r', linestyle='--')
    ax.grid(True)
    ax.set_yscale('log')  # Definindo o eixo Y como escala logarítmica
    plt.xlabel("Time (years)")
    plt.ylabel(r"Minimum distance (m)")
    plt.tight_layout()
    figure_name = plot_folder / f"dist_{solar_system_objects[n]}.png"
    plt.savefig(figure_name, bbox_inches="tight")
    ax.clear()































