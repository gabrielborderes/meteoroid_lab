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


n_file = 0

p = int(3)
q = int(2)

#for fl in range(factor_flush):
for fl in range(2):
    print(f"Part_{fl} running")
    ephem_file = data_folder / f"part_{fl}/ephemerides_21P_all.h5"

    if first:
        first = False
        with h5py.File(ephem_file, "r") as hf:
            t = np.concatenate([
                hf["t"][()],
            ], axis=0)
            O = np.concatenate([
                hf["body_y"][()],
                hf["O"][()],
            ], axis=0)
            o = np.concatenate([
                hf["body_z"][()],
                hf["o"][()],
            ], axis=0)
            M = np.concatenate([
                hf["body_Per"][()],
                hf["M"][()],
            ], axis=0)

        phi = np.zeros((O.shape[0],len(t)))
        n_sam_fl = O.shape[1]

    else:
        with h5py.File(ephem_file, "r") as hf:
            O = np.concatenate([
                hf["body_y"][()],
                hf["O"][()],
            ], axis=0)
            o = np.concatenate([
                hf["body_z"][()],
                hf["o"][()],
            ], axis=0)
            M = np.concatenate([
                hf["body_Per"][()],
                hf["M"][()],
            ], axis=0)

    for nt in tqdm(range(O.shape[1])):
        l_J = O[5,nt] + o[5,nt] + M[5,nt]
        for n in range(nb,O.shape[0]):
            l_bd = O[n,nt] + o[n,nt] + M[n,nt]
            varpi = O[n,nt] + o[n,nt]
            phi0 = ( p + q ) * l_J - p * l_bd - q * varpi
            phi[n][nt] = np.arctan2(np.sin(phi0),np.cos(phi0))






fig, ax = plt.subplots(figsize=(7, 6))
for n in range(nb,O.shape[0]):
    dphi_dt = np.gradient(phi[n,:], t)
    ax.plot(t/year, dphi_dt, color="blue")

ax.grid(True)
ax.set_yscale('log')  # Definindo o eixo Y como escala logarítmica
plt.xlabel("Time (years)")
plt.ylabel(r"$\dv{}{t}\phi$")
plt.tight_layout()
figure_name = plot_folder / f"dif.png"
plt.savefig(figure_name, bbox_inches="tight")
ax.clear()































