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
)
nb = len(solar_system_objects)
year = (units.year).to(units.second)

data_folder = pathlib.Path("data_tmp").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

particle_file = data_folder / "particles_21P.h5"
ephem_file = data_folder / f"ephemerides_21P_all.h5"


# sim = rebound.Simulation()
# sim.units = ("m", "s", "kg")
# date = "2000-01-01 00:00"
# for i in range(len(solar_system_objects)):
#     sim.add(solar_system_objects[i], date=date, hash=solar_system_objects[i])
#
# sim.move_to_com()
# Configure simulation
# sim.integrator = "ias15"
# Let's save it for next time
# Note: sim.save_to_file() only saves the particle data, not the integrator settings, etc.
# sim.collision = "direct"
# sim.collision = "none"
# ps = sim.particles


extra_sim = np.arange(0,59900.*24.*3600.,24.*3600.)

x_bodies = np.zeros((nb, len(extra_sim)))
y_bodies = np.zeros((nb, len(extra_sim)))




# for ti, t in tqdm(enumerate(extra_sim), total=len(extra_sim)):
#     sim.integrate(t)
#
#     for bind in range(nb):
#         p = ps[f"{solar_system_objects[bind]}"]
#         x_bodies[bind][ti] = p.x
#         y_bodies[bind][ti] = p.y
#




OB_id = 1
part_i = 24


my_cmap = plt.get_cmap("viridis")


plot_folder = data_folder / f"21P_plot/ress_anfgle"
#if comm.rank == 0:
if not plot_folder.is_dir():
    plot_folder.mkdir()

#comm.barrier()



with h5py.File(particle_file, 'r') as f:
    shell_part = f["shell_part"][:]
    N_activity = f["N_activity"][:]




N_part = shell_part.shape[0]
nb = len(solar_system_objects)




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
    Per = np.concatenate([
        hf["body_Per"][()],
        hf["Per"][()],
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
    e = np.concatenate([
        hf["body_Per"][()],
        hf["e"][()],
    ], axis=0)


dist = np.sqrt(x**2+y**2+z**2)
d_mar_comet = np.sqrt((x[4,:]-x[9,:])**2+(y[4,:]-y[9,:])**2+(z[4,:]-z[9,:])**2)

e[e==0.0] = np.nan

x[x==0.0] = np.nan
y[y==0.0] = np.nan
z[z==0.0] = np.nan
Per[Per==0.0] = np.nan
dist[dist==0.0] = np.nan
# O[e>1.] = np.nan
# o[e>1.] = np.nan
# M[e>1.] = np.nan

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
n_file = 0





fig = plt.figure()
fig.set_figheight(14)
fig.set_figwidth(21)

ax1 = plt.subplot2grid(shape=(2, 1), loc=(0, 0))
ax2 = plt.subplot2grid(shape=(2, 1), loc=(1, 0))

ax1.plot(t/year, d_mar_comet/au.value, color="blue")
ax1.set_xlabel("time (year)")
ax1.set_ylabel("d(mars)x(21P/Giacobini–Zinner) (au)")
ax1.set_yscale('log')
ax1.set_xlim([0,200])


ax2.plot(t/year, Per[nb-1]/year, color="blue", label=f"21P/Giacobini–Zinner")
ax2.axhline(y=7.11552, color='r', linestyle='--')
ax2.set_xlabel("time (year)")
ax2.set_ylabel("T (year)", color='blue')
ax3 = ax2.twinx()
ax3.plot(t/year, dist[nb-1]/au.value, color="green", label="d (au)")
ax3.set_ylabel("d (au)", color='green')
ax2.set_xlim([0,200])
ax2.set_ylim([5.9,7.7])



fig.savefig( f"mars_comet.png", bbox_inches="tight")

p = int(3)
q = int(2)

l_J = O[5,:] + o[5,:] + M[5,:]

for ii in tqdm(range(len(x[nb:, 0]))):
    i = ii + nb
    fig = plt.figure()
    fig.set_figheight(14)
    fig.set_figwidth(21)

    ax1 = plt.subplot2grid(shape=(4, 1), loc=(0, 0))
    ax2 = plt.subplot2grid(shape=(4, 1), loc=(1, 0))
    ax3 = plt.subplot2grid(shape=(4, 1), loc=(2, 0))
    ax4 = plt.subplot2grid(shape=(4, 1), loc=(3, 0))

    l_bd = O[i,:] + o[i,:] + M[i,:]
    varpi = O[i,:] + o[i,:]
    phi0 = ( p + q ) * l_J - p * l_bd - q * varpi
    phi = np.arctan2(np.sin(phi0),np.cos(phi0))
    d_earth = np.sqrt((x[3,:]-x[i,:])**2+(y[3,:]-y[i,:])**2+(z[3,:]-z[i,:])**2)

    ax1.plot(t/year, np.degrees(phi), color="blue")
    ax1.set_ylim([-180,180])
    ax1.set_xlim([0,200])
    ax1.set_xlabel("time (year)")
    ax1.set_ylabel(r"$\varphi$ (degrees)")
    ax2.plot(t/year, d_earth/au.value, color="blue")
    ax2.set_yscale('log')
    ax2.set_xlim([0,200])
    ax2.set_xlabel("time (year)")
    ax2.set_ylabel(r"$d_{Earth}$ (au)")
    ax3.plot(t/year, Per[i,:]/year, color="blue")
    ax3.set_xlabel("time (year)")
    ax3.set_ylabel("T (year)")
    ax3.set_xlim([0,200])
    ax4.plot(t/year, e[i,:], color="blue")
    ax4.set_ylabel("e")
    ax4.set_xlabel("time (year)")
    ax4.set_xlim([0,200])
    fig.savefig(plot_folder / f"t_int{ii}.png", bbox_inches="tight")
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()





















