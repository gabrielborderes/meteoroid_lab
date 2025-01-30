# import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
import h5py

import rebound
import numpy as np
import scipy.constants as constants
import matplotlib.font_manager as font_manager
from tqdm import tqdm

#try:
#    from mpi4py import MPI

#    comm = MPI.COMM_WORLD
#except ImportError:
#    class COMM_WORLD:
#        rank = 0
#        size = 1

#        def barrier(self):
#            pass


from config_3200Phaethon import (
    data_folder,
    solar_system_objects,
    particle_file,
    ephem_file,
    Nog,
    init_sim_file,
)
nb = len(solar_system_objects)



sim = rebound.Simulation()
sim.units = ("km", "s", "kg")
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


extra_sim = np.arange(0,59900.*24.*3600.,24.*3600.)

x_bodies = np.zeros((nb, len(extra_sim)))
y_bodies = np.zeros((nb, len(extra_sim)))




for ti, t in tqdm(enumerate(extra_sim), total=len(extra_sim)):
    sim.integrate(t)

    for bind in range(nb):
        p = ps[f"{solar_system_objects[bind]}"]
        x_bodies[bind][ti] = p.x
        y_bodies[bind][ti] = p.y








au = constants.au*1e-3


OB_id = 1
part_i = 24


my_cmap = plt.get_cmap("viridis")


plot_folder = data_folder / f"3200_Phaethon_plot"
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

print(x.shape)

size = shell_part[:, 0]

fig = plt.figure()
fig.set_figheight(7)
fig.set_figwidth(21)

ax1 = plt.subplot2grid(shape=(1, 3), loc=(0, 0))  # Sphere
ax2 = plt.subplot2grid(shape=(1, 3), loc=(0, 1))  # Orbit
ax3 = plt.subplot2grid(shape=(1, 3), loc=(0, 2))
############################################################################
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


for ti in tqdm(range(iter_inds)):
    #pbar.update(1)
    zoom = 1.8
    # dx = x[nb:, ti] - x[nb - 1, ti]
    # dy = y[nb:, ti] - y[nb - 1, ti]
    # new_particle_span_x = np.max(np.abs(dx[np.logical_not(np.isnan(dx))])) / au
    # new_particle_span_x = np.floor(new_particle_span_x/axis_granuralrity)*axis_granuralrity
    # new_particle_span_y = np.max(np.abs(dy[np.logical_not(np.isnan(dy))])) / au
    # new_particle_span_y = np.floor(new_particle_span_y/axis_granuralrity)*axis_granuralrity
    # if new_particle_span_x > particle_span_x:
    #     particle_span_x = new_particle_span_x
    # if new_particle_span_y > particle_span_y:
    #     particle_span_y = new_particle_span_y
    # particle_span = np.max([particle_span_x, particle_span_y])
    particle_span = 0.1
    #tt = 2007.843836+t[ti]/(3600.*24.*365.25)
    tt = 49.5+t[ti]/(3600.*24.*365.25)
    ax1.set_xlabel("x (AU)", fontsize=sizefont)
    ax1.set_ylabel("y (AU)", fontsize=sizefont)
    #ax1.text(0.14,0.13, f"time = {tt:.2f} years", fontproperties=font_manager.FontProperties(size=14))
    plt.title(f"time = {tt:.2f} years")
    ax2.set_xlabel("x (AU)", fontsize=sizefont)
    ax2.set_ylabel("y (AU)", fontsize=sizefont)
    ax3.set_xlabel("x (AU)", fontsize=sizefont)
    ax3.set_ylabel("y (AU)", fontsize=sizefont)
    #n_part = total_bodies - nb

    ax1.xaxis.set_tick_params(labelsize=sizefont)
    ax1.yaxis.set_tick_params(labelsize=sizefont)
    ax2.xaxis.set_tick_params(labelsize=sizefont)
    ax2.yaxis.set_tick_params(labelsize=sizefont)
    ax3.xaxis.set_tick_params(labelsize=sizefont)
    ax3.yaxis.set_tick_params(labelsize=sizefont)
    # ax2.yticks(fontsize=sizefont)
    # ax1.set_aspect("equal", "box")
    ax2.set_aspect("equal", "box")
    ax3.set_aspect("equal", "box")
    ax1.scatter(
        x=x[:nb, ti] / au, y=y[:nb, ti] / au, c="grey"
    )
    ax2.scatter(
        x=x[:nb, ti] / au, y=y[:nb, ti] / au, c="grey"
    )
    ax3.scatter(
        x=x[:nb, ti] / au, y=y[:nb, ti] / au, c="grey"
    )
    
    ax1.scatter(
        x=x[nb-1, 0] / au, y=y[nb-1, 0] / au,  c="grey"
    )
    ax2.scatter(
        x=x[nb-1, 0] / au, y=y[nb-1, 0] / au,  c="grey"
    )
    ax3.scatter(
        x=x[nb-1, 0] / au, y=y[nb-1, 0] / au, c="grey"
    )
    
    ax1.scatter(
        x[nb:, ti] / au,
        y[nb:, ti] / au,
        s=2.0,
        c=np.log10(shell_part[:, 0]),
        cmap=my_cmap,
    )
    im = ax2.scatter(
        x[nb:, ti] / au,
        y[nb:, ti] / au,
        s=2.0,
        c=np.log10(shell_part[:, 0]),
        cmap=my_cmap,
    )
    im = ax3.scatter(
        x[nb:, ti] / au,
        y[nb:, ti] / au,
        s=2.0,
        c=np.log10(shell_part[:, 0]),
        cmap=my_cmap,
    )
    
    # if (
    #     x[nb-1][0] / au * x[nb-1][0] / au + y[nb-1][0] / au * y[nb-1][0] / au > 0.01 * zoom * zoom or k == 0
    # ):
    #     ax1.annotate(
    #         "Outburst",
    #         xy=(x[nb-1][0] / au, y[nb-1][0] / au),
    #         horizontalalignment="center",
    #     )
    # ax2.annotate(
    #     "Outburst",
    #     xy=(x[nb-1][0] / au, y[nb-1][0] / au),
    #     horizontalalignment="center",
    # )
    # ax3.annotate(
    #     "Outburst",
    #     xy=(x[nb-1][0] / au, y[nb-1][0] / au),
    #     horizontalalignment="center",
    # )

    Nog = 90
    ax1.plot(x_bodies[1][0:Nog] / au, y_bodies[1][0:Nog] / au)
    ax2.plot(x_bodies[1][0:Nog] / au, y_bodies[1][0:Nog] / au)
    ax3.plot(x_bodies[1][0:Nog] / au, y_bodies[1][0:Nog] / au)


    Nog = 230
    ax1.plot(x_bodies[2][0:Nog] / au, y_bodies[2][0:Nog] / au)
    ax2.plot(x_bodies[2][0:Nog] / au, y_bodies[2][0:Nog] / au)
    ax3.plot(x_bodies[2][0:Nog] / au, y_bodies[2][0:Nog] / au)



    Nog = 370
    ax1.plot(x_bodies[3][0:Nog] / au, y_bodies[3][0:Nog] / au)
    ax2.plot(x_bodies[3][0:Nog] / au, y_bodies[3][0:Nog] / au)
    ax3.plot(x_bodies[3][0:Nog] / au, y_bodies[3][0:Nog] / au)



    Nog = 700
    ax1.plot(x_bodies[5][0:Nog] / au, y_bodies[5][0:Nog] / au)
    ax2.plot(x_bodies[5][0:Nog] / au, y_bodies[5][0:Nog] / au)
    ax3.plot(x_bodies[5][0:Nog] / au, y_bodies[5][0:Nog] / au)


    Nog = 4400
    ax1.plot(x_bodies[6][0:Nog] / au, y_bodies[6][0:Nog] / au)
    ax2.plot(x_bodies[6][0:Nog] / au, y_bodies[6][0:Nog] / au)
    ax3.plot(x_bodies[6][0:Nog] / au, y_bodies[6][0:Nog] / au)


    Nog = 11000
    ax1.plot(x_bodies[7][0:Nog] / au, y_bodies[7][0:Nog] / au)
    ax2.plot(x_bodies[7][0:Nog] / au, y_bodies[7][0:Nog] / au)
    ax3.plot(x_bodies[7][0:Nog] / au, y_bodies[7][0:Nog] / au)

    Nog = 30600
    ax1.plot(x_bodies[8][0:Nog] / au, y_bodies[8][0:Nog] / au)
    ax2.plot(x_bodies[8][0:Nog] / au, y_bodies[8][0:Nog] / au)
    ax3.plot(x_bodies[8][0:Nog] / au, y_bodies[8][0:Nog] / au)

    Nog = 59900
    ax1.plot(x_bodies[9][0:Nog] / au, y_bodies[9][0:Nog] / au)
    ax2.plot(x_bodies[9][0:Nog] / au, y_bodies[9][0:Nog] / au)
    ax3.plot(x_bodies[9][0:Nog] / au, y_bodies[9][0:Nog] / au)




    Nog = 520
    ax1.plot(x_bodies[10][0:Nog] / au, y_bodies[10][0:Nog] / au)
    ax2.plot(x_bodies[10][0:Nog] / au, y_bodies[10][0:Nog] / au)
    ax3.plot(x_bodies[10][0:Nog] / au, y_bodies[10][0:Nog] / au)


    for k in range(0, nb):
        if (
            x[k][-1] / au * x[k][-1] / au + y[k][-1] / au * y[k][-1] / au > 0.01 * zoom * zoom or k == 0
        ):
            ax1.annotate(
                solar_system_objects[k],
                xy=(x[k][ti] / au, y[k][ti] / au),
                horizontalalignment="center",
            )
        if k > 4:
            ax2.annotate(
                solar_system_objects[k],
                xy=(x[k][ti] / au, y[k][ti] / au),
                horizontalalignment="center",
            )
            ax3.annotate(
                solar_system_objects[k],
                xy=(x[k][ti] / au, y[k][ti] / au),
                horizontalalignment="center",
            )

    if cbar is None:
        cbar = fig.colorbar(im, ax=ax3)
        cbar.set_label("$r_{part}=10^{x}$ m", fontsize=sizefont)
        im.set_clim(vmin=-6, vmax=-2)
    ax2.set_xlim([-1.5 * zoom, 1.5 * zoom])
    ax2.set_ylim([-1.5 * zoom, 1.5 * zoom])
    ax3.set_xlim([-3.5 * zoom, 3.5 * zoom])
    ax3.set_ylim([-3.5 * zoom, 3.5 * zoom])
    ax1.set_xlim([
        x[nb - 1, ti] / au - particle_span,
        x[nb - 1, ti] / au + particle_span,
    ])
    ax1.set_ylim([
        y[nb - 1, ti] / au - particle_span,
        y[nb - 1, ti] / au + particle_span,
    ])
    
    fig.savefig(plot_folder / f"t_int{ti}.png", bbox_inches="tight")

    ax1.clear()
    ax2.clear()
    ax3.clear()
#pbar.close()
plt.close(fig)
