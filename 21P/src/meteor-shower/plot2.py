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


from config_l import (
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


extra_sim = np.arange(0,59900.*24.*3600.,24.*3600.)

x_bodies = np.zeros((nb, len(extra_sim)))
y_bodies = np.zeros((nb, len(extra_sim)))




for ti, t in tqdm(enumerate(extra_sim), total=len(extra_sim)):
    sim.integrate(t)

    for bind in range(nb):
        p = ps[f"{solar_system_objects[bind]}"]
        x_bodies[bind][ti] = p.x
        y_bodies[bind][ti] = p.y





my_cmap = plt.get_cmap("viridis")


plot_folder = data_folder / f"21P_plot"
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


    for ti in tqdm(range(0,iter_inds,30)):
        fig = plt.figure()
        fig.set_figheight(14)
        fig.set_figwidth(21)


        ax1 = plt.subplot2grid(shape=(2, 3), loc=(0, 0))
        ax2 = plt.subplot2grid(shape=(2, 3), loc=(0, 1))
        ax3 = plt.subplot2grid(shape=(2, 3), loc=(0, 2))
        ax4 = plt.subplot2grid(shape=(2, 3), loc=(1, 0), colspan=3)
        ############################################################################
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
        ax4.plot(t[:Per.shape[1]]/year, Per[nb-1]/year, color="blue", label=f"21P/Giacobini–Zinner")
        ax4.scatter(np.ones_like(Per[nb:,ti])*t[fl*n_sam_fl + ti]/year, Per[nb:,fl*n_sam_fl + ti]/year,  c=np.log10(shell_part[:, 0]), cmap=my_cmap)
        ax4.scatter(t[fl*n_sam_fl + ti]/year, Per[nb-1,fl*n_sam_fl + ti]/year,  c="red")
        #ax4.scatter(np.ones_like(Per_n[nb:,ti])*t[fl*n_sam_fl + ti]/year, Per_n[nb:, ti]/year,  c=np.log10(shell_part[:, 0]), cmap=my_cmap)
        #ax4.scatter(np.ones_like(Per[:,ti])*t[fl*n_sam_fl + ti]/year, Per[:,fl*n_sam_fl + ti]/year,  c="gray")
        ax4.axhline(y=7.11552, color='r', linestyle='--')
        #ax4.set_xlim([ti, tf])
        ax4.set_xlabel("time (year)")
        ax4.set_ylabel("T (year)", color='blue')

        #ax5 = ax4.twinx()
        #ax5.plot(t[:Per.shape[1]]/year, dist[nb-1]/au.value, color="green", label="d (au)")
        #ax5.set_ylabel("d (au)", color='green')

        tt = 49.5+t[ti]/(3600.*24.*365.25)
        ax1.set_xlabel("x (AU)", fontsize=sizefont)
        ax1.set_ylabel("y (AU)", fontsize=sizefont)
        #ax1.text(0.14,0.13, f"time = {tt:.2f} years", fontproperties=font_manager.FontProperties(size=14))
        ax4.set_title(f"time = {tt:.2f} years")
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
            x=x[:nb, ti] / au.value, y=y[:nb, ti] / au.value, c="grey"
        )
        ax2.scatter(
            x=x[:nb, ti] / au.value, y=y[:nb, ti] / au.value, c="grey"
        )
        ax3.scatter(
            x=x[:nb, ti] / au.value, y=y[:nb, ti] / au.value, c="grey"
        )

        # ax1.scatter(
        #     x=x[nb-1, 0] / au.value, y=y[nb-1, 0] / au.value,  c="grey"
        # )
        # ax2.scatter(
        #     x=x[nb-1, 0] / au.value, y=y[nb-1, 0] / au.value,  c="grey"
        # )
        # ax3.scatter(
        #     x=x[nb-1, 0] / au.value, y=y[nb-1, 0] / au.value, c="grey"
        # )

        ax1.scatter(
            x[nb:, ti] / au.value,
            y[nb:, ti] / au.value,
            s=2.0,
            c=np.log10(shell_part[:, 0]),
            cmap=my_cmap,
        )
        im = ax2.scatter(
            x[nb:, ti] / au.value,
            y[nb:, ti] / au.value,
            s=2.0,
            c=np.log10(shell_part[:, 0]),
            cmap=my_cmap,
        )
        im = ax3.scatter(
            x[nb:, ti] / au.value,
            y[nb:, ti] / au.value,
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
        ax1.plot(x_bodies[1][0:Nog] / au.value, y_bodies[1][0:Nog] / au.value)
        ax2.plot(x_bodies[1][0:Nog] / au.value, y_bodies[1][0:Nog] / au.value)
        ax3.plot(x_bodies[1][0:Nog] / au.value, y_bodies[1][0:Nog] / au.value)


        Nog = 230
        ax1.plot(x_bodies[2][0:Nog] / au.value, y_bodies[2][0:Nog] / au.value)
        ax2.plot(x_bodies[2][0:Nog] / au.value, y_bodies[2][0:Nog] / au.value)
        ax3.plot(x_bodies[2][0:Nog] / au.value, y_bodies[2][0:Nog] / au.value)



        Nog = 370
        ax1.plot(x_bodies[3][0:Nog] / au.value, y_bodies[3][0:Nog] / au.value)
        ax2.plot(x_bodies[3][0:Nog] / au.value, y_bodies[3][0:Nog] / au.value)
        ax3.plot(x_bodies[3][0:Nog] / au.value, y_bodies[3][0:Nog] / au.value)



        Nog = 700
        ax1.plot(x_bodies[4][0:Nog] / au.value, y_bodies[4][0:Nog] / au.value)
        ax2.plot(x_bodies[4][0:Nog] / au.value, y_bodies[4][0:Nog] / au.value)
        ax3.plot(x_bodies[4][0:Nog] / au.value, y_bodies[4][0:Nog] / au.value)


        Nog = 4400
        ax1.plot(x_bodies[5][0:Nog] / au.value, y_bodies[5][0:Nog] / au.value)
        ax2.plot(x_bodies[5][0:Nog] / au.value, y_bodies[5][0:Nog] / au.value)
        ax3.plot(x_bodies[5][0:Nog] / au.value, y_bodies[5][0:Nog] / au.value)


        Nog = 11000
        ax1.plot(x_bodies[6][0:Nog] / au.value, y_bodies[6][0:Nog] / au.value)
        ax2.plot(x_bodies[6][0:Nog] / au.value, y_bodies[6][0:Nog] / au.value)
        ax3.plot(x_bodies[6][0:Nog] / au.value, y_bodies[6][0:Nog] / au.value)

        Nog = 3000
        #ax1.plot(x_bodies[9][0:Nog] / au.value, y_bodies[9][0:Nog] / au.value)
        ax2.plot(x_bodies[9][0:Nog] / au.value, y_bodies[9][0:Nog] / au.value)
        ax3.plot(x_bodies[9][0:Nog] / au.value, y_bodies[9][0:Nog] / au.value)


        for k in range(0, nb):
            if (
                x[k][-1] / au.value * x[k][-1] / au.value + y[k][-1] / au.value * y[k][-1] / au.value > 0.01 * zoom * zoom or k == 0
            ):
                ax1.annotate(
                    solar_system_objects[k],
                    xy=(x[k][ti] / au.value, y[k][ti] / au.value),
                    horizontalalignment="center",
                )
            if k > 4:
                ax2.annotate(
                    solar_system_objects[k],
                    xy=(x[k][ti] / au.value, y[k][ti] / au.value),
                    horizontalalignment="center",
                )
                ax3.annotate(
                    solar_system_objects[k],
                    xy=(x[k][ti] / au.value, y[k][ti] / au.value),
                    horizontalalignment="center",
                )

        if cbar is None:
            cbar = fig.colorbar(im, ax=ax3)
            cbar.set_label("$r_{part}=10^{x}$ m", fontsize=sizefont)
            im.set_clim(vmin=min_size_log, vmax=max_size_log)

        ax2.set_xlim([-1.5 * zoom, 1.5 * zoom])
        ax2.set_ylim([-1.5 * zoom, 1.5 * zoom])
        ax3.set_xlim([-3.5 * zoom, 3.5 * zoom])
        ax4.set_xlim([0,t[-1]/year])
        ax4.set_ylim([5.9,7.7])
        ax3.set_ylim([-3.5 * zoom, 3.5 * zoom])
        ax1.set_xlim([
            x[nb - 1, ti] / au.value - particle_span,
            x[nb - 1, ti] / au.value + particle_span,
        ])
        ax1.set_ylim([
            y[nb - 1, ti] / au.value - particle_span,
            y[nb - 1, ti] / au.value + particle_span,
        ])

        fig.savefig(plot_folder / f"t_int{n_file}.png", bbox_inches="tight")

        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        #ax5.clear()
        plt.close(fig)
        cbar = None

        n_file = n_file + 1


#pbar.close()

