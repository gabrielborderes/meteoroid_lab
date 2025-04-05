import rebound
import reboundx
import scipy.constants as constants

import h5py
import numpy as np

import sys

import os
import re
import configparser

from tqdm import tqdm

def record_coll(sim_pointer, collision):
    from os import path

    sim = sim_pointer.contents

    ps = sim.particles

    if (code_sim[ps[collision.p1].hash.value] in solar_system_objects ):
        planet = collision.p1
        particle = collision.p2
    else:
        planet = collision.p2
        particle = collision.p1
    print(sim.t,'close enconteur ', code_sim[ps[planet].hash.value],code_sim[ps[particle].hash.value])

    osPart = ps[particle].orbit(primary=ps[planet])

    aC,eC,iC = osPart.a , osPart.e, osPart.inc

    rp =  aC*(1-eC)

    deltaP = ps[particle] - sim.particles[planet]
    deltaR = np.sqrt(deltaP.x**2 + deltaP.y**2 + deltaP.z**2)

    with open(collisionFile,'a') as f:
        f.write(f'{sim.t};{deltaR};{rp};{code_sim[ps[planet].hash.value]};{code_sim[ps[particle].hash.value]}\n')

    if (code_sim[sim.particles[collision.p1].hash.value] in solar_system_objects):
        return 2
    else:
        return 1






try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
except ImportError:
    rank = 0
    size = 1

    class COMM_WORLD:
        rank = 0
        size = 1

        def barrier(self):
            pass


from config_3200Phaethon import (
    data_folder,
    checkpoint,
    solar_luminosity,
    solar_system_objects,
    init_sim_file,
    particle_bulk_density,
    sim_time,
    Nout,
    body_radius,
    Eje_p_R,
    T_comet,
    RotVel,
    particle_file,
    ephem_file,
)

LONG_ARCH = True


if size <= 1:
    filename = data_folder / "cache_3200_Phaethon_all.bin"
    save_file = ephem_file
    collisionFile = data_folder / 'closeRegister.txt'

else:
    filename = data_folder / f"cache_3200_Phaethon_{rank}.bin"
    save_file = f"ephemerides_3200_Phaethon_{rank}.h5"
    collisionFile = data_folder / f'closeRegister_{rank}.txt'

with open(collisionFile,'w') as f:
    f.write(f'time;deltaR;rp;p1;p2\n')


with h5py.File(particle_file, 'r') as f:
    shell_part = f["shell_part"][:]
    N_activity = f["N_activity"][:]
    particle_velocities = f["particle_velocities"][:]
    key_sim = f["key_sim"][:]

code_sim = dict()

for i, key in enumerate(key_sim):
    code_sim[str(key)] = solar_system_objects[i]



N_part = shell_part.shape[0]
nb = len(solar_system_objects)


if os.path.exists(checkpoint):
    print("not implemented yet")
    sys.exit()

else:
    sim = rebound.Simulation(str(init_sim_file))
    count_file = int(0)

sim.move_to_com()
# Configure simulation
sim.integrator = "ias15"
# Let's save it for next time
# Note: sim.save_to_file() only saves the particle data, not the integrator settings, etc.
# sim.collision = "direct"
# sim.collision_resolve = record_coll

# Set the integration time
ps = sim.particles


rebx = reboundx.Extras(sim)
rf = rebx.load_force("radiation_forces")
rebx.add_force(rf)
rf.params["c"] = 2.99792458e8
ps["Sun"].params["radiation_source"] = 1




iter_inds = np.arange(rank, N_part, size, dtype=np.int64)
N_process_part = len(iter_inds)

flush = int(Nout/100)

x = np.zeros((sim.N + N_process_part, flush))
y = np.zeros((sim.N + N_process_part, flush))
z = np.zeros((sim.N + N_process_part, flush))
vx = np.zeros((sim.N + N_process_part, flush))
vy = np.zeros((sim.N + N_process_part, flush))
vz = np.zeros((sim.N + N_process_part, flush))
id_part = 0

count = int(0)
pbar = tqdm(total=len(sim_time), position=rank)
for ti, t in enumerate(sim_time):
    sim.integrate(t)
    pbar.update(1)

    if ti == N_activity.shape[0]:
        sim.collision = "direct"
        sim.collision_resolve = record_coll

    if ti < N_activity.shape[0]:
        pbar.set_description(f"added {N_activity[ti]} particles")
        for pi in range(N_activity[ti]):
            if id_part not in iter_inds:
                id_part += 1
                continue

            psc = ps[nb - 1]
            part_vel = particle_velocities[id_part]
            # Consider the velocity as rotating in plane of the inertial frame
            CRotVel_x = -RotVel * body_radius * shell_part[id_part, 2]
            CRotVel_y = RotVel * body_radius * shell_part[id_part, 1]
            #      = orbital pos. + Ejection pos. (outburst)
            part_x = psc.x + Eje_p_R * shell_part[id_part, 1]
            part_y = psc.y + Eje_p_R * shell_part[id_part, 2]
            part_z = psc.z + Eje_p_R * shell_part[id_part, 3]
            #       = orbital vel.  + Ejection vel. (outburst)       + velocity due to the comet rotation
            part_vx = psc.vx + part_vel * shell_part[id_part, 1] + CRotVel_x
            part_vy = psc.vy + part_vel * shell_part[id_part, 2] + CRotVel_y
            part_vz = psc.vz + part_vel * shell_part[id_part, 3]
            sim.add(
                x=part_x,
                y=part_y,
                z=part_z,
                vx=part_vx,
                vy=part_vy,
                vz=part_vz,
                hash=f"{id_part}",
            )

            grain_radius = shell_part[id_part, 0]  # grain radius in m
            Q_pr = 1.0
            density = particle_bulk_density
            # grain_radius = 1.e-5 # grain radius in m
            # density = 1000. # kg/m^3 = 1g/cc
            beta = rebx.rad_calc_beta(
                constants.G,
                rf.params["c"],
                ps[0].m,
                solar_luminosity,
                grain_radius,
                density,
                Q_pr,
            )
            ps[f"{id_part}"].params["beta"] = beta
            # = 3.*luminosity*Q_pr/(16.*np.pi*sim.G*1.e9*ps[0].m*rf.params["c"]*density*grain_radius)

            id_part += 1

    for bind in range(nb):
        p = ps[f"{solar_system_objects[bind]}"]
        x[bind][count] = p.x
        y[bind][count] = p.y
        z[bind][count] = p.z
        vx[bind][count] = p.vx
        vy[bind][count] = p.vy
        vz[bind][count] = p.vz

    for arr_ind, pind in enumerate(iter_inds):
        if pind >= id_part:
            break
        p = ps[f"{pind}"]
        x[nb + arr_ind][count] = p.x
        y[nb + arr_ind][count] = p.y
        z[nb + arr_ind][count] = p.z
        vx[nb + arr_ind][count] = p.vx
        vy[nb + arr_ind][count] = p.vy
        vz[nb + arr_ind][count] = p.vz


    count =  int(count+1)
    if count == flush:
        count = int(0)

        sim.save_to_file(str(checkpoint))

        sv_folder = data_folder / f"part_{count_file}"
        sv_folder.mkdir(parents=True, exist_ok=True)

        file_name = sv_folder / save_file
        # Save data
        with h5py.File(str(file_name), "w") as hf:
            hf.create_dataset("index", data=iter_inds)
            hf.create_dataset("t", data=sim_time)
            hf.create_dataset("x", data=x)
            hf.create_dataset("y", data=y)
            hf.create_dataset("z", data=z)
            hf.create_dataset("vx", data=vx)
            hf.create_dataset("vy", data=vy)
            hf.create_dataset("vz", data=vz)

        count_file = int(count_file+1)

pbar.close()
# Save data
if count != 0:
    sim.save_to_file(str(checkpoint))

    sv_folder = data_folder / f"part_{count_file}"
    sv_folder.mkdir(parents=True, exist_ok=True)

    file_name = sv_folder / save_file
    # Save data
    with h5py.File(str(file_name), "w") as hf:
        hf.create_dataset("index", data=iter_inds)
        hf.create_dataset("t", data=sim_time)
        hf.create_dataset("x", data=x)
        hf.create_dataset("y", data=y)
        hf.create_dataset("z", data=z)
        hf.create_dataset("vx", data=vx)
        hf.create_dataset("vy", data=vy)
        hf.create_dataset("vz", data=vz)




