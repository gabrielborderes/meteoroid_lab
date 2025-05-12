import rebound
import reboundx
import scipy.constants as constants
import os

import h5py
import numpy as np


from tqdm import tqdm

def record_coll(sim_pointer, collision):
    from os import path

    sim = sim_pointer.contents



    #with open(collisionFile,'a') as f:
        #f.write(f'{sim.t}-{sim.particles[collision.p1].hash}-{sim.particles[collision.p2].hash}\n')


    if (sim.particles[collision.p1].m > 0.0):
        ps = sim.particles[collision.p2]
        if (sim.particles[collision.p2].m > 0.0):
            ps2 = sim.particles[collision.p1]
            R = np.sqrt((ps.x - ps2.x)**2+(ps.y - ps2.y)**2+(ps.z - ps2.z)**2)
            #with open(collisionFile,'a') as f:
                #f.write(f'main body collision- distance- {R} {ps.a} {ps.e} {ps.inc} {ps.Omega} {ps.omega} {ps.f}  \n')

            return 0
        else:
            with open(f"data/{sim.particles[collision.p1].hash.value}.txt",'a') as f:
                f.write(f'{sim.t} {ps.a} {ps.e} {ps.inc} {ps.Omega} {ps.omega} {ps.f} {ps.r} \n')

            return 2
    else:
        with open(f"data/{sim.particles[collision.p2].hash.value}.txt",'a') as f:
            ps = sim.particles[collision.p1]
            f.write(f'{sim.t} {ps.a} {ps.e} {ps.inc} {ps.Omega} {ps.omega} {ps.f} {ps.r} \n')

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


from config import (
    data_folder,
    checkpoint,
    dictionary_file,
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
    factor_flush,
    srp_active

)

#LONG_ARCH = True
collisionFile = data_folder / 'closeRegister.txt'

if size <= 1:
    filename = data_folder / "cache_21P_all.bin"
    save_file = ephem_file

else:
    filename = data_folder / f"cache_21P_{rank}.bin"
    save_file = f"ephemerides_21P_{rank}.h5"

with open(collisionFile,'w') as f:
    f.write(f'\n')


with h5py.File(particle_file, 'r') as f:
    shell_part = f["shell_part"][:]
    #radius = f["radius"][:]
    N_activity = f["N_activity"][:]
    particle_velocities = f["particle_velocities"][:]


N_part = shell_part.shape[0]
nb = len(solar_system_objects)


if os.path.exists(checkpoint):
    print("checkpoint reboot not implemented yet")
    sys.exit()
else:
    sim = rebound.Simulation(str(init_sim_file))
    count_file = int(0)



with open(dictionary_file,'w') as f:
    for i in range(len(solar_system_objects)):
        f.write(f"{solar_system_objects[i]}, hash: {sim.particles[solar_system_objects[i]].hash}\n")



sim.move_to_com()
# Configure simulation
sim.integrator = "ias15"
# Let's save it for next time
# Note: sim.save_to_file() only saves the particle data, not the integrator settings, etc.


# Set the integration time
ps = sim.particles



if srp_active:
    rebx = reboundx.Extras(sim)
    rf = rebx.load_force("radiation_forces")
    rebx.add_force(rf)
    rf.params["c"] = 2.99792458e8
    ps["Sun"].params["radiation_source"] = 1




iter_inds = np.arange(rank, N_part, size, dtype=np.int64)
N_process_part = len(iter_inds)


C_per = 2.*np.pi/np.sqrt(sim.G)
M_sun = sim.particles["Sun"].m

flush = int(Nout/factor_flush)

x = np.zeros((sim.N + N_process_part, flush))
y = np.zeros((sim.N + N_process_part, flush))
z = np.zeros((sim.N + N_process_part, flush))
vx = np.zeros((sim.N + N_process_part, flush))
vy = np.zeros((sim.N + N_process_part, flush))
vz = np.zeros((sim.N + N_process_part, flush))

oe_a = np.zeros((sim.N + N_process_part, flush))
oe_e = np.zeros((sim.N + N_process_part, flush))
oe_i = np.zeros((sim.N + N_process_part, flush))
oe_O = np.zeros((sim.N + N_process_part, flush))
oe_o = np.zeros((sim.N + N_process_part, flush))
oe_f = np.zeros((sim.N + N_process_part, flush))
oe_M = np.zeros((sim.N + N_process_part, flush))

oe_Per = np.zeros((sim.N + N_process_part, flush))


id_part = 0

pbar = tqdm(total=len(sim_time), position=rank)

# print(solar_system_objects[-1])
#
# print(sim.particles[solar_system_objects[-1]])

count = int(0)
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

            # print(">>>>>",solar_system_objects[-1])
            # print(sim.particles[solar_system_objects[-1]])

            psc = sim.particles[solar_system_objects[-1]]
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
            ps[f"{id_part}"].r = shell_part[id_part, 0]
            with open(dictionary_file,'a') as f:
                f.write(f"{id_part}, hash: {sim.particles[f"{id_part}"].hash}\n")

            if srp_active:
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
        if bind != 0:
            oe_a[bind][count] = p.a
            oe_e[bind][count] = p.e
            oe_i[bind][count] = p.inc
            oe_O[bind][count] = p.Omega
            oe_o[bind][count] = p.omega
            oe_f[bind][count] = p.f
            oe_M[bind][count] = p.M
            oe_Per[bind][count] = C_per*np.sqrt(p.a**3 / (p.m + M_sun))

    for arr_ind, pind in enumerate(iter_inds):
        # if pind >= id_part:
        #     break

        try:
            p = ps[f"{pind}"]  # Tentando acessar a partícula
            p_exist = True
        except rebound.ParticleNotFound:
            p_exist = False


        # if f"{pind}" in ps:
        #     p = ps[f"{pind}"]
        # else:
        #     break
        if p_exist:
            x[nb + arr_ind][count] = p.x
            y[nb + arr_ind][count] = p.y
            z[nb + arr_ind][count] = p.z
            vx[nb + arr_ind][count] = p.vx
            vy[nb + arr_ind][count] = p.vy
            vz[nb + arr_ind][count] = p.vz
            oe_a[nb + arr_ind][count] = p.a
            oe_e[nb + arr_ind][count] = p.e
            oe_i[nb + arr_ind][count] = p.inc
            oe_O[nb + arr_ind][count] = p.Omega
            oe_o[nb + arr_ind][count] = p.omega
            oe_f[nb + arr_ind][count] = p.f
            oe_M[nb + arr_ind][count] = p.M
            oe_Per[nb + arr_ind][count] = C_per*np.sqrt(p.a**3 / (p.m + M_sun))



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
            hf.create_dataset("a", data=oe_a)
            hf.create_dataset("e", data=oe_e)
            hf.create_dataset("i", data=oe_i)
            hf.create_dataset("O", data=oe_O)
            hf.create_dataset("o", data=oe_o)
            hf.create_dataset("f", data=oe_f)
            hf.create_dataset("M", data=oe_M)
            hf.create_dataset("Per", data=oe_Per)

        count_file = int(count_file+1)




pbar.close()

if count != 0:
    sim.save_to_file(str(checkpoint))

    sv_folder = data_folder / f"part_{count_file}"
    sv_folder.mkdir(parents=True, exist_ok=True)

    file_name = sv_folder / save_file
    with h5py.File(str(file_name), "w") as hf:
        hf.create_dataset("index", data=iter_inds)
        hf.create_dataset("t", data=sim_time)
        hf.create_dataset("x", data=x)
        hf.create_dataset("y", data=y)
        hf.create_dataset("z", data=z)
        hf.create_dataset("vx", data=vx)
        hf.create_dataset("vy", data=vy)
        hf.create_dataset("vz", data=vz)
        hf.create_dataset("a", data=oe_a)
        hf.create_dataset("e", data=oe_e)
        hf.create_dataset("i", data=oe_i)
        hf.create_dataset("O", data=oe_O)
        hf.create_dataset("o", data=oe_o)
        hf.create_dataset("f", data=oe_a)
        hf.create_dataset("M", data=oe_a)
        hf.create_dataset("Per", data=oe_Per)

    count_file = int(count_file+1)



sim.save_to_file(str(filename))


