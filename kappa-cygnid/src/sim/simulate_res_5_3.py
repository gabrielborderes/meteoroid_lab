import rebound
from tqdm import tqdm
import numpy as np
import h5py
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import os

import re
import configparser

from astropy.constants import GM_sun, au
from astropy import units

year = (units.year).to(units.second)
day = (units.day).to(units.second)

def record_coll(sim_pointer, collision):
    from os import path

    sim = sim_pointer.contents
    with open(file_coll,'a') as f:
        f.write(f'Collision at {sim.t} between {sim.particles[collision.p1].hash} and {sim.particles[collision.p2].hash}\n')
    return 0                             # Don't remove either particle


def read_meteors(file_name):
    with open(file_name, "r") as file:
        lines = file.readlines()

    header = re.split(r'\s+', lines[0].strip())
    data = {col: [] for col in header}

    for line in lines[1:]:
        values = re.split(r'\s+', line.strip())
        for i, value in enumerate(values):
            data[header[i]].append(value)

    # Lista de colunas a serem convertidas para float
    float_columns = {
        "JD",
        "AXIS",
        "AXISER",
        "ECC",
        "ECCERR",
        "INCL",
        "INCLER",
        "NODE",
        "NODEER",
        "ARGUP",
        "ARGUPER",
    }
    for col in float_columns:
        if col in data:  # Evita erro caso a coluna não exista
            data[col] = [float(x) for x in data[col]]

    return data

def OE_2_position(a, e, i, Omega, omega, f):

    r = a * (1 - e**2) / (1 + e * np.cos(f))
    varpi = omega + f

    x = r * (np.cos(Omega) * np.cos(varpi) - np.sin(Omega) * np.sin(varpi) * np.cos(i))
    y = r * (np.sin(Omega) * np.cos(varpi) + np.cos(Omega) * np.sin(varpi) * np.cos(i))
    z = r * (np.sin(varpi) * np.sin(i))

    return np.array([x, y, z])


def func_dist(f, a, e, i, Omega, omega, EP):
    POO = OE_2_position(a, e, i, Omega, omega, f)
    dist = np.sqrt((POO[0] - EP[0]) ** 2 + (POO[1] - EP[1]) ** 2 + (POO[2] - EP[2]) ** 2)

    return dist



config = configparser.ConfigParser()
#config.read("param.config")

config.read("param_ress.config")

# READ SIMULATION PARAMETERS
file_in_meteor = config["system"]["file_in_meteor"]
file_coll = config["system"]["file_coll"]
save_file = config["system"]["save_file"]
checkpoint = config["system"]["checkpoint"]

sim_time = config["sim_param"].getfloat("sim_time")
sim_step = config["sim_param"].getfloat("sim_step")

major_bodies = config["sim_param"]["major_bodies"].split(", ")
minor_bodies = config["sim_param"]["minor_bodies"].split(", ")

active_cl = config["clones"].getboolean("active")
n_clones = config["clones"].getint("n_clones")
seed_rd = config["clones"].getint("seed_rd")
np.random.seed(seed_rd)

index = []


# READ THE METEORS SAMPLE DATA
data = read_meteors(file_in_meteor)
n_G1 = data["GROUP"].count("G1")
n_G1A = data["GROUP"].count("G1A")



if os.path.exists(checkpoint):
    sim = rebound.Simulation(str(checkpoint))
    for i in range(len(major_bodies)):
        index.append(major_bodies[i])

    for i in range(len(minor_bodies)):
        index.append(minor_bodies[i])

    sor_ind = sorted(range(len(data['JD'])), key=lambda i: data['JD'][i], reverse=True)

    for i in tqdm(sor_ind):

        hash_meteor = f"{data['CODE'][i]}"
        index.append(hash_meteor)
        if active_cl:
            for j in range(n_clones):
                    hash_cl = f"{data['CODE'][i]}_cl_{j}"
                    index.append(hash_cl)


    time_array = np.arange(sim.t,-1.*sim_time*year,-1.*sim_step*year)
    dir, name= os.path.split(save_file)

    padrao = re.compile(rf"^{re.escape(name)}_(\d+)\.h5$")

    count_files_list = []


    for root, dirs, files in os.walk(dir):
        for arquivo in files:
            match = padrao.match(arquivo)
            if match:
                count_files_list.append(int(match.group(1)))

    if count_files_list:
        count_file = int(max(count_files_list) + 1)
    else:
        count_file = int(0)


else:
    # Starting simulation
    sim = rebound.Simulation()
    sim.units = ("m", "s", "kg")

    # Starting the system on the latest meteor time
    date_ini = f"JD{max(data['JD']):.10f}"

    #kappa_cygnid = tuple(map(float, config["sim_param"]["kappa_cygnid"].split(", ")))

    with open(file_coll,'w') as f:
        f.write(f'Dictionary\n')

    for i in range(len(major_bodies)):
        sim.add(
            major_bodies[i],
            date = date_ini,
            hash = major_bodies[i]
        )
        index.append(major_bodies[i])
        with open(file_coll,'a') as f:
            f.write(f'Body : {major_bodies[i]}, hash {sim.particles[major_bodies[i]].hash}\n')

    for i in range(len(minor_bodies)):
        sim.add(
            minor_bodies[i],
            date = date_ini,
            hash = minor_bodies[i]
        )
        index.append(minor_bodies[i])
        with open(file_coll,'a') as f:
            f.write(f'Body : {minor_bodies[i]}, hash {sim.particles[minor_bodies[i]].hash}\n')

    # SORT THE METEOR SEQUENCE FROM THE LATEST (MOST RECENT) TO THE OLDEST
    sor_ind = sorted(range(len(data['JD'])), key=lambda i: data['JD'][i], reverse=True)

    sim.move_to_com()
    sim.integrator = "ias15"
    # sim.collision = "direct"
    # sim.collision_resolve = record_coll

    # ADD METEORS IN THE SIMULATION
    for i in tqdm(sor_ind):

        # IF IS NOT THE MOST RECENT, IT EVOLVE THE SYSTEM
        if i > 0:
            new_time_stop = (data['JD'][i] - data['JD'][0]) * day
            sim.integrate(new_time_stop)

        #sma = data['AXIS'][i] * au.value
        sma = sim.particles['Jupiter'].a*(3./5.)**(2./3.)

        if sma > (data['AXIS'][i] - data['AXISER'][i]) * au.value:
            if sma < (data['AXIS'][i] + data['AXISER'][i]) * au.value:
                ecc = data['ECC'][i]
                incl = np.radians(data['INCL'][i])
                Ome = np.radians(data['NODE'][i])
                ome = np.radians(data['ARGUP'][i])

                # FIND EARTH POSSITION
                pe = sim.particles["399"]
                EP = np.array([pe.x, pe.y, pe.z])

                # FIND TRUE ANOMALY
                seek_min = minimize(
                    func_dist,
                    x0=np.radians(0.0e0),
                    args=
                        (
                        sma,
                        ecc,
                        incl,
                        Ome,
                        ome,
                        EP,
                        ),
                    bounds=[(-2.*np.pi, 2.*np.pi)],
                    method='Nelder-Mead',
                )

                nu = seek_min.x[0]

                hash_meteor = f"{data['CODE'][i]}"
                # ADD PARTICLE
                sim.add(
                    m = 0.0e0,
                    a = sma,
                    e = ecc,
                    inc= incl,
                    Omega = Ome,
                    omega = ome,
                    f = nu,
                    hash = hash_meteor
                )
                index.append(hash_meteor)
                with open(file_coll,'a') as f:
                    f.write(f'Body : {data['CODE'][i]}, hash {sim.particles[hash_meteor].hash}\n')
                if active_cl:
                    for j in range(n_clones):
                            #sma = (data['AXIS'][i] + \
                            #    (2.*np.random.rand() - 1)*data['AXISER'][i] )* au.value
                            ecc = data['ECC'][i] + \
                                (2.*np.random.rand() - 1)*data['ECCERR'][i]
                            incl = np.radians(data['INCL'][i] + \
                                (2.*np.random.rand() - 1)*data['INCLER'][i])
                            Ome = np.radians(data['NODE'][i] + \
                                (2.*np.random.rand() - 1)*data['NODEER'][i])
                            ome = np.radians(data['ARGUP'][i] + \
                                (2.*np.random.rand() - 1)*data['ARGUPER'][i])

                            seek_min = minimize(
                                func_dist,
                                x0=np.radians(0.0e0),
                                args=
                                    (
                                    sma,
                                    ecc,
                                    incl,
                                    Ome,
                                    ome,
                                    EP,
                                    ),
                                bounds=[(-2.*np.pi, 2.*np.pi)],
                                method='Nelder-Mead',
                            )

                            nu = seek_min.x[0]
                            hash_cl = f"{data['CODE'][i]}_cl_{j}"
                            sim.add(
                                m = 0.0e0,
                                a = sma,
                                e = ecc,
                                inc= incl,
                                Omega = Ome,
                                omega = ome,
                                f = nu,
                                hash = hash_cl
                            )
                            index.append(hash_cl)
                            with open(file_coll,'a') as f:
                                f.write(f'Body : {hash_cl}, hash {sim.particles[hash_cl].hash}\n')
                                f.write(f'a = {sma}, e = {ecc}, i = {incl}, O = {Ome}, o = {ome}, f = {nu}\n')

    time_array = np.arange(0,-1.*sim_time*year,-1.*sim_step*year)
    count_file = int(0)


print("Total number of bodies:",len(index))
ps = sim.particles



nbd = len(index)
flush = int(len(time_array)/100)

x = np.zeros((nbd, flush))
y = np.zeros((nbd, flush))
z = np.zeros((nbd, flush))
vx = np.zeros((nbd, flush))
vy = np.zeros((nbd, flush))
vz = np.zeros((nbd, flush))
a = np.zeros((nbd, flush))
e = np.zeros((nbd, flush))
I = np.zeros((nbd, flush))
O = np.zeros((nbd, flush))
o = np.zeros((nbd, flush))
f = np.zeros((nbd, flush))
M = np.zeros((nbd, flush))


count = int(0)
pbar = tqdm(total=len(time_array))
for ti, t in enumerate(time_array):
    sim.integrate(t)
    pbar.update(1)


    for j, idx in enumerate(index):
        p = ps[f"{idx}"]
        x[j][count] = p.x
        y[j][count] = p.y
        z[j][count] = p.z
        vx[j][count] = p.vx
        vy[j][count] = p.vy
        vz[j][count] = p.vz
        if j > 0:
            a[j][count] = p.a
            e[j][count] = p.e
            I[j][count] = p.inc
            O[j][count] = p.Omega
            o[j][count] = p.omega
            f[j][count] = p.f
            M[j][count] = p.M

    count =  int(count+1)
    if count == flush:
        count = int(0)

        sim.save_to_file(str(checkpoint))
        file_name = f"{save_file}_{count_file}.h5"
        with h5py.File(file_name, "w") as hf:
            hf.create_dataset("index", data=index)
            hf.create_dataset("time", data=time_array)
            hf.create_dataset("x", data=x)
            hf.create_dataset("y", data=y)
            hf.create_dataset("z", data=z)
            hf.create_dataset("vx", data=vx)
            hf.create_dataset("vy", data=vy)
            hf.create_dataset("vz", data=vz)
            hf.create_dataset("a", data=a)
            hf.create_dataset("e", data=e)
            hf.create_dataset("i", data=I)
            hf.create_dataset("Omega", data=O)
            hf.create_dataset("omega", data=o)
            hf.create_dataset("f", data=f)
            hf.create_dataset("M", data=M)
            hf.create_dataset("met_code", data=data['CODE'])
            hf.create_dataset("met_group", data=data["GROUP"])

        count_file = int(count_file+1)




pbar.close()


