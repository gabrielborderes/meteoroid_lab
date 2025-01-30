import numpy as np
import rebound
import socket
import sys
import h5py
from tqdm import tqdm
import scipy.constants as constants

from dasst.ejection_models.comets import sublimation

from config_3200Phaethon import (
    data_folder,
    init_sim_file,
    solar_system_objects,
    outgass_time,
    body_radius,
    albedo,
    solar_luminosity,
    sublimation_heat,
    body_mass,
    particle_bulk_density,
    gas_molecule_mass,
    surface_temperature_coeff,
    K_drag,
    N_part,
    min_size_log,
    max_size_log,
    particle_file,
)

np.random.seed(426347637)


#  SETING PARTICLE SPHERICAL SHELL
#####################
# TEST PARTILCES


if not init_sim_file.is_file():
    sim = rebound.Simulation()
    # Get data from NASA Horizons
    try:
        sim.units = ("km", "s", "kg")
        date = "2000-01-01 00:00"
        for i in range(len(solar_system_objects)):
            sim.add(solar_system_objects[i], date=date, hash=solar_system_objects[i])
    except socket.error:
        print("A socket error occured. Maybe Horizons is down?")
        sys.exit(0)  # we ignore the error and exit


    sim.N_active = sim.N - 1
    t = -1850.5*365.25*24.*3600.
    sim.integrate(t)
    sim.save_to_file(str(init_sim_file))

sim = rebound.Simulation(str(init_sim_file))

n_bd = len(solar_system_objects)
ps = sim.particles

ejection_possible = np.full((len(outgass_time), ), False, dtype=bool)

masslosses = np.empty_like(outgass_time)
helio_distances = np.empty_like(outgass_time)
for ti, t in tqdm(enumerate(outgass_time), total=len(outgass_time)):
    sim.integrate(t)

    psc = ps[solar_system_objects[-1]]
    pss = ps["Sun"]
    helio_distance = np.sqrt(
        (psc.x - pss.x) ** 2 + (psc.y - pss.y) ** 2 + (psc.z - pss.z) ** 2
    )
    helio_distances[ti] = helio_distance * 1.0e3  # in (m)
    masslosses[ti] = sublimation.whipple_1951.dMdt(
        helio_distances[ti],
        body_radius,
        albedo,
        solar_luminosity,
        sublimation_heat,
    )

    min_size = 10.**min_size_log
    min_mass = (
        particle_bulk_density * (4.0 / 3.0) * np.pi * min_size ** 3
    )
    max_vel = sublimation.whipple_1951.velocity(
        helio_distances[ti],
        body_radius,
        body_mass,
        particle_bulk_density,
        min_mass,
        masslosses[ti],
        gas_molecule_mass,
        surface_temperature_coeff,
        K_drag,
    )
    ejection_possible[ti] = np.logical_and(np.imag(max_vel) == 0, np.logical_not(np.isnan(max_vel)))
    if ejection_possible[ti]:
        print(f"max vel={max_vel} m/s @ {helio_distances[ti]/constants.au} AU")

trunc_masslosses = masslosses.copy()
trunc_masslosses[np.logical_not(ejection_possible)] = 0
N_activity = (trunc_masslosses/np.sum(trunc_masslosses))*N_part
N_activity = np.round(N_activity).astype(np.int64)

N_part = np.sum(N_activity)

shell_part = np.zeros((N_part, 4))
for i in range(N_part):
    x = (max_size_log - min_size_log)*np.random.rand() + min_size_log
    shell_part[i, 0] = 10.**x
    theta = 2.*np.pi*np.random.rand()
    phi = np.arccos(1.-2.*np.random.rand())
    shell_part[i, 1] = np.cos(theta)*np.sin(phi)
    shell_part[i, 2] = np.sin(theta)*np.sin(phi)
    shell_part[i, 3] = np.cos(phi)


particle_velocities = np.full((N_part, ), np.nan, dtype=np.float64)
particle_distances = np.full((N_part, ), np.nan, dtype=np.float64)

part_ind = 0
max_tries = 100
for step_i, num in tqdm(enumerate(N_activity), total=len(outgass_time)):
    for pi in range(num):
        particle_mass = (
            particle_bulk_density * (4.0 / 3.0) * np.pi * shell_part[part_ind, 0] ** 3
        )
        particle_distances[part_ind] = helio_distances[step_i]
        particle_velocities[part_ind] = sublimation.whipple_1951.velocity(
            helio_distances[step_i],
            body_radius,
            body_mass,
            particle_bulk_density,
            particle_mass,
            masslosses[step_i],
            gas_molecule_mass,
            surface_temperature_coeff,
            K_drag,
        )
        tries = 0
        while not np.logical_and(np.imag(particle_velocities[part_ind]) == 0, np.logical_not(np.isnan(particle_velocities[part_ind]))):
            x = (max_size_log - min_size_log)*np.random.rand() + min_size_log
            shell_part[part_ind, 0] = 10.**x
            particle_mass = (
                particle_bulk_density * (4.0 / 3.0) * np.pi * shell_part[part_ind, 0] ** 3
            )
            particle_velocities[part_ind] = sublimation.whipple_1951.velocity(
                helio_distances[step_i],
                body_radius,
                body_mass,
                particle_bulk_density,
                particle_mass,
                masslosses[step_i],
                gas_molecule_mass,
                surface_temperature_coeff,
                K_drag,
            )
            tries += 1
            if tries > max_tries:
                particle_velocities[part_ind] = np.nan
                break
        part_ind += 1

print(f"{N_part} particles generated")
print(f"{np.sum(np.isnan(particle_velocities))} particles failed ejecting")



with h5py.File(str(particle_file), "w") as hf:
    hf.create_dataset("times", data=outgass_time)
    hf.create_dataset("shell_part", data=shell_part)
    hf.create_dataset("radius", data=shell_part[:, 0])
    hf.create_dataset("N_activity", data=N_activity)
    hf.create_dataset("particle_velocities", data=particle_velocities)
    hf.create_dataset("particle_distances", data=particle_distances)
    hf.create_dataset("masslosses", data=masslosses)
    hf.create_dataset("helio_distances", data=helio_distances)
    hf.create_dataset("ejection_possible", data=ejection_possible)





