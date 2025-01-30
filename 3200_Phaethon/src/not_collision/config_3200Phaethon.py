import pathlib
import numpy as np

data_folder = pathlib.Path("data").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

init_sim_file = data_folder / "cache_3200Phaethon_init.bin"
particle_file = data_folder / "particles_3200_Phaethon.h5"
ephem_file = data_folder / f"ephemerides_3200_Phaethon_all.h5"



# SIMULATION Parameters
day = 24.0 * 3600.0
year = 365.25 * day

#sim_tf = 2100. * year
sim_tf = 21. * year
outgass_tf = 1.35 * year

sim_delta = 2. * year
outgass_delta = 1. * day

outgass_time = np.arange(0,outgass_tf ,outgass_delta)
sim_wo_outgass = np.arange(outgass_tf + day, sim_tf ,sim_delta)
Nog = len(outgass_time)

sim_time = np.concatenate([outgass_time,sim_wo_outgass])
Nout = len(sim_time)


# WHYPLE Model Parameter ####
body_radius = 3.125e3  # (m)
albedo = 0.1066
solar_luminosity = 4e26
sublimation_heat = 1.88e6
body_mass = 1.4e+14
particle_bulk_density = 1.5e3
gas_molecule_mass = 20 * 1.661e-24 * 1e-3
surface_temperature_coeff = 300
K_drag = 26.0 / 9.0

# Particle parameters
N_part = 10
min_size_log = -4
max_size_log = -1

# Comet physical parameters (dist-km)
Eje_p_R = 10.0 * body_radius
T_comet = 3.604  # hour (https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=3200)
RotVel = 2.0 * np.pi / (T_comet * 3600.0)



#THE COMET MUST BE THE LAST IN THE LIST
solar_system_objects = [
    "Sun",
    "Mercury",
    "Venus",
    "399",
    "301",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "1983 TB",
]
