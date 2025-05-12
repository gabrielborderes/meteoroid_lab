import pathlib
import numpy as np

data_folder = pathlib.Path("data").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

dictionary_file = data_folder / "dictionary.txt"
init_sim_file = data_folder / "cache_21P_init.bin"
particle_file = data_folder / "particles_21P.h5"
ephem_file = data_folder / f"ephemerides_21P_all.h5"
checkpoint = data_folder / "reboot.bin"


factor_flush = 100


srp_active = True

# SIMULATION Parameters
day = 24.0 * 3600.0
year = 365.25 * day

#sim_tf = 10000. * year #long
sim_tf = 168. * year
#sim_tf = 25. * year
outgass_tf = 7. * year

sim_delta = 1 * year #long
#sim_delta = 1 * day
outgass_delta = 1 * day

outgass_time = np.arange(0,outgass_tf ,outgass_delta)
sim_wo_outgass = np.arange(outgass_tf + day, sim_tf ,sim_delta)
Nog = len(outgass_time)

sim_time = np.concatenate([outgass_time,sim_wo_outgass])
Nout = len(sim_time)


# WHYPLE Model Parameter ####
body_radius = 1.0e3  # (m)
albedo = 0.1066
solar_luminosity = 4e26
sublimation_heat = 1.88e6
body_mass = 2.51e12
particle_bulk_density = 1.5e3 # Using Halley densiity
gas_molecule_mass = 20 * 1.661e-24 * 1e-3
surface_temperature_coeff = 300
K_drag = 26.0 / 9.0

# Particle parameters
N_part = 10000

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
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "1900 Y1",
]
