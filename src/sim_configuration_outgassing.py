import pathlib
import numpy as np

data_folder = pathlib.Path("./data").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

init_sim_file = data_folder / "cache_12p_init.bin"

# SIMULATION Parameters
sim_years = 750
Frame_betw_days = 4
day = 24.0 * 3600.0
Frame_betw = Frame_betw_days * day
Delta_T = sim_years * 365.25 * day
c2 = 629*day #Outgassing
c1 = c2 - day
a = np.arange(0, c1, day)
b = np.arange(c2, Delta_T, Frame_betw)
times = np.concatenate([a,b])
Nout = len(times)

# WHYPLE Model Parameter ####
body_radius = 34e3  # (m)
albedo = 0.2
solar_luminosity = 4e26
sublimation_heat = 1.88e6
body_mass = 2.0579526276115536e+16
particle_bulk_density = 3.5e3
gas_molecule_mass = 20 * 1.661e-24 * 1e-3
surface_temperature_coeff = 300
K_drag = 26.0 / 9.0

# THE COMET MUST BE THE LAST IN THE LIST
solar_system_objects = [
    "Sun",
    "Mercury",
    "Venus",
    "Earth",
    "LUNA",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "PONS-BROOKS",
]
