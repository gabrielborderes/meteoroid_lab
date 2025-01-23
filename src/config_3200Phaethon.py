import pathlib
import numpy as np

data_folder = pathlib.Path("../data").resolve()
if not data_folder.is_dir():
    data_folder.mkdir()

init_sim_file = data_folder / "cache_3200Phaenton_init.bin"

# SIMULATION Parameters
sim_years = 1.35
Frame_betw_days = 1
day = 24.0 * 3600.0
Frame_betw = Frame_betw_days * day
Delta_T = sim_years * 365.25 * day

times = np.arange(0, Delta_T, Frame_betw)
Nout = len(times)

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

# THE COMET MUST BE THE LAST IN THE LIST
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
