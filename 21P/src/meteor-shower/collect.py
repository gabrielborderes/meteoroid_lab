import h5py
import numpy as np
from tqdm import tqdm

from config import (
    data_folder,
    solar_system_objects,
    particle_file,
    factor_flush,
)

with h5py.File(particle_file, 'r') as f:
    shell_part = f["shell_part"][:]

N_part = shell_part.shape[0]

nb = len(solar_system_objects)

for fl in range(factor_flush+1):

    h5_files = list(data_folder.glob(f"part_{fl}/ephemerides_21P_*.h5"))
    ephem_file = data_folder / f"part_{fl}/ephemerides_21P_all.h5"



    data_sets = {}

    state_keys = [
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "a",
        "e",
        "i",
        "O",
        "o",
        "f",
        "M",
        "Per",
    ]


    t = None

    for file in tqdm(h5_files, total=len(h5_files)):
        with h5py.File(file, "r") as hf:
            index = hf["index"][()]
            if t is None:
                t = hf["t"][()]

            if "index" not in data_sets:
                data_sets["index"] = [index]
            else:
                data_sets["index"].append(index)
            # print(index.size - hf["x"].shape[0])
            for key in state_keys:
                sim_data = hf[key][()]

                pdata = sim_data[nb:, :]
                if key not in data_sets:
                    data_sets[key] = [pdata]
                else:
                    data_sets[key].append(pdata)

                bkey = "body_" + key
                bdata = sim_data[:nb, :]
                if bkey not in data_sets:
                    data_sets[bkey] = bdata

    for key in state_keys:
        data_sets[key] = np.concatenate(data_sets[key], axis=0)
        # print(data_sets[key].shape)
    data_sets["index"] = np.concatenate(data_sets["index"], axis=0)

    sort_ind = np.argsort(data_sets["index"])
    data_sets["index"] = data_sets["index"][sort_ind]
    for key in state_keys:
        data_sets[key] = data_sets[key][sort_ind]

    with h5py.File(str(ephem_file), "w") as hf:
        for key in data_sets:
            hf.create_dataset(key, data=data_sets[key])
        hf.create_dataset("t", data=t)


    for file in h5_files:
        file.unlink()
