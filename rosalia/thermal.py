import os
import numpy as np
from scipy.ndimage import zoom

class ThermalBackground:
    thermal = {"F062": [0.011, 0.011, 0.011, 0.012, 0.012, 0.012, 0.012, 0.012, 0.011, 0.012, 0.011, 0.012, 0.012, 0.012, 0.011, 0.012, 0.012, 0.012],
               "F087": [0.012, 0.012, 0.011,0.012, 0.012, 0.012, 0.013, 0.012, 0.012, 0.011, 0.012, 0.013, 0.011, 0.012, 0.011, 0.012, 0.012, 0.013], 
               "F106": [0.014, 0.014, 0.015, 0.015, 0.014, 0.014, 0.015, 0.013, 0.013, 0.013, 0.014, 0.014, 0.012, 0.014, 0.012, 0.015, 0.015, 0.015],
               "F129": [0.010,	0.010, 0.010, 0.010, 0.010,	0.011, 0.010,0.010,	0.010, 0.010, 0.010, 0.012, 0.010, 0.011, 0.009, 0.010, 0.010, 0.010],
               "F146": [0.15, 0.15, 0.16, 0.15, 0.16, 0.17, 0.15, 0.16, 0.16, 0.15, 0.15, 0.17, 0.15, 0.15, 0.15, 0.15, 0.15, 0.16],
               "F158": [0.016, 0.017, 0.016, 0.017, 0.017, 0.018, 0.016, 0.016, 0.017, 0.016, 0.017, 0.019, 0.016, 0.017, 0.016, 0.016, 0.016, 0.017],
               "F184": [0.112, 0.112, 0.106, 0.118, 0.114, 0.114, 0.122, 0.113, 0.113, 0.116, 0.113, 0.115, 0.118, 0.118, 0.110, 0.121, 0.121, 0.123],
               "F213": [2.57, 2.46, 2.36, 2.66, 2.51, 2.59, 2.67, 2.52, 2.56, 2.65, 2.53, 2.65, 2.69, 2.66, 2.47, 2.74, 2.79, 2.80],
               "grism":[0.076, 0.069, 0.068, 0.099, 0.075, 0.078, 0.101, 0.092, 0.080, 0.075, 0.070, 0.079, 0.096, 0.079, 0.070, 0.097, 0.098, 0.094],
               "prism":[0.098, 0.101, 0.072, 0.132, 0.138, 0.100, 0.160, 0.181, 0.151, 0.099, 0.110, 0.089, 0.125, 0.128, 0.078, 0.141, 0.171, 0.139]}


def read_sca_tables(filename):
    sca_data = {}
    current_key = None
    current_rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Detect SCA headers (e.g., "SCA1", "SCA 8", "SCA 10")
            if line.startswith("SCA"):
                # Save previous block if present
                if current_key is not None and current_rows:
                    sca_data[current_key] = np.array(current_rows, dtype=float)
                    current_rows = []

                # Normalize SCA name (remove spaces)
                
                current_key = line.replace(" ", "")
                SCAi = str(int(current_key.replace("SCA",""))).zfill(2)
                current_key = "SCA"+SCAi
                continue

            # Otherwise it's a row of numbers
            numbers = line.split()
            current_rows.append([float(x) for x in numbers])

        # Store last table after loop ends
        # Flip the y-axis. origin=bottom.
        if current_key is not None and current_rows:
            # print(np.array(current_rows, dtype=float).shape)
            sca_data[current_key] = np.array(current_rows, dtype=float) # np.flip(np.array(current_rows, dtype=float), axis=0)
            # sca_data[current_key] = np.flip(np.array(current_rows, dtype=float), axis=0)

        for i in range(18):
            key = "SCA" + str(i+1).zfill(2)
            sca_data[key] = np.flip(sca_data[key], axis=0)
    return sca_data


def get_thermal_background(bandpass, SCA):
    # Example usage:

    if bandpass is "F146" or bandpass is "F213" or bandpass is "F184":
        filename = os.environ["ROSALIACACHE"] + "CORE/THERMAL/"+bandpass+"_thermal.txt"   # your text file
        sca_tables = read_sca_tables(filename)
        target_size = 4088
        SCAi = str(SCA).zfill(2)
        SCAkey = "SCA" + SCAi
        zoom_factor = target_size / sca_tables[SCAkey].shape[0]

        upscaled = zoom(sca_tables[SCAkey], zoom_factor, order=3)  # cubic interpolation
        return(upscaled)   # (4088, 4088)
    
    else:
        thermal_per_sca = ThermalBackground.thermal[bandpass]
        thermal_background = np.zeros((4088,4088)) + thermal_per_sca[SCA-1]
        return(thermal_background)
        



