# %%
import os
import random
import time
import numpy as np
from p_tqdm import p_umap
import csv
from constants import *
# %%
def random_distribution_over_sphere(lower, upper):
    """
    Generate a random distribution over a sphere.
    """
    # Generate a random point on the surface of a sphere.
    v = random.uniform(-1, 1)
    
    return np.arcsin(v)

# %% 

def getXYArray(num):
    x_array = []
    y_array = []
    with open(f'{PATH_TO_LOST}/{num}.txt') as file:
        for line in file: 
            if line.find('expected_id') == -1:
                value = float(line.split(' ')[1])
                if line.find('x') != -1:
                    x_array.append(value)
                elif line.find('y') != -1: 
                    y_array.append(value)
    return x_array, y_array 

# %%
def generateCentroid(num):
    #./lost pipeline generate 1 1000 1000 0 95 50 0.1 0.001 0 0 0 0 0 4 0 0 0 
    random.seed(time.time())
    RESOLUTION = 1000
    PIXEL_SIZE = 0
    FOCAL_LENGTH = random.randint(95, 100)
    OBSERVED_REFERENCE_STAR_BRIGHTNESS = random.randint(50, 400)
    STAR_SPREAD_STDDEV = random.uniform(0.1, 2.0)
    CAMERA_SENSITIVITY = random.uniform(0.001, 0.05)
    DARK_CURRENT = random.uniform(0.0, 0.1)
    NOISE_STDDEV = random.uniform(0.0, 0.3)
    EXPOSURE_TIME = random.uniform(0.0, 0.5)
    READOUT_TIME = 0 
    ENABLE_SHOT_NOISE = random.randint(0, 1)
    OVERSAMPLING = 4
    BORESIGHT_RIGHT_ASCENSION = random.randint(0, 359)
    BORESIGHT_DECLINATION = random_distribution_over_sphere(-180, 180)
    BORESIGHT_ROLL = random.randint(0, 359)
    MOTION_BLUR_DIRECTION_RIGHT_ASCENSION = random.uniform(0.001, 1)
    MOTION_BLUR_DIRECTION_DECLINATION = random.uniform(0.001, 1)
    MOTION_BLUR_DIRECTION_ROLL = random.uniform(0.001, 5)

    CENTROID_ALGO = 'cog'

    # Run the command to generate the data.
    os.system(f'cd {PATH_TO_LOST} '+ 
    '&& ./lost pipeline generate 1 ' + 
    f'{RESOLUTION} {RESOLUTION} {PIXEL_SIZE} {FOCAL_LENGTH} {OBSERVED_REFERENCE_STAR_BRIGHTNESS} '+
    f'{STAR_SPREAD_STDDEV} {CAMERA_SENSITIVITY} {DARK_CURRENT} {NOISE_STDDEV} {EXPOSURE_TIME} '+ 
    f'{READOUT_TIME} {ENABLE_SHOT_NOISE} {OVERSAMPLING} {BORESIGHT_RIGHT_ASCENSION} '+ 
    f'{BORESIGHT_DECLINATION} {BORESIGHT_ROLL} '+
    f'centroid {CENTROID_ALGO} done print_centroids'+ 
    f'{num}.txt done')

    # Move the data to the data folder.

    x_vals, y_vals = getXYArray(num)

    os.system(f'cd {PATH_TO_LOST} && rm -rf {num}.txt {DATA_PATH} >/dev/null 2>&1')    

    with open('data.csv', 'a') as f:
        csv_writer = csv.writer(f, delimiter = ',')
        csv_writer.writerow([x_vals, y_vals, BORESIGHT_RIGHT_ASCENSION, BORESIGHT_DECLINATION, BORESIGHT_ROLL])

# %%
def main():
    num_images = int(input("How many images do you want to generate? "))
    f = open('data.csv', 'w')
    with open('data.csv', 'a') as f:
        f.write('x_vals,y_vals,ra,de,roll\n')
    generateCentroid(0)
    # generated = p_umap(generateCentroid, range(num_images))

  # %%  
if __name__ == "__main__":
    os.system(f'cd {PATH_TO_LOST} && yes | make clean && CXXFLAGS=-O3 make -j{len(os.sched_getaffinity(0))}')
    if not os.path.isdir(DATA_PATH):
        os.system(f'mkdir {DATA_PATH}')
    main()
    os.system(f'mv data.csv {DATA_PATH}')