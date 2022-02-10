# %%
import os
import random
import time
import numpy as np
from p_tqdm import p_umap
# %%
def random_distribution_over_sphere(lower, upper):
    """
    Generate a random distribution over a sphere.
    """
    # Generate a random point on the surface of a sphere.
    v = random.uniform(-1, 1)
    
    return np.arcsin(v)

# %%
PATH_TO_LOST = '/home/sathvikc/lost'
DATA_PATH = '/home/sathvikc/lost/ML/data/'
print(os.path.exists(DATA_PATH))

if not os.path.isdir(DATA_PATH):
    os.system(f'mkdir {DATA_PATH}')
# %%
os.system(f'cd {PATH_TO_LOST} && yes | make clean && CXXFLAGS=-O3 make -j{len(os.sched_getaffinity(0))}')

# %% 
def generateImage(num):
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
    BORESIGHT_RIGHT_ASCENSION = random.randint(0, 360)
    BORESIGHT_DECLINATION = random_distribution_over_sphere(-180, 180)
    BORESIGHT_ROLL = random.randint(0, 360)
    MOTION_BLUR_DIRECTION_RIGHT_ASCENSION = random.randint(0, 1)
    MOTION_BLUR_DIRECTION_DECLINATION = random.randint(0, 1)
    MOTION_BLUR_DIRECTION_ROLL = random.randint(0, 1)
    FILE_NAME = str(num) + str(time.time()).replace('.', '')

    # Run the command to generate the data.
    os.system(f'cd {PATH_TO_LOST} '+ 
    '&& ./lost pipeline generate 1 ' + 
    f'{RESOLUTION} {RESOLUTION} {PIXEL_SIZE} {FOCAL_LENGTH} {OBSERVED_REFERENCE_STAR_BRIGHTNESS} '+
    f'{STAR_SPREAD_STDDEV} {CAMERA_SENSITIVITY} {DARK_CURRENT} {NOISE_STDDEV} {EXPOSURE_TIME} '+ 
    f'{READOUT_TIME} {ENABLE_SHOT_NOISE} {OVERSAMPLING} {BORESIGHT_RIGHT_ASCENSION} '+ 
    f'{BORESIGHT_DECLINATION} {BORESIGHT_ROLL} {MOTION_BLUR_DIRECTION_RIGHT_ASCENSION} '+ 
    f'{MOTION_BLUR_DIRECTION_DECLINATION} {MOTION_BLUR_DIRECTION_ROLL} done plot_raw_input '+ 
    f'{FILE_NAME}.png done >/dev/null 2>&1')

    # Move the data to the data folder.
    os.system(f'cd {PATH_TO_LOST} && mv {FILE_NAME}.png {DATA_PATH} >/dev/null 2>&1')

    f = open('data.txt', 'a')
    f.write(f'{FILE_NAME}.png, ({BORESIGHT_RIGHT_ASCENSION},{BORESIGHT_DECLINATION},{BORESIGHT_ROLL})\n')
    f.close()



# %%
def main():
    num_images = int(input("How many images do you want to generate? "))
    f = open('data.txt', 'w')
    generated = p_umap(generateImage, range(num_images))

  # %%  
if __name__ == "__main__":
    main()