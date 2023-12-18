'''
Python interface to the LOST Open-source Star Tracker.

Depends on the LOST binary being properly built & in the same folder as this
Python module.

# Usage
-------

For general usage, first run :func:`lost.database` to initialize the star
identification database. Then, run :func:`lost.identify` to identify
images, or :func:`lost.generate` to generate simulated images. Example:

```
import lost
lost.database()

# read in the test image as a numpy array
import imageio.v3 as imageio
im = imageio.imread('img_7660.png')

# identify attitude of satellite using default parameters
result = lost.identify(im)

# pretty print attitude info with JSON module
import json
print(json.dumps(result, indent=True))

# generate image using default parameters
raw, annotated = lost.generate()
show(raw)
show(annotated)
```

# Detailed documentation
------------------------

For detailed usage documentation, run the LOST binary's help commands:

* `./lost database --help`: database command
* `./lost pipeline --help`: image generation & identification pipeline command
'''


import imageio
import numpy as np
import subprocess
import os
import pathlib


# TODO:
# [x] Bare bones MVP of database, generation, and identification working
# [x] Organize this file
# [x] Solidify interfaces (especially around setting parameters)
# [x] Properly bundle into a python package (wheel, .whl)
# [x] Bundle package to minimal set of files for install
# [x] Review package files and consider adding symlink for lost binary
# [x] Change from `X_default_args` to `X_args(overrides: dict)`
# [ ] ** Example code/usage, incl. in docstrings
# [ ] Ensure all parameters are available
# [ ] Test that all parameters work correctly
# [ ] Thorough docstrings
# [ ] Properly handle CLI print output
# [ ] Error trapping/self-consistency checking
# [ ] Type hints
# [ ] Add LOST help functionality
# [ ] Update readme (or make python-specific readme?)
# [ ] Filesystem pipes instead of files
# [ ] Plotting helper commands using matplotlib
# [ ] Propose future work/splitting things out (python vs cli vs others)
# [ ] Incorporate wheel check https://github.com/jwodder/check-wheel-contents
# [ ] Architecture packages: https://github.com/python-poetry/poetry/issues/5205
# [ ] Investigate "Package would be ignored" warning on wheel build ('lost.tmp')

# paths for important things
LOST_DIR_PATH = os.path.dirname(os.path.realpath(__file__))
LOST_EXECUTABLE_PATH = f"{LOST_DIR_PATH}/lost"
TEMP_DIR_PATH = f'{LOST_DIR_PATH}/tmp'


# paths for temporary files
RAW_INPUT_PATH = f'{TEMP_DIR_PATH}/raw-input.png'
ANNOTATED_INPUT_PATH = f'{TEMP_DIR_PATH}/input.png'
ANNOTATED_OUTPUT_PATH = f'{TEMP_DIR_PATH}/annotated_output.png'
ATTITUDE_PATH = f'{TEMP_DIR_PATH}/attitude.txt'
DATABASE_PATH = f'{TEMP_DIR_PATH}/tmp_database.dat'


# make temporary directory if it doesn't exist
pathlib.Path(TEMP_DIR_PATH).mkdir(exist_ok=True)


# TODO: revisit return type (string output from LOST?)
def lost(args: dict) -> None:
    '''Call LOST, passing a `dict` of arguments.'''
    lost_cli_list(flatten_dict_to_list(args))


def lost_cli_list(args: list) -> None:
    '''Call LOST command line interface, passing a `list` of arguments.'''
    pass_args = [str(arg) for arg in args]
    subprocess.run([LOST_EXECUTABLE_PATH] + pass_args, cwd=LOST_DIR_PATH)


#######################
# DATABASE GENERATION #
#######################


def database_args(overrides: dict = {}) -> dict:
    '''
    Returns dictionary of default arguments for :func:`lost.database`.
    
    Applies `overrides` dict over generated/default values. For example,
    `database_args({'--max-stars': 4000})` will result in '--max-stars' mapping
    to 4000 in the returned dict.
    '''
    args = {
        'database': None,
        '--max-stars': 5000,
        '--kvector': None,
        '--kvector-min-distance': 0.2,
        '--kvector-max-distance': 15.0,
        '--kvector-distance-bins': 10_000,
        '--output': DATABASE_PATH,
    }
    args.update(overrides)
    return args


def database(args: dict = database_args()) -> None:
    '''
    Calls LOST's database generation command.

    Must be called before :func:`lost.identify` to initialize LOST.

    See :func:`lost.database_args` for arguments.
    '''
    # TODO: validate sanity of database parameters
    lost(args)


####################
# IMAGE GENERATION #
####################


def generate_args(overrides: dict = {},
                  generate_raw: bool = True,
                  generate_annotated: bool = True) -> dict:
    '''
    Returns `dict` of default arguments for :func:`lost.generate`.

    Applies `overrides` dict over generated/default values. For example,
    `generate_args({'--generate-de': 8})` will result in '--generate-de' mapping
    to 8 in the returned dict.

    If `generate_raw` is `True`, include command to generate raw input image.

    If `generate_annotated` is `True`, include command to generate annotated
        input image.
    '''
    args = {
        'pipeline': None,
        '--generate': '1',
        '--generate-x-resolution': 1024,
        '--generate-y-resolution': 1024,
        '--fov': 30,
        # what to do with read_noise_stddev=0.05 ?
        # this is an unrecognized parameter?
        # '--generate-reference-brightness', str(reference_brightness),
        '--generate-spread-stddev': 1,
        '--generate-ra': 88,
        '--generate-de': 7,
        '--generate-roll': 0,
    }

    # add image generation arguments as applicable
    if generate_raw:
        args['--plot-raw-input'] = RAW_INPUT_PATH
    if generate_annotated:
        args['--plot-input'] = ANNOTATED_INPUT_PATH

    args.update(overrides)
    return args


def generate(args: dict = generate_args()) -> \
        tuple[np.ndarray, np.ndarray]:
    '''
    Calls LOST's image generation command, returning generated images.

    See :func:`lost.generate_args` for parameters.

    Returns `(raw_result: np.ndarray, annotated_result: np.ndarray)`.
    '''
    lost(args)

    raw_result = None
    if '--plot-raw-input' in args:
        # TODO: use appropriate raw input path
        raw_result = imread(RAW_INPUT_PATH)
        delete_file(RAW_INPUT_PATH)

    annotated_result = None
    if '--plot-input' in args:
        # TODO: use appropriate annotated input path
        annotated_result = imread(ANNOTATED_INPUT_PATH)
        delete_file(ANNOTATED_INPUT_PATH)

    return (raw_result, annotated_result)


########################
# IMAGE IDENTIFICATION #
########################


def identify_args(overrides: dict = {}) -> dict:
    '''
    Returns `dict` of default arguments for :func:`lost.identify`.
    
    Applies `overrides` dict over generated/default values. For example,
    `identify_args({'--star-id-algo': 'tetra'})` will result in
    '--star-id-algo' mapping to 'tetra' in the returned dict.
    '''
    args = {
        'pipeline': None,
        '--png': RAW_INPUT_PATH,
        '--focal-length': 49,
        '--pixel-size': 22.2,
        '--centroid-algo': 'cog',  # 'cog', 'dummy', 'iwcog'
        '--centroid-mag-filter': 5,
        '--database': DATABASE_PATH,
        '--star-id-algo': 'py',  # 'dummy', 'gv', 'py', 'tetra'
        '--angular-tolerance': 0.05,
        '--false-stars': 1000,
        '--max-mismatch-prob': 0.0001,
        '--attitude-algo': 'dqm',  # 'dqm' (Davenport Q), 'triad', 'quest'
        '--print-attitude': ATTITUDE_PATH,
    }
    args.update(overrides)
    return args


def identify(image: np.ndarray, args: dict = identify_args()) -> dict:
    '''
    Identifies `image: np.ndarray`, returning attitude information as `dict`.

    Running :func:`lost.database` is a prerequisite.

    See :func:`lost.identify_args` for parameters.

    Returns dictionary of attitude information:
    ```
    {
        "attitude_known": int,   # 1 if identified successfully
        "attitude_ra": float,    # right ascension, degrees
        "attitude_de": float,    # declination, degrees
        "attitude_roll": float,  # roll, degrees
        "attitude_i": float,     # attitude quaternion i
        "attitude_j": float,     # attitude quaternion j
        "attitude_k": float,     # attitude quaternion k
        "attitude_real": float,  # attitude quaternion real part
    }
    ```
    '''
    # save the given image to disk so LOST can use it
    imwrite(RAW_INPUT_PATH, image)

    # identify image
    lost(args)

    # parse/load attitude file
    with open(ATTITUDE_PATH) as f:
        read_data = f.read()

    result = {}
    rows = read_data.split('\n')
    for row in rows:
        if row == '':
            continue
        sp = row.split(' ')
        if sp[0] == 'attitude_known':
            result[sp[0]] = int(sp[1])
        else:
            result[sp[0]] = float(sp[1])

    # clean up leftover junk (image, attitude file)
    delete_file(ATTITUDE_PATH)
    delete_file(RAW_INPUT_PATH)

    # return attitue information, optionally output
    return result


###########################
# MISCELLANEOUS UTILITIES #
###########################


def flatten_dict_to_list(dictionary: dict) -> list:
    '''
    'flattens' a dictionary into a list, skipping None values.

    {'a': 'b', 'c': None, 'd': 3.14} -> ['a', 'b', 'c', 'd', 3.14]
    '''
    arr = []
    for key, value in dictionary.items():
        arr.append(key)
        if not (value is None):
            arr.append(value)
    return arr


def imread(path: str) -> np.ndarray:
    return imageio.imread(path)


def imwrite(path: str, image: np.ndarray) -> None:
    imageio.imwrite(path, image)


def delete_file(path: str) -> None:
    pathlib.Path(path).unlink()
