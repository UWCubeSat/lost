'''
LOST Open-source Star Tracker (Python Interface)

This module is a Python interface to the LOST Open-source Star Tracker. It is
based on the Command Line Interface (CLI) for LOST.

It depends on the LOST binary being properly built & bundled in the same folder
with this Python module.

To report issues and learn more, visit https://github.com/UWCubeSat/lost.


## Usage
-------

For general usage, first run :func:`lost.database` to initialize the star
identification database.

Then, run :func:`lost.identify` to identify images, or :func:`lost.generate`
to generate simulated images.

In general, `database`, `identify`, and `generate` work like the CLI
equivalents, and arguments are specified in a similar way, too. Methods ending
in `_args` are helper methods that create reasonable arguments for each command,
with the option to override or add parameters:

```Python
args = lost.X_args({ '--foo': 'override value' })
lost.X(args)
```

For detailed LOST documentation, run the LOST binary's help commands:

* `./lost database --help`: database command
* `./lost pipeline --help`: image generation & identification pipeline command

For more usage details, see the following tutorial, docstrings for methods, and
the source code in `__init__.py`.

## Tutorial
-----------

### Setup

Set up LOST by generating the database. You can add overrides to the `args`
`dict` by passing in a `dict` to `database_args`. In this case, the override is
redundant (5000 is the default), but it illustrates the approach:

```Python
import lost

args = lost.database_args({ '--max-stars': 5000 })
lost.database(args)
```

### Load Image

Load our test image (downloadable at https://markasoftware.com/img_7660.png).
It's a PNG loaded from disk, which results in an `np.ndarray` of shape
`(667, 1000, 3)` (height, width, spectra) and data type `uint8`.

```Python
import imageio.v3 as imageio
im = imageio.imread('img_7660.png')
```

### Identify Image

Identify the image. Overrides specified the same way as for `database_args`.

```Python
# identify attitude of satellite
args = lost.identify_args(algo='py')
result = lost.identify(im)

# pretty print attitude info using JSON module
import json
print(json.dumps(result, indent=True))
```

### Generate Images

Generate simulated images. Overrides work as before.

```Python
import matplotlib as mpl
from matplotlib import pyplot as plt

# Show some number of np.ndarray images side-by-side using pyplot.
def show(*ims) -> None:
    mpl.rcParams['figure.dpi'] = 600
    fig, axes = plt.subplots(1, len(ims), tight_layout=True, squeeze=False)
    for i, ax in enumerate(axes.flat):
        ax.imshow(ims[i])
        ax.axis('off')
    plt.show()

# Generate images using LOST.
args = lost.generate_args({ '--generate-de': 8 })
raw1, annotated1 = lost.generate(args)

args = lost.generate_args({ '--generate-de': 5 })
raw2, annotated2 = lost.generate(args)

show(raw1, annotated1, raw2, annotated2)
```
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
# [x] Type hints
# [x] Debug flag to print CLI args
# [x] ** Example code/usage, incl. in docstrings
# [x] Thorough docstrings
# [x] Update readme (or make python-specific readme?)
# [ ] Filesystem pipes instead of files
#     https://tutorialspoint.com/How-to-create-and-use-a-named-pipe-in-Python
# [ ] Error trapping/self-consistency checking (presence of LOST bin, params...)
# [ ] Ensure all parameters are available
# [ ] Test that all parameters work correctly
# [ ] Properly handle CLI print output
# [ ] Add LOST help functionality
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
TETRA_DATABASE_PATH = f'{TEMP_DIR_PATH}/tetra_database.dat'
PY_DATABASE_PATH = f'{TEMP_DIR_PATH}/py_database.dat'  # pyramidal, not python


# module configuration flags
debug_print_cli_args = False


# make temporary directory if it doesn't exist
pathlib.Path(TEMP_DIR_PATH).mkdir(exist_ok=True)


# TODO: revisit return type (string output from LOST?)
def lost(args: dict) -> None:
    '''Call LOST, passing a `dict` of arguments.'''
    lost_cli_list(flatten_dict_to_list(args))


def lost_cli_list(args: list) -> None:
    '''Call LOST command line interface, passing a `list` of arguments.'''
    pass_args = [str(arg) for arg in args]
    if debug_print_cli_args:
        print('calling lost CLI with args:', pass_args)
    subprocess.run([LOST_EXECUTABLE_PATH] + pass_args, cwd=LOST_DIR_PATH)


#######################
# DATABASE GENERATION #
#######################


def database_args(overrides: dict = {}, algo: str = 'py') -> dict:
    '''
    Returns dictionary of default arguments for :func:`lost.database`.

    Applies `overrides` dict over generated/default values. For example,
    `database_args({'--max-stars': 4000})` will result in '--max-stars' mapping
    to 4000 in the returned dict.

    Sets up for pyramidal if `algo` is `'py'`, or tetra if `algo` is `'tetra'`.
    '''
    if algo == 'py':
        args = {
            'database': None,
            '--max-stars': 5000,
            '--kvector': None,
            '--kvector-min-distance': 0.2,
            '--kvector-max-distance': 15.0,
            '--kvector-distance-bins': 10_000,
            '--output': PY_DATABASE_PATH,
        }
    elif algo == 'tetra':
        args = {
            'database': None,
            '--min-mag': 7,
            '--tetra': None,
            '--tetra-max-angle': 12,
            '--output': TETRA_DATABASE_PATH,
        }
    else:
        raise f"Invalid database algo {algo}. Must be 'py' or 'tetra'."
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
                  generate_annotated: bool = True,
                  gen_type: str = 'default') -> dict:
    '''
    Returns `dict` of default arguments for :func:`lost.generate`.

    Applies `overrides` dict over generated/default values. For example,
    `generate_args({'--generate-de': 8})` will result in '--generate-de' mapping
    to 8 in the returned dict.

    If `generate_raw` is `True`, include command to generate raw input image.

    If `generate_annotated` is `True`, include command to generate annotated
        input image.

    If `gen_type` is `'default'`, generates according to example in README. If
        it's `'oresat'`, uses OreSat-like image generation.
    '''
    if gen_type == 'default':
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
    elif gen_type == 'oresat':
        args = {
            'pipeline': None,
            '--generate': 1,
            '--fov': 17,
            '--generate-x-resolution': 1280,
            '--generate-y-resolution': 960,
            '--generate-ra': 79.4232,
            '--generate-de': 46.2072,
            '--generate-roll': 78.2978,
            '--generate-perturb-centroids': 0,
            '--generate-shot-noise': 'true',
            '--generate-read-noise-stddev': 0.01,
            '--generate-dark-current': 0.07,
        }
    else:
        raise f"Invalid gen_type {gen_type}. Must be 'default' or 'oresat'."

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


def identify_args(overrides: dict = {}, algo: str = 'py') -> dict:
    '''
    Returns `dict` of default arguments for :func:`lost.identify`.

    Applies `overrides` dict over generated/default values. For example,
    `identify_args({'--fov': 18})` will result in
    `'--fov'` mapping to `18` in the returned dict.

    Sets up for pyramidal if `algo` is `'py'`,s or tetra if `algo` is `'tetra'`.
    '''
    if algo == 'py':
        args = {
            'pipeline': None,
            '--png': RAW_INPUT_PATH,
            '--focal-length': 49,
            '--pixel-size': 22.2,
            '--centroid-algo': 'cog',  # 'cog', 'dummy', 'iwcog'
            '--centroid-mag-filter': 5,
            '--database': PY_DATABASE_PATH,
            '--star-id-algo': 'py',  # 'dummy', 'gv', 'py', 'tetra'
            '--angular-tolerance': 0.05,
            '--false-stars': 1000,
            '--max-mismatch-prob': 0.0001,
            '--attitude-algo': 'dqm',  # 'dqm' (Davenport Q), 'triad', 'quest'
            '--print-attitude': ATTITUDE_PATH,
        }
    elif algo == 'tetra':
        args = {
            'pipeline': None,
            '--png': RAW_INPUT_PATH,
            '--fov': 17,
            '--centroid-algo': 'cog',
            '--centroid-filter-brightest': 4,
            '--database': TETRA_DATABASE_PATH,
            '--star-id-algo': 'tetra',
            '--false-stars': 0,
            '--attitude-algo': 'dqm',
            '--print-attitude': ATTITUDE_PATH,
        }
    else:
        raise f"Invalid identification algo {algo}. Must be 'py' or 'tetra'."
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
