import imageio
import subprocess
import os
import pathlib


# TODO:
# [x] Bare bones MVP of database, generation, and identification working
# [ ] Organize this file
# [ ] Solidify interfaces (especially around setting parameters)
# [ ] Ensure all parameters are available
# [ ] Test that all parameters work correctly
# [ ] Add docstrings/documentation
# [ ] Example code/usage
# [ ] Properly bundle into a python package


lost_dir = os.path.dirname(os.path.realpath(__file__))

lost_path = f"{lost_dir}/lost"
raw_input_path = f'{lost_dir}/tmp/raw-input.png'
annotated_input_path = f'{lost_dir}/tmp/input.png'
annotated_output_path = f'{lost_dir}/tmp/annotated_output.png'
attitude_path = f'{lost_dir}/tmp/attitude.txt'
database_path = f'{lost_dir}/tmp/tmp_database.dat'


# runs lost with the provided array of command line arguments
# TODO: investigate return/command line output
def lost(args):
    subprocess.run([lost_path] + args, cwd=lost_dir)


# generates kvector database (required before identifying images)
def database(max_stars=5000, kvector=True, kvector_min_distance=0.2,
             kvector_max_distance=15.0, kvector_distance_bins=10000):
    lost(['database',
          '--max-stars', str(max_stars),
          '--kvector',
          '--kvector-min-distance', str(kvector_min_distance),
          '--kvector-max-distance', str(kvector_max_distance),
          '--kvector-distance-bins', str(kvector_distance_bins),
          '--output', database_path])


def generate(x_resolution=1024, y_resolution=1024, fov=30,
             reference_brightness=100, spread_stddev=1, read_noise_stddev=0.05,
             ra=88, de=7, roll=0, generate_raw=True, generate_annotated=True):

    # general arguments
    args = [
        'pipeline',
        '--generate', '1',
        '--generate-x-resolution', str(x_resolution),
        '--generate-y-resolution', str(y_resolution),
        '--fov', str(fov),
        # this is an unrecognized parameter?
        # '--generate-reference-brightness', str(reference_brightness),
        '--generate-spread-stddev', str(spread_stddev),
        '--generate-ra', str(ra),
        '--generate-de', str(de),
        '--generate-roll', str(roll)
    ]

    # add image generation arguments as applicable
    if generate_raw:
        args.extend(['--plot-raw-input', raw_input_path])
    if generate_annotated:
        args.extend(['--plot-input', annotated_input_path])

    # execute LOST
    lost(args)

    # load image to memory
    raw_result = None
    if generate_raw:
        raw_result = imread(raw_input_path)
        delete_file(raw_input_path)

    annotated_result = None
    if generate_annotated:
        annotated_result = imread(annotated_input_path)
        delete_file(annotated_input_path)

    # return image as array, optionally annotated too
    return (raw_result, annotated_result)


focal_length = 49

centroid_algo = 'cog'  # 'cog', 'dummy', 'iwcog'
star_id_algo = 'py'  # 'dummy', 'gv', 'py', 'tetra'
attitude_algo = 'dqm'  # 'dqm' (Davenport Q), 'triad', 'quest'

centroid_mag_filter = 5

pixel_size = 22.2
angular_tolerance = 0.05

false_stars = 1000

max_mismatch_prob = 0.0001


def identify(image, plot_output=False):
    # save the given image to disk so LOST can use it
    imwrite(raw_input_path, image)

    args = [
        'pipeline',
        '--png', raw_input_path,
        '--focal-length', str(focal_length),
        '--pixel-size', str(pixel_size),
        '--centroid-algo', centroid_algo,
        '--centroid-mag-filter', str(centroid_mag_filter),
        '--database', database_path,
        '--star-id-algo', star_id_algo,
        '--angular-tolerance', str(angular_tolerance),
        '--false-stars', str(false_stars),
        '--max-mismatch-prob', str(max_mismatch_prob),
        '--attitude-algo', attitude_algo,
        '--print-attitude', attitude_path,
        # '--plot-output', annotated_output_path,  # TODO: option to return
    ]

    if plot_output:
        args += ['--plot-output', annotated_output_path]

    # identify image
    lost(args)

    # parse/load attitude file
    with open(attitude_path) as f:
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
    delete_file(attitude_path)
    delete_file(raw_input_path)

    if plot_output:
        result['annotated_output_image'] = imread(annotated_output_path)
        delete_file(annotated_output_path)

    # return attitue information, optionally output
    return result


def imread(path):
    im = imageio.imread(path)
    return im


def imwrite(path, image):
    imageio.imwrite(path, image)


def delete_file(path):
    pathlib.Path(path).unlink()
