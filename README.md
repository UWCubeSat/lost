# LOST: Open-source Star Tracker

This is a framework for developing high-accuracy star trackers on small, low-power, low-cost
satellites. It is being developed in the Husky Satellite Lab, a cubesat team at the University of
Washington.

# Building LOST
When making an actual, physical star tracker, you will most likely need to pick out the specific
parts of LOST you want and make lots of changes. However, there's a LOST binary you can use to
test and benchmark things quickly.

**Requirements:**
- Linux or Mac. If you have Windows, I recommend installing the Windows Subsystem for Linux.
- A C compiler, such as GCC. On Debian, `apt install gcc`
- GNU Make. On Debian, `apt install make`
- Git. On Debian, `apt install git`
- Cairo, to read and write PNGs as well as draw on them for debugging and demo purposes. You won't
  need this on your cubesat. On debian, `apt install libcairo2-dev`. Elsewhere, follow the
  instructions here: https://www.cairographics.org/download/ under the Binary Packages section.
- Eigen3, a header-only linear algebra and numerical methods library. On debian, `apt install
  libeigen3-dev`. Elsewhere, download the latest stable release from https://eigen.tuxfamily.org/.
  You can install it system-wide, or just extract it so that the main folder is in
  `vendor/eigen3/Eigen` under the LOST repository.

Then, clone this repository (`git clone https://github.com/uwcubesat/lost`), then `cd lost`, then
`make` will compile everything. Then you can just run `./lost` and play around with the options!

If you're developing LOST, you need to re-run `make` every time you edit any of the source code
before running `./lost`.

# Usage (OUTDATED: Will be updated soon)
We are currently in the process of rewriting our command-line interface, which will have more
thorough documentation. The current command-line interface prompts for parameters as needed.
Parameters can also be specified on the command line for non-interactive use.

## Generating a false image
Here's an example command line to generate a false image, plot that image to `raw-input.png`, and
all the stars in the catalog to `input.png`:

```shell
./lost pipeline generate 1 1024 1024 30 8000 1 35 88 7 0 done plot_raw_input raw-input.png plot_input annotated-input.png done
```

You can also just run `./lost` and enter the parameters one at a time, to learn what they are all
for.

## Identifying a real image
Here's how to identify a real image:
1. Download https://markasoftware.com/img_7660.png
2. Generate a /database/ by running:
```shell
./lost build_database 5000 kvector 0.2 15 10000 done my-database.dat
```
3. we can actually identify the image by running:
```shell
./lost pipeline png img_7660.png 22.2 49 \
       centroid cog \
       centroid_magnitude_filter 5 \
       database my-database.dat \
       starid pyramid .05 1000 .0001 \
       attitude dqm \
       done \
       print_attitude - \
       plot_output annotated-7660.png
```

This will print the attitude (right ascension, declination, roll) to standard output, and write an
annotated image to `annotated-7660.png`.

When identifying a different image, some parameters must be changed. Most important is the
camera's pixel size and focal length, which quantify how "zoomed-in" the camera is. The pixel size
(22.2μm in the example) is the side length, in micrometers, of each pixel. This is commonly given
on a sensor datasheet. Usually in the 1μm-10μm range (The example uses 22.2 because it was
resized). Next is the focal length, in millimeters, which is usually found in the lens
documentation (49mm in the example). If either of these parameters is wrong, the star
identification has a high probability of mismatch. Any statistical guarantees of the algorithms
break down completely when the camera parameters are not known accurately.

The centroid magnitude filter may also need to be adjusted, depending on the resolution and
noisy-ness of your images. If the output file has many centroids (red boxes) where there are no
visible stars, then the filter should be increased. If there are many stars without centroids, the
filter should be decreased.


# Parts of a Star Tracking System
- **Undistortion or cropping:** It's critical for captured images to be "flat". Unfortunately, real
  world lenses make things look a little less than flat. Algorithms can undistort images or simply
  crop out the edges to remove the areas where distortion is the worst.

  **Our framework does...**
    - [ ] Undistortion Routines (probably should happen after centroiding)
    - [ ] Cropping Routines (that keep track of how FOV changes due to crop)
    - [ ] Noise removal (median of images from many angles)
- **Centroiding:** Each star in the photo should be reduced to a single point, with sub-pixel
  accuracy. Auxiliary data, such as star brightness or the likelihood it is a binary star, may be
  collected too.

  **Our framework does...**
    - [X] Simple centroiding
    - [ ] Iterative weighted centroiding
    - [ ] 2D Gaussian fit centroiding
    - [ ] Gaussian Grid centroiding
    - [X] Coordinate Conversion (between pixel and spherical/angular)
- **Catalog Building:** This happens on the ground. The format of this catalog depends a lot on the
  Star Identification algorithm used. It might contain information about distances to adjacent
  stars, expected brightness, etc. One thing all catalogs have in common is the actual spherical
  coordinates of the stars, so that once the stars have been identified, the spacecraft's actual
  attitude can be determined.

  **Our framework does...**
    - [X] Downloading a star catalog and converting to an internal format
    - [X] Database-building routines for common Star-ID algorithms.
- **Star Identification:** The "main" step of star tracking: Going from a list star positions (and
  possibly magnitudes or other info)

  **Our framework does...**
    - [ ] Padgett Grid identification method with variations:
        - [ ] Flower Method
        - [ ] Sequential Sum method
    - [X] Pyramid method
    - [X] Geometric Voting
    - [ ] Uncalibrated K-vectory method
    - [ ] Star-ND or Liebe
    - [ ] LIS, Tracking, and Uncalibrated modes
- **Attitude Determination:** Once enough stars have been identified, they can be combined with
  information on camera parameters to determine the attitude.

  **Our framework does...**
    - [X] Davenport Q method
    - [X] TRIAD
    - [ ] QUEST
    - [ ] ESOQ

## Other things LOST can do
  Other parts of our framework that are not essential parts of a star-tracking system:
- [X] Centroiding and Star-ID Results Visualization
- [X] Simulated image generation
- [ ] Re-projection of stars after fix visualization
- [ ] Benchmarking against reference images
- [ ] Demo App that acquires a fix, calibrates the camera, then tracks, with real-time
  visualizations. Would be really cool to make it work on a phone!
- [ ] Automatic corruption of images to test noise tolerance. Star removal, false star adding,
  moon and earth and sun adding, optical offset, focal length mismatch.

