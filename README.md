# LOST: Open-source Star Tracker

[![Build, Test, and Lint](https://github.com/UWCubeSat/lost/actions/workflows/build-test-lint.yml/badge.svg)](https://github.com/UWCubeSat/lost/actions/workflows/build-test-lint.yml)
[![Publish Docs](https://github.com/UWCubeSat/lost/actions/workflows/publish-docs.yml/badge.svg)](https://github.com/UWCubeSat/lost/actions/workflows/publish-docs.yml)

LOST is star tracker software for small, low-power, low-cost satellites. It is being developed in
the Husky Satellite Lab, a CubeSat team at the University of Washington.

In depth code documentation is available on this repo's [Doxygen Page](https://uwcubesat.github.io/lost/)

# Building LOST

When making an actual, physical star tracker, you will most likely need to pick out the specific
parts of LOST you want and make lots of changes. However, there's a LOST binary you can use to
test and benchmark things quickly.

## Local installation

- Linux or Mac. If you have Windows, I recommend installing the Windows Subsystem for Linux.
- A C++ compiler, such as g++. On Debian, `apt install g++`
- GNU Make. On Debian, `apt install make`
- Groff, to generate help text. `apt install groff`
- (Recommended) ASAN (Address Sanitizer). We use this to catch memory errors early. `apt install libasan`. (If you wish to use LOST without the address sanitizer, build it with `make LOST_DISABLE_ASAN=1`)
- Git. On Debian, `apt install git`
- Cairo, to read and write PNGs as well as draw on them for debugging and demo purposes. You won't
  need this on your CubeSat. On Debian, `apt install libcairo2-dev`. Elsewhere, follow the
  instructions here: https://www.cairographics.org/download/ under the Binary Packages section.
- Eigen3, a header-only linear algebra and numerical methods library. On Debian, `apt install
  libeigen3-dev`. Elsewhere, download the latest stable release from https://eigen.tuxfamily.org/.
  You can install it system-wide, or just extract it so that the main folder is in
  `vendor/eigen3/Eigen` under the LOST repository.
  
For mac users:

- Download [Homebrew](https://brew.sh/), run `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"` in terminal
- Install [cairo](https://formulae.brew.sh/formula/cairo#default) via homebrew `brew install cairo`
- Install [Eigen]("https://formulae.brew.sh/formula/eigen#default") via homebrew `brew install eigen`
    - Locate Eigen files and move them to `vendor/eigen3/Eigen` under the LOST repository

Clone this repository (`git clone https://github.com/uwcubesat/lost`), then `cd lost`, then
`make` will compile everything. Then you can just run `./lost` and play around with the options!

If you're developing LOST, you need to re-run `make` every time you edit any of the source code
before running `./lost`.

## Using Docker

This option is best for Mac, non-Debian Linux users, or anyone who wants to keep LOST and the development dependencies in a container.

Docker is best if you just want to *run* LOST. If you want to help develop LOST, local installation is preferred.

- Get [Docker](https://www.docker.com/get-started/) onto your system (Docker Desktop, an installation that provides a
  graphical interface, is recommended for most users)
- Clone this repository (`git clone https://github.com/uwcubesat/lost`), then `cd lost`
- Run `docker-compose up` to build and start the container. Use this command or `docker-compose start` any time you want
  to start it up again later.
- Connect to the container with `docker attach {name of docker process}` or use the terminal button in Docker Desktop.
- Once inside the container, cd into the `/lost` directory (if necessary) and use `make` to compile LOST and run it
  with `./lost`
  - Be sure to run `make` again after you edit the source code
- Close the container with `docker-compose stop`. You can also close and remove the container with `docker-compose down`

# Usage

Run LOST via the `./lost` executable followed by the appropriate arguments. Executing`./lost` with no arguments will
bring up a usage guide and options to see possible arguments.

## Generating a false image

Here's an example command line to generate a false image, plot it to `raw-input.png`, and all the stars in the
catalog to `input.png`:

```shell
./lost pipeline \
  --generate 1 \
  --plot-raw-input raw-input.png \
  --plot-input input.png
```

We name the image outputs of this command as `input.png` because they will be used as image inputs to the star
identification pipeline. The above command utilizes many of the default values for the parameters. All of these
parameters can be explicitly set by using the appropriate flags:

```shell
./lost pipeline \
  --generate 1 \
  --generate-x-resolution 1024 \
  --generate-y-resolution 1024 \
  --fov 30 \
  --generate-reference-brightness 100 \
  --generate-spread-stddev 1 \
  --generate-read-noise-stddev 0.05 \
  --generate-ra 88 \
  --generate-de 7 \
  --generate-roll 0 \
  --plot-raw-input raw-input.png \
  --plot-input annotated-input.png
```

## Identifying a real image

Here's how to identify a real image:

1. Download https://markasoftware.com/img_7660.png
2. Generate a *database* named `my-database.dat` by running:

```shell
./lost database \
  --max-stars 5000 \
  --kvector \
  --kvector-min-distance 0.2 \
  --kvector-max-distance 15 \
  --kvector-distance-bins 10000 \
  --output my-database.dat
```

3. Identify the image by running:

```shell
./lost pipeline \
  --png img_7660.png \
  --focal-length 49 \
  --pixel-size 22.2 \
  --centroid-algo cog \
  --centroid-mag-filter 5 \
  --database my-database.dat \
  --star-id-algo py \
  --angular-tolerance 0.05 \
  --false-stars 1000 \
  --max-mismatch-prob 0.0001 \
  --attitude-algo dqm \
  --print-attitude attitude.txt \
  --plot-output annotated-7660.png
```

This will print the attitude (right ascension, declination, roll) to the file `attitude.txt`, and write an annotated
image to `annotated-7660.png`.

When identifying a different image, some parameters must be changed. Most important is the
camera's pixel size and focal length, which quantify how "zoomed-in" the camera is. The pixel size
(22.2μm in the example) is the side length, in micrometers, of each pixel. This is commonly given
on a sensor datasheet. Usually in the 1μm-10μm range (The example uses 22.2 because it was
resized). Next is the focal length, in millimeters, which is usually found in the lens
documentation (49mm in the example). If either of these parameters is wrong, the star
identification has a high probability of mismatch. Any statistical guarantees of the algorithms
break down completely when the camera parameters are not known accurately.

The centroid magnitude filter may also need to be adjusted, depending on the resolution and
noisiness of your images. If the output file has many centroids (red boxes) where there are no
visible stars, then the filter should be increased. If there are many stars without centroids, the
filter should be decreased.

# Parts of a Star Tracking System

- **Undistortion or cropping:** It's critical for captured images to be "flat". Unfortunately, real-world lenses make
  things look a little less than flat. Algorithms can undistort images or simply
  crop out the edges to remove the areas where distortion is the worst.

  **Our framework does...**
  - [ ] Undistortion Routines (probably should happen after centroiding)
  - [ ] Cropping Routines (that keep track of how FOV changes due to crop)
  - [ ] Noise removal (median of images from many angles)
- **Centroiding:** Each star in the photo should be reduced to a single point, with sub-pixel
  accuracy. Auxiliary data, such as star brightness or the likelihood it is a binary star, can be
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
- **Star Identification:** The "main" step of star tracking: Going from a list of star positions (and
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
  - [X] QUEST
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

