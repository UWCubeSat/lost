#!/bin/bash


# Issue #44: Throw an error if there are multiple ways to compute the FOV
# but first, a normally running pipeline to make sure that the syntax doesn't break in the future
./lost pipeline --generate 1 || exit 1
./lost pipeline --generate 1 --pixel-size 5 --focal-length 100 --fov 17 2>/dev/null && exit 1
# also throw an error if we get pixel size or focal length without the other:
./lost pipeline --generate 1 --pixel-size 5 2>/dev/null && exit 1
./lost pipeline --generate 1 --focal-length 100 2>/dev/null && exit 1

# Issue #73: Allow printing compare outputs to terminal
# General test that it only allow non-binary outputs to be printed to terminal
# macOS and Linux have different `script` syntax.
if script -q /dev/null -c "true" > /dev/null 2>&1; then
  run_script() { script -q /dev/null -c "$*"; }
else
  run_script() { script -q /dev/null "$@"; }
fi
run_script ./lost pipeline --generate 1 --centroid-algo cog --compare-centroids - | grep -q centroids_mean_error || exit 1
run_script ./lost pipeline --generate 1 --centroid-algo cog --compare-centroids stdout | grep -q centroids_mean_error || exit 1
# `script` captures stderr by default
run_script ./lost pipeline --generate 1 --plot-raw-input - | grep -q 'WARNING' || exit 1

# Issue #32: Unenlightening error message when neither --generate nor --png is passed to pipeline
./lost pipeline 2>&1 | grep -q ERROR || exit 1

# Comparator assertions
./lost pipeline --generate 1 --plot-output /dev/null 2>&1 | grep -qFe '--plot-output' || exit 1
./lost pipeline --generate 1 --centroid-algo cog --plot-output /dev/null 2>&1 | grep -qFe '--plot-output' && exit 1

# Issue #36: Cog and Attitude without Star-ID
./lost pipeline --generate 1 --centroid-algo cog --attitude-algo triad 2>/dev/null && exit 1

# Run the generator without centroids a whole bunch and make sure no assertions go off
./lost pipeline --generate 200 --generate-perturb-centroids 5 --generate-centroids-only || exit 1

# Error message for database that does not exist
no_database_out=$(./lost pipeline --generate 1 --database 'does not.exist' --star-id-algo py 2>&1) && exit 1
[[ $no_database_out == 'Error reading database!'* ]] || exit 1

# Speed 95-th percentile should be different than max for 20 but not 19 trials
nineteen_out=$(./lost pipeline --generate 19 --generate-centroids-only --attitude-algo quest --print-speed -) || exit 1
nineteen_max_ns=$(echo "$nineteen_out" | grep total_max_ns | cut -d' ' -f2)
nineteen_95_ns=$(echo "$nineteen_out" | grep 'total_95%_ns' | cut -d' ' -f2)
[[ -n "$nineteen_max_ns" && -n "$nineteen_95_ns" ]] || exit 1
(( nineteen_max_ns == nineteen_95_ns )) || exit 1
twenty_out=$(./lost pipeline --generate 20 --generate-centroids-only --attitude-algo quest --print-speed -) || exit 1
twenty_max_ns=$(echo "$twenty_out" | grep total_max_ns | cut -d' ' -f2)
twenty_95_ns=$(echo "$twenty_out" | grep 'total_95%_ns' | cut -d' ' -f2)
[[ -n "$twenty_max_ns" && -n "$twenty_95_ns" ]] || exit 1
(( twenty_max_ns > twenty_95_ns )) || exit 1

echo 'TEST: [SUCCESS]: All clear!'
exit 0
