#!/usr/bin/env bash

# The pyramid algorithm should /never/ misidentify the main, 4-star pyramid in a generated image
# (when not using centroiding). That's because it only accepts the main pyramid when it is the
# uniquely correct pyramid. Since it knows the centroids exactly, there is no chance of some pyramid
# being identified as a different pyramid.

# In light of this fact, here we test randomly generating a number of images, and ensuring pyramid
# doesn't fuck it up.

# params: low and hi
function rand_int {
  lo=$1
  hi=$2
  echo $((lo + (RANDOM % (hi - lo + 1))))
}

# create database if not exists
test -e pyramid-incorrect.dat || ./lost database \
  --max-stars 9999 \
  --kvector \
  --kvector-min-distance 0.5 \
  --kvector-max-distance 10 \
  --kvector-distance-bins 10000 \
  --output pyramid-incorrect.dat

for _ in $(seq "${1:-100}"); do
  # x_res=$(rand_int 100 1000)
  # y_res=$(rand_int 100 1000)
  x_res=1000
  y_res=1000
  fov=$(rand_int 5 70)
  ra=$(rand_int 0 359)
  de=$(rand_int -89 89)
  roll=$(rand_int 0 359)
  tolerance=$(awk "BEGIN{print $RANDOM%10/100+.002}") # floating point numbers are hard, alright?
  set -x
  lost_output=$(
    ./lost pipeline \
      --generate 1 \
      --generate-x-resolution "$x_res" \
      --generate-y-resolution "$y_res" \
      --fov "$fov" \
      --generate-reference-brightness 0 \
      --generate-spread-stddev 0 \
      --generate-read-noise-stddev 0 \
      --generate-ra "$ra" \
      --generate-de "$de" \
      --generate-roll "$roll" \
      --database pyramid-incorrect.dat \
      --star-id-algo py \
      --angular-tolerance "$tolerance" \
      --false-stars 100 \
      --max-mismatch-prob 0.001 \
      --compare-star-ids
  )
  set +x
  num_incorrect_stars=$(grep -oP "(?<=starid_num_incorrect )\\d+" <<<"$lost_output")
  num_correct_stars=$(grep -oP "(?<=starid_num_correct )\\d+" <<<"$lost_output")
  if ((num_correct_stars > 0 && num_incorrect_stars == 0)); then
    echo "TEST: [SUCCESS]: All clear!"
  elif ((num_correct_stars == 0)); then
    echo "TEST: [ERROR]: No correct stars identified!"
    exit 1
  elif ((num_incorrect_stars > 0)); then
    echo "TEST: [ERROR]: Incorrect stars identified!"
    exit 1
  else
    echo "TEST: [ERROR]: Unknown error!"
    exit 1
  fi
done
