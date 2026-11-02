#!/usr/bin/env bash
# End-to-end tests for the Tetra star-id algorithm: database generation, accuracy on generated
# images, a real image, and error handling. The database is regenerated on every run.
set -e

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

fail() {
  echo "TEST: [ERROR]: $1"
  exit 1
}

# Value of a "key value" line in LOST output
get_value() {
  awk -v key="$1" '$1 == key {print $2}' <<<"$2"
}

# Database generation
./lost database \
  --tetra \
  --tetra-max-fov 25 \
  --output "$tmp_dir/tetra.dat" > /dev/null 2>&1 || fail "Tetra database generation failed"

# Generated images (fixed seed, so results are deterministic). Fail if fewer than half are
# identified; some images legitimately have too few stars or no matching pattern.
num_images=20
min_correct=$((num_images / 2))
generated_flags=(
  --generate "$num_images"
  --generate-random-attitudes
  --generate-seed 42
  --generate-x-resolution 1024
  --generate-y-resolution 1024
  --fov 20
  --generate-exposure 1.0
  --generate-spread-stddev 1
  --generate-read-noise-stddev 0.05
  --database "$tmp_dir/tetra.dat"
  --centroid-algo cog
  --centroid-mag-filter 3
  --star-id-algo tetra
  --attitude-algo dqm
  --compare-star-ids
  --compare-attitudes
)

# Default (TRIAD) pairing must never misidentify a star
lost_output=$(./lost pipeline "${generated_flags[@]}" 2>/dev/null)
num_correct=$(get_value starid_num_images_correct "$lost_output")
num_incorrect=$(get_value starid_num_images_incorrect "$lost_output")
error_rate=$(get_value attitude_error_rate "$lost_output")
echo "Tetra (triad): $num_correct/$num_images generated images correct, $num_incorrect incorrect"
((num_incorrect == 0)) || fail "Tetra misidentified stars in $num_incorrect generated images"
((num_correct >= min_correct)) || fail "Tetra identified only $num_correct/$num_images generated images"
[[ $error_rate == 0 ]] || fail "Tetra attitude error rate is $error_rate"

# Centroid-radius pairing can swap stars in symmetric patterns, so only check it still identifies
lost_output=$(./lost pipeline "${generated_flags[@]}" --tetra-pairing centroid-radius 2>/dev/null)
num_correct=$(get_value starid_num_images_correct "$lost_output")
num_incorrect=$(get_value starid_num_images_incorrect "$lost_output")
echo "Tetra (centroid-radius): $num_correct/$num_images generated images correct, $num_incorrect incorrect"
((num_correct >= min_correct)) \
  || fail "Tetra (centroid-radius) identified only $num_correct/$num_images generated images"

# Real image: same attitude as the README Pyramid example
curl -sS https://markasoftware.com/img_7660.png -o "$tmp_dir/img_7660.png"
lost_output=$(
  ./lost pipeline \
    --png "$tmp_dir/img_7660.png" \
    --focal-length 49 \
    --pixel-size 22.2 \
    --centroid-algo cog \
    --centroid-mag-filter 5 \
    --database "$tmp_dir/tetra.dat" \
    --star-id-algo tetra \
    --attitude-algo dqm \
    --print-attitude 2>/dev/null
)
ra=$(awk '/attitude_ra/ {print int($2)}' <<<"$lost_output")
de=$(awk '/attitude_de/ {print int($2)}' <<<"$lost_output")
roll=$(awk '/attitude_roll/ {print int($2)}' <<<"$lost_output")
((ra == 17 && de == 63 && roll == 12)) || fail "Tetra on real image gave attitude ($ra, $de, $roll)"

# Error handling
./lost pipeline --generate 1 --database "$tmp_dir/tetra.dat" --centroid-algo cog \
  --star-id-algo tetra --tetra-pairing bogus > /dev/null 2>&1 \
  && fail "Illegal --tetra-pairing was accepted"

./lost database --max-stars 5000 --kvector --output "$tmp_dir/kvector.dat" > /dev/null 2>&1
no_tetra_out=$(./lost pipeline --generate 1 --database "$tmp_dir/kvector.dat" --centroid-algo cog \
  --star-id-algo tetra 2>&1) || fail "Tetra crashed on a database without a Tetra section"
[[ $no_tetra_out == *'Could not get to Tetra database'* ]] \
  || fail "No error message for database without a Tetra section"

echo "TEST: [SUCCESS]: All clear!"
