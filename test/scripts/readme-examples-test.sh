#!/usr/bin/env bash

tmp_dir=$(mktemp -d)

# Generated image test
./lost pipeline \
  --generate 1 \
  --plot-raw-input "$tmp_dir/raw-input.png" \
  --plot-input "$tmp_dir/input.png"
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
  --plot-raw-input "$tmp_dir/raw-input.png" \
  --plot-input "$tmp_dir/annotated-input.png"

# Real image test
curl https://markasoftware.com/img_7660.png -o "$tmp_dir/img_7660.png"

./lost database \
  --max-stars 5000 \
  --kvector \
  --kvector-min-distance 0.2 \
  --kvector-max-distance 15 \
  --kvector-distance-bins 10000 \
  --output "$tmp_dir/my-database.dat"

set -x
lost_output=$(
  ./lost pipeline \
    --png "$tmp_dir/img_7660.png" \
    --focal-length 49 \
    --pixel-size 22.2 \
    --centroid-algo cog \
    --centroid-mag-filter 5 \
    --database "$tmp_dir/my-database.dat" \
    --star-id-algo py \
    --angular-tolerance 0.05 \
    --false-stars 1000 \
    --max-mismatch-prob 0.0001 \
    --attitude-algo dqm \
    --plot-output "$tmp_dir/annotated-7660.png" \
    --print-attitude
)
set +x

ra=$(grep -oP "(?<=attitude_ra )[0-9]+" <<<"$lost_output")
de=$(grep -oP "(?<=attitude_de )[0-9]+" <<<"$lost_output")
roll=$(grep -oP "(?<=attitude_roll )[0-9]+" <<<"$lost_output")

if ((ra == 17 && de == 63 && roll == 12)); then
  echo "TEST: [SUCCESS]: All clear!"
else
  echo "TEST: [ERROR]: README example with real image failed!"
  exit 1
fi

rm -rf $tmp_dir