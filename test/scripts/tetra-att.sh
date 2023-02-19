#!/usr/bin/env bash

tmp_dir=$(mktemp -d)

function rand_int {
  lo=$1
  hi=$2
  echo $((lo + (RANDOM % (hi - lo + 1))))
}

th="0.5"

# compare float to int
function cmpFI() {
  local ours="$1"
  local actual="$2"
  # if ! [[ "$ours" =~ ^-?[0-9]*\.?[0-9]+$ ]]; then
  #   exit 1
  # fi
  # if ! [[ "$actual" =~ ^-?[0-9]+$ ]]; then
  #   exit 1
  # fi

  set +e
  if [[ $(echo "if (${ours}-${actual} < 0) (${actual}-${ours}) else (${ours}-${actual})" | bc) -le "$th" ]]; then
    return 1
  else
    return 0
  fi
  set -e
}

for _ in $(seq "${1:-200}"); do

  ra=$(rand_int 0 359)
  de=$(rand_int -89 89)
  roll=$(rand_int 0 359)
  # set -x
  set +x

  ./lost pipeline \
    --generate 1 \
    --generate-x-resolution 1024 \
    --generate-y-resolution 1024 \
    --fov 12 \
    --generate-reference-brightness 100 \
    --generate-spread-stddev 1 \
    --generate-read-noise-stddev 0.05 \
    --generate-ra "$ra" \
    --generate-de "$de" \
    --generate-roll "$roll" \
    --plot-raw-input "$tmp_dir/raw-input.png" \


  # set +x

  lost_output=$(
    ./lost pipeline \
      --png "$tmp_dir/raw-input.png" \
      --fov 12 \
      --centroid-algo cog \
      --centroid-mag-filter 5 \
      --database my-database-small.dat \
      --star-id-algo tetra \
      --angular-tolerance 0.05 \
      --false-stars 1000 \
      --attitude-algo dqm \
      --print-attitude)

  # echo "LOST OUTPUT: $lost_output"

  # raOurs=$(grep -oP "(?<=attitude_ra )[0-9]+" <<<"$lost_output")
  # deOurs=$(grep -oP "(?<=attitude_de )[0-9]+" <<<"$lost_output")
  # rollOurs=$(grep -oP "(?<=attitude_roll )[0-9]+" <<<"$lost_output")

  # raOurs=$(grep -oP "(?<=attitude_ra )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output" | awk '{print int($1 + 0.5)}')
  # deOurs=$(grep -oP "(?<=attitude_de )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output" | awk '{print int($1 + 0.5)}')
  # rollOurs=$(grep -oP "(?<=attitude_roll )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output" | awk '{print int($1 + 0.5)}')

  raOurs=$(grep -oP "(?<=attitude_ra )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output")
  deOurs=$(grep -oP "(?<=attitude_de )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output")
  rollOurs=$(grep -oP "(?<=attitude_roll )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output")

  echo "Real: $ra, $de, $roll"
  echo "Calculated: $raOurs, $deOurs, $rollOurs"

  # if (( raOurs == ra && deOurs == de && rollOurs == roll )); then
  if cmpFI "raOurs" "ra" && cmpFI "deOurs" "de" && cmpFI "rollOurs" "roll"; then
    echo "TEST: [SUCCESS]: All clear!"
  else
    echo "TEST: [ERROR]: Incorrect attitude calculation!"
    # exit 1
  fi

  # set +x


done

rm -rf $tmp_dir
