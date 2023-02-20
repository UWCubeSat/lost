#!/usr/bin/env bash

# Tetra testing

th="0.5"

# compare actual to expected
function cmpAE() {
  ours="$1"
  actual="$2"

  # echo "ours is $ours, actual is $actual"

  if (( $(echo "$ours < $actual" | bc -l) )); then
    diff=$(echo "$actual - $ours" | bc -l )
  elif (( $(echo "$ours > $actual" | bc -l) )); then
    diff=$(echo "$ours - $actual" | bc -l )
  else
    # Because 0=true, 1=false in bash rip
    return 0
  fi

  if (( $(echo "$diff < $th" | bc -l) )); then
    return 0
  else
    return 1
  fi
}


# params: low and hi
function rand_int {
  lo=$1
  hi=$2
  echo $((lo + (RANDOM % (hi - lo + 1))))
}

# create database if not exists
# test -e tetra-incorrect.dat || ./lost database \
#   --min-mag 7 \
#   --tetra \
#   --tetra-max-angle 12 \
#   --output tetra-incorrect.dat

for _ in $(seq "${1:-1000}"); do
  ra=$(rand_int 0 359)
  de=$(rand_int -89 89)
  roll=$(rand_int 0 359)

  # set -x
  lost_output=$(
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
      --database my-database-small.dat \
      --centroid-algo cog \
      --centroid-mag-filter 5 \
      --star-id-algo tetra \
      --compare-star-ids \
      --attitude-algo dqm \
      --print-attitude
  )
  # set +x
  num_incorrect_stars=$(grep -oP "(?<=starid_num_incorrect )\\d+" <<<"$lost_output")
  num_correct_stars=$(grep -4oP "(?<=starid_num_correct )\\d+" <<<"$lost_output")

  raOurs=$(grep -oP "(?<=attitude_ra )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output")
  deOurs=$(grep -oP "(?<=attitude_de )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output")
  rollOurs=$(grep -oP "(?<=attitude_roll )[-+]?[0-9]*\.?[0-9]+" <<<"$lost_output")

  echo "Real: $ra, $de, $roll"
  echo "Calculated: $raOurs, $deOurs, $rollOurs"

  if cmpAE "$raOurs" "$ra" && cmpAE "$deOurs" "$de" && cmpAE "$rollOurs" "$roll"; then
    echo "TEST: [SUCCESS]: ATT All clear!"
  else
    echo "TEST: [ERROR ATT]: Incorrect attitude calculation!"
    # exit 1
  fi

  if ((num_correct_stars > 0 && num_incorrect_stars == 0)); then
    echo "TEST: [SUCCESS]: All clear!"
  elif ((num_correct_stars == 0)); then
    echo "TEST: [ERROR]: No correct stars identified!"
    echo "Values: $ra, $de, $roll"
    # exit 1
  elif ((num_incorrect_stars > 0)); then
    echo "TEST: [ERROR]: Incorrect stars identified!"
    echo "Values: $ra, $de, $roll"
    # exit 1
  else
    echo "TEST: [ERROR]: Unknown error!"
    exit 1
  fi
done
