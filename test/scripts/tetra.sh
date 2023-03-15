#!/usr/bin/env bash

# Tetra testing

th="0.5"

totalTime=0
# CHANGE THIS, n = number of images to try
n=1000
totalCorrect=0
totalIncorrect=0
totalNoStars=0

totalAttCorrect=0
totalAttIncorrect=0


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

  diff=$(echo "$diff % 360" | bc)

  if (( $(echo "$diff < $th" | bc -l) || $(echo "(360 - $diff) < $th" | bc -l) )); then
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

# This works fine
# create database if not exists
# test -e tetra-incorrect.dat || ./lost database \
#   --min-mag 7 \
#   --tetra \
#   --tetra-max-angle 12 \
#   --output tetra-incorrect.dat

for _ in $(seq "${1:-$n}"); do
  ra=$(rand_int 0 359)
  de=$(rand_int -89 89)
  roll=$(rand_int 0 359)
  fov=$(rand_int 10 60)

  # ra=141
  # de=-48
  # roll=122

  # TODO: separate image generation and actual pipeline run for accurate timing
  # set -x
  lost_output=$(
    ./lost pipeline \
      --generate 1 \
      --generate-x-resolution 1024 \
      --generate-y-resolution 1024 \
      --fov "$fov" \
      --generate-reference-brightness 100 \
      --generate-spread-stddev 1 \
      --generate-read-noise-stddev 0.05 \
      --generate-ra "$ra" \
      --generate-de "$de" \
      --generate-roll "$roll" \
      --database my-database-small-3.dat \
      --centroid-mag-filter 5 \
      --star-id-algo tetra \
      --compare-star-ids \
      --attitude-algo dqm \
      --print-attitude
  )

  # --centroid-algo cog \
  # set +x
  num_incorrect_stars=$(grep -oP "(?<=starid_num_incorrect )\\d+" <<<"$lost_output")
  num_correct_stars=$(grep -4oP "(?<=starid_num_correct )\\d+" <<<"$lost_output")

  raOurs=$(grep -oP "(?<=attitude_ra )[-+]?[0-9]*\.?[0-9]*[eE]?[-+]?[0-9]*" <<<"$lost_output")
  deOurs=$(grep -oP "(?<=attitude_de )[-+]?[0-9]*\.?[0-9]*[eE]?[-+]?[0-9]*" <<<"$lost_output")
  rollOurs=$(grep -oP "(?<=attitude_roll )[-+]?[0-9]*\.?[0-9]*[eE]?[-+]?[0-9]*" <<<"$lost_output")

  echo "Fov: $fov"
  echo "Real: $ra, $de, $roll vs Calculated: $raOurs, $deOurs, $rollOurs"

  # TODO: 360 and 0 should be deemed equivalent

  if cmpAE "$raOurs" "$ra" && cmpAE "$deOurs" "$de" && cmpAE "$rollOurs" "$roll"; then
    echo "TEST: [SUCCESS]: ATT All clear!"
    totalAttCorrect=$(( $totalAttCorrect + 1 ))
  else
    echo "TEST: [ERROR ATT]: Incorrect attitude calculation!"
    totalAttIncorrect=$(( $totalAttIncorrect + 1 ))
  fi

  if ((num_correct_stars > 0 && num_incorrect_stars == 0)); then
    echo "TEST: [SUCCESS]: All clear!"
    totalCorrect=$(( $totalCorrect + 1 ))
  elif ((num_correct_stars == 0)); then
    echo "TEST: [ERROR]: No correct stars identified!"
    # echo "Values: $ra, $de, $roll"
    totalNoStars=$(( $totalNoStars + 1 ))
    # exit 1
  elif ((num_incorrect_stars > 0)); then
    echo "TEST: [ERROR]: Incorrect stars identified!"
    # echo "Values: $ra, $de, $roll"
    totalIncorrect=$(( $totalIncorrect + 1 ))
    # exit 1
  else
    echo "TEST: [ERROR]: Unknown error!"
    exit 1
  fi
done

echo "Total n = $n"
echo "Total correct: $totalCorrect"
echo "Total incorrect: $totalIncorrect"
echo "Total no stars: $totalNoStars"

echo "Total attitude correct: $totalAttCorrect"
echo "Total attitude incorrect: $totalAttIncorrect"
