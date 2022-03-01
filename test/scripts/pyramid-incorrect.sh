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
    echo $((lo + (RANDOM % (hi-lo+1))))
}

# create database if not exists
test -e pyramid-incorrect.dat || ./lost build_database 9999 kvector .5 10 10000 done pyramid-incorrect.dat

for i in $(seq ${1:-10})
do
    x_res=$(rand_int 100 1000)
    y_res=$(rand_int 100 1000)
    fov=$(rand_int 5 70)
    ra=$(rand_int 0 359)
    de=$(rand_int -89 89)
    roll=$(rand_int 0 359)
    tolerance=$(awk "BEGIN{print $RANDOM%10/100+.002}") # floating point numbers are hard, alright?
    set -x
    lost_output=$(./lost pipeline generate 1 $x_res $y_res 0 $fov 0 0 0 $ra $de $roll database pyramid-incorrect.dat starid pyramid $tolerance 100 .001 done compare_stars - done)
    set +x
    num_incorrect_stars=$(grep -oP "(?<=starid_num_incorrect )\\d+" <<< "$lost_output")
    if (( num_incorrect_stars >= 4 ))
    then
	echo "TEST: Too many incorrect stars: $num_incorrect_stars!"
	exit 1
    else
	echo "TEST: All clear."
    fi
done
