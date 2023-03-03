#!/bin/bash

set -x # print commands being run

echo 'Issue #44: Throw an error if there are multiple ways to compute the FOV'
# but first, a normally running pipeline to make sure that the syntax doesn't break in the future
./lost pipeline --generate 1 || exit 1
# not sure why the parens are needed here, but they are
./lost pipeline --generate 1 --pixel-size 5 --focal-length 100 --fov 17 2>/dev/null && exit 1
# also throw an error if we get pixel size or focal length without the other:
./lost pipeline --generate 1 --pixel-size 5 2>/dev/null && exit 1
./lost pipeline --generate 1 --focal-length 100 2>/dev/null && exit 1

exit 0
