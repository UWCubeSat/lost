# Script to check performance on images in given folder and prevent regression testing
# Go over all images in some folder and check performance, log it to stdout (for now)
# TODO: figure out some storage plan

# CLI:
# test directory name, output log name
# (optional) attitude estimator
# Example usage: python reg.py samples/ log.txt -att_estimator quest

import subprocess
import argparse
import os
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("test_dir", type=str)
parser.add_argument("log", type=str, default="log.txt")
parser.add_argument("-att_estimator", type=str, nargs="?", default="dqm")
args = parser.parse_args()

print(f"Testing images in {args.test_dir}")
print(f"attitude estimator: {args.att_estimator}")
print(f"Logging outputs to {args.log}")


def get_diff(expected, actual):
    """Get element-wise angle difference between expected and actual (what we got)"""
    return [expected[i] - actual[i] for i in range(len(actual))]


output_log = open(args.log, "a+")  # append to end of log, don't overwrite
for img_name in os.listdir(args.test_dir):
    cmd = (
        f"./lost pipeline \
      --png {img_name} \
      --fov 17 \
      --centroid-algo cog \
      --centroid-filter-brightest 6 \
      --star-id-algo tetra \
      --database tetra-algo-12.dat \
      --false-stars 0 \
      --attitude-algo {args.att_estimator} \
      --print-attitude attitude.txt \
      --plot-output annotated-res.png",
    )
    subprocess.call(cmd, shell=True)
    # log attitude.txt
    # Log:
    # which image was tested
    # Output from attitude.txt
    dt = datetime.datetime.now().isoformat()
    output_log.write("----------------------------\n")
    output_log.write(f"{img_name}-{dt}-{args.att_estimator}\n")
    output_log.write("----------------------------\n")
    with open("attitude.txt") as att_log:
        targets = ["attitude_ra", "attitude_de", "attitude_roll"]
        for line in att_log:
            line = line.split(" ")
            if line[0] in targets:
                output_log.write(line[1])

output_log.close()
