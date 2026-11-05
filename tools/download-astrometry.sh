#!/bin/sh

# Downloads results from astrometry.net for a given submission, to compare
# astrometry's plate-solving against LOST's.

set -e

REFERER="Referer: https://nova.astrometry.net/api/login"
HOST=https://nova.astrometry.net

if [ -z "$1" ]; then
    echo "USAGE: ./download-astrometry.sh <submission-id>"
    echo ""
    echo "The submission ID is shown on your astrometry.net results page"
    echo "(e.g. https://nova.astrometry.net/status/SUBMISSION_ID)"
    exit 1
fi

sub_id=$1

echo "Fetching submission $sub_id..."
sub_json=$(curl -sL "$HOST/api/submissions/$sub_id")

job_id=$(echo "$sub_json" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['jobs'][0])" 2>/dev/null) || {
    echo "Error: could not find a job for submission $sub_id"
    echo "Response: $sub_json"
    exit 1
}

image_id=$(echo "$sub_json" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['images'][0])" 2>/dev/null) || {
    echo "Warning: could not find image ID, skipping original image download"
    image_id=""
}

echo "Found job $job_id, image $image_id"

echo "Fetching calibration info..."
curl -sL "$HOST/api/jobs/$job_id/info" | python3 -m json.tool

mkdir -p "$sub_id"

echo "Downloading files into $sub_id/..."
curl -sL -H "$REFERER" -o "$sub_id/axy.fits" "$HOST/axy_file/$job_id"
curl -sL -H "$REFERER" -o "$sub_id/wcs.fits" "$HOST/wcs_file/$job_id"
curl -sL -H "$REFERER" -o "$sub_id/annotated.jpg" "$HOST/annotated_display/$job_id"

if [ -n "$image_id" ]; then
    curl -sL -H "$REFERER" -o "$sub_id/original.png" "$HOST/image/$image_id"
fi

echo "Done. Files saved to $sub_id/"
ls -lh "$sub_id/"
