#!/bin/sh

# Downloads an image and metadata from astrometry, to check astrometry's (supposedly good) results
# against our (supposedly crappy) ones.

set -e

if [ -z $1 ]
then
    echo 'USAGE: ./download-astrometry.sh <user-image-id>'
    exit 1
fi

if ! command -v convert >/dev/null
then
    echo 'Please install imagemagick for jpeg->png conversion!'
    exit 1
fi

astrometry_host=http://nova.astrometry.net

user_image_page=$(curl ""$astrometry_host/user_images/$1)
axy_link=$astrometry_host$(echo "$user_image_page" | grep -oP "/axy_file/[0-9]+")
jpeg_link=$astrometry_host$(echo "$user_image_page" | grep -oP '(?<="original": ")/image/[0-9]+')

mkdir "$1"
curl -o "$1/axy.fits" "$axy_link"
curl -o "$1/img.jpg" "$jpeg_link"
convert "$1/img.jpg" "$1/img.png"

