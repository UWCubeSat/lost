./lost pipeline \
  --generate 1 \
  --plot-raw-input raw-input.png \
  --plot-input input.png
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
  --plot-raw-input raw-input.png \
  --plot-input annotated-input.png

wget https://markasoftware.com/img_7660.png

./lost database \
  --max-stars 5000 \
  --kvector \
  --kvector-min-distance 0.2 \
  --kvector-max-distance 15 \
  --kvector-distance-bins 10000 \
  --output my-database.dat
./lost pipeline \
  --png img_7660.png \
  --focal-length 49 \
  --pixel-size 22.2 \
  --centroid-algo cog \
  --centroid-mag-filter 5 \
  --database my-database.dat \
  --star-id-algo py \
  --angular-tolerance 0.05 \
  --false-stars 1000 \
  --max-mismatch-prob 0.0001 \
  --attitude-algo dqm \
  --print-attitude attitude.txt \
  --plot-output annotated-7660.png