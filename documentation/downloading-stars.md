# Downloading a Star Catalog

# Automatically
Run `./download-bright-star-database` from the base directory of LOST to generate
`bright-star-database.tsv`.

# Manually
I recommend using http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=I/50. VizieR is a web
interface for accessing about a bajillion different catalogs of stellar objects. This specific
catalog is the "Brightest Stars Catalog". It contains about 9,000 stars, which are all that any
relatively cheap camera is likely to pick up.

First, go over to the little "Preferences" box on the left. Near the top, set the max rows to
"unlimited" and the output style to `|-separated values`. Make sure that `J2000` is checked and that
the `Decimal` option is selected. This tells the database to convert into a standard stellar
coordinate system.

Then, see the `"target"` box near the top of the center. Write `"0+0"`, `"J2000"`, `360 deg`. This field is
for people who want to only fetch stars from a certain slice of the sky, but we want to download the
whole thing!

In the main panel, you should check the `HR`, `Multiple`, and `VMag` boxes, which correspond to a
unique ID for the star, an indicator whether the star is binary (or more), and star brightness.
Uncheck the others.

Then hit "Submit" and download your stars!

After downloading, remove all the lines from the file that aren't data; there's some header-y stuff
near the top.
