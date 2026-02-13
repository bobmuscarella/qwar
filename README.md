# qwar
An R package for Quantitative Wood Anatomy

The `qwar` package uses as input annotations of microscopic images in the form of .svg (vector image) files.  The package converts these annotations to `sf` objects for spatial analyses.  The functionality of `qwar` includes:
- Initial checks of annotations and conversion to `sf` objects
- Computing distance of features (so far normally vessels) to the cambium line
- Computing a variety of vessel characteristics including area, grouping indices, hydraulically-weighted diameter, theoretical conductivity
- Computing vessel fraction (% area) in a sample
- Computing all parameters in bins of a given distance from the cambium (for time-series analyses)

Future functionality will include other parameters as well models to evaluate changes in QWA parameters in a sample.

Install the package from Github like this:
```{r}
require(devtools)
devtools::install_github("bobmuscarella/qwar")
library(qwar)
```


## News
- v 0.0.7: Added `colorID2` function to deal with new structure of svg file data from the `grImports2` package. Added option for no messages (`verbose=F`) during vessel area calculations.
- v 0.0.6: Changed `thumbnail_check` function to allow for missing `bsf` object.
- v 0.0.5: Changed `read_svg` function to allow different names for temporary svg file.
- v 0.0.4: Added `annotations_to_sf2` function to deal with new structure of svg file data from the `grImports2` package.
- v 0.0.3: Changed subsetting script in `cam_dist` function because it wasn't giving error - not sure what changed...?
- v 0.0.2: Added and edited functionality to deal with whole branch samples (with polygons defining outer cambium and pith)
- v 0.0.1: Initial release.
