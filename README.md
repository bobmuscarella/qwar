# qwar
An R package for Quantitative Wood Anatomy

The `qwar` package uses as input annotations of microscopic images in the form of .svg (vector image) files.  The package converts these annotations to `sf` objects for spatial analyses.  The functionality of `qwar` includes:
- Initial checks of annotations and conversion to `sf` objects
- Computing distance of features (so far normally vessels) to the cambium line
- Computing a variety of vessel characteristics including area, grouping indices, hydraulically-weighted diameter, theoretical conductivity
- Computing vessel fraction (% area) in a sample
- Computing all parameters in bins of a given distance from the cambium (for time-series analyses)

Future functionality will include other parameters as well models to evaluate changes in QWA parameters in a sample.
