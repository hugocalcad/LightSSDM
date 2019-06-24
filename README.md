LightSSDM: A lighter version of SSDM(Stacked species distribution modelling) packages

==================

LightSSDM is a lighter implementation of SSDM package in R, written by Sylvain Schmitt et.al. (https://cran.r-project.org/package=SSDM), the purpose is to execute the creation of species distribution models (SDM) described in the article of the Schmitt SSDM package, besides bugs were corrected.

Some algorithms were changed by others that process extensive data more efficiently in time and computer resources, in others some input data was changed so that processing is more efficient; The algorithms that handle the LightSSDM are: Generalized linear model (GLM), Generalized additive model for large datasets (BAM), Multivariate adaptive regression splines (MARS), Generalized boosted regressions model (GBM), Classification tree analysis (CTA), Random forest (RF), Maximum entropy (MAXENT), Artificial neural network (ANN), and Support vector machines (KSVM).

For a full list of changes see [`NEWS`](./NEWS.md).

# Installation

Please be aware that SSDM package use a lot of dependencies (see [`DESCRIPTION`](./DESCRIPTION))

### Install from Github

You can install the latest version of **LightSSDM** from Github using the [`devtools`](https://github.com/hadley/devtools) package:

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github('hugocalcad/LightSSDM')
```

# Usage

After installing, **LightSSDM** package, you can launch the graphical user interface by typing gui() in the console, it has no changes at the original package.

<center>
[**Click to enlarge**](https://raw.githubusercontent.com/sylvainschmitt/SSDM/master/examples/SSDM.gif)<center>
![Screenshot](https://raw.githubusercontent.com/sylvainschmitt/SSDM/master/examples/SSDM.gif)

# Functionnalities

Has the same structure as SSDM package only it was changed to support other functions, lighter that in the original packages her is the same resume that you haver inthe SSDM package for Sylvain Schmitt.

### Data preparation

* `load_occ`: Load occurrence data
* `load_var`: Load environmental variables

### Modelling main functions

* `modelling`: Build an SDM using a single algorithm
* `ensemble_modelling`: Build an SDM that assembles multiple algorithms
* `stack_modelling`: Build an SSDMs that assembles multiple algorithms and species

### Model main methods

* `ensemble,Algorithm.SDM-method`: Build an ensemble SDM
* `stacking,Ensemble.SDM-method`: Build an SSDM
* `update,Stacked.SDM-method`: Update a previous SSDM with new occurrence data

### Model classes

* `Algorithm.SDM`: S4 class to represent SDMs
* `Ensemble.SDM`: S4 class to represent ensemble SDMs
* `Stacked.SDM`: S4 class to represent SSDMs

### Miscellanous

* `gui`: user-friendly interface for SSDM package
* `plot.model`: Plot SDMs
* `save.model`: Save SDMs
* `load.model`: Load SDMs
