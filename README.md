VISION [![Travis-CI Build Status](https://travis-ci.org/YosefLab/VISION.svg?branch=master)](https://travis-ci.org/YosefLab/VISION)
===========
Here we present VISION, a module that can sit downstream of other common analyses such as clustering, dimensionality reduction, and trajectory inference. Specifically, VISION aids in the interpretation of scRNA-seq data, with or without predetermined labels, or stratifications of the data (e.g. clusterings) using the notion of cell-cell similarity maps (as interpreted from some latent space) and biological signatures (functional sets of genes that can be obtained online from, for example, MSigDB). Finally, VISION can evaluate the effect of cell- level meta data, such as library quality, batch, clinical information, or additional experimental readouts (e.g., protein levels from a Cite-Seq experiment). Importantly, the use of VISION can greatly facilitate collaborative projects, as it offers a low- latency interactive report for the end- user, which can be hosted online and viewed on a web browser without the need for installing developer-grade software.


Installing VISION
-----------------------

We recommend installing VISION via github using devtools:

```r
require(devtools)
install_github("YosefLab/VISION")
```

The VISION Pipeline
-----------------------
VISION generally follows the same pipeline from iteration to iteration, where minor differences can be specified via the various parameters in a VISION object. On a typical VISION run:

- If a latent space is not specified, PCA is performed and the top 30 components are retained.
- A KNN graph is constructed from the latent space, named the cell-cell similarity map
- Signature scores are computed using the expression matrix
- Signature local “consistencies” on the cell-cell similarity map are computed using the Geary-C statistic, an auto-correlation statistic.
- An interactive web-based report is generated that can be used to explore and interpret the dataset.

How to run VISION
-----------------------

You can refer to the [vignettes](/vignettes) to run VISION. To note, there is an extra vignette detailing how
to properly interface with [Dynverse](https://github.com/dynverse) for incorporating VISION into your
trajectory inference pipeline.

Sample Output
-------------
![Link to an example output report of ~9,000 CBMC's sequenced with the CITE-seq protocol]("http://s124.millennium.berkeley.edu:7703/")
