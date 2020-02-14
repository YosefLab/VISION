# VISION 2.1.0

Added parameter `sig_gene_threshold` with **changed default behavior**

* Before it was a guideline to filter lowly expressed genes before running VISION
* Now, by default, genes expressed in fewer than 0.1% of cells will be filtered automatically

Bug Fixes:

* Better colors in output when more than 10 categories
* Errors with certain output object accessors
* Selections not saving when loading results server
* Crashes when running on more than 1200 signatures

# VISION 2.0.0

Lots of changes for this version.

However the object structure has been refactored.  This means that any old saved RDS objects will not work with the viewResults function once the code has been upgraded.  To resolve this, either re-create the objects or use the latest v1.x release of VISION.  Incompatible changes like this will be avoided in the future, but a refactor was long needed.

* Added support for surface protein data (e.g., CITE-seq)
    * This data is entered separately using the `proteinData` argument
    * Autocorrelation scores will be computed on protein data
    * Differential expression tests will also test protein data
    * Can view protein vs. protein (FACS-style) in output report
* Improved handling of sparse data
    * When inputing sparse expression data, the processing pipeline will no longer expand this data to dense at any stage
* Additional improvements for performance and memory usage
* New exported pipeline functions for custom workflows
* UI Improvements
    * Can see multiple, different views simultaneously in the output (by unchecking the 'Update All?' changes only affect the selected plot
    * Y-axis for feature histograms can be log-scaled
    * Filtering (FDR) and Export for DE results table
    * Caching of DE results to avoid re-running
    * Additional DE options (subsampling, gene filtering) for reduced run times
    * Improved signature-gene heatmap

Many thanks to @Yanay1 for improvements to the DE results and heatmap!

# VISION 1.1.1

* Bugfixes
    * Fixed issue with selection on high-dpi displays
    * Better Encoding of URI Components
* Removed gzmem as a suggested dependency
* Various style fixes
* Added support for Seurat 3.X objects

# VISION 1.1.0

* Added differential expression testing to the output report
    * Thanks to @Yanay1 (Yanay Rosen)
* Added gene-signature importance calculation
* More flexible projection plotting
* Added UMAP as a projection method

# VISION 1.0.1

* Bugfixes related to caching in the output report
* Change default `cellsPerPartition` value to 10
* Change default `filterThreshold` parameter value to 5% in applyMicroClustering (for consistency with Vision constructor)

# VISION 1.0.0

* Added Layout Options to the output report
* Selected cells can now be saved/loaded
* Addition of "get*"-style accessor methods to the VISION object
* Added an interface to Seurat objects
* Added convenience methods for working with outputs from 10x Genomic's CellRanger
* **Note**: Some changes were made to the Vision object structure and how results are stored.  If updating to 1.0.0, older Vision objects may need to be re-created and re-run before using `viewResults`
