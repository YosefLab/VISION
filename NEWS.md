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
