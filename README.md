FastProjectR
===========
Here we present FastProjectR, an extension of FastProject for R, enabling the two-dimensional visualization of single-cell data. The software package utilizes a wide spectrum of projection methods and provides functionality for systematically investigating the biological relevance of these low dimensional represenations by incorporating domain knowledge. While preserving the spirit of the first version, this R extension most importantly now supports the analysis of massive amounts of cell data - theoretically up to a million cells sequenced at any depth - and an improved visualization platform that allows for subset analysis based on the two dimensional representation.

FastProjectR analyzes a gene expressin matrix and produces a dynamic output report in which two-dimensional projections of the         data can be explored. Annotated gene sets (referred to as "signatures") are incorporated so that features in the projections can be understood in relation to the biological processes they might represent. FastProjectR provides a novel method of scoreing each cell against a gene signature so as to minimize the effect of missed transcripts as well as a method to rank signature-projection pairings so that meaningful associations can be quickly identified. Additionally, FastProjectR is written with a modular architecture and designed to serve as a platform for incorporating and comparing new projection methods and gene selection algorithms.

Installing FastProjectR
-----------------------

FastProjectR is currently not hosted on CRAN nor Bioconductor because it is still in development. To install and use the package, first ensure that an R version >= 3.4 is installed and that the package **devtools** is installed.

To simplify the install process, first run the following command from within an R session whose working direcory is that where the pacakge is located:
``` r
  devtools::install(dependencies=T)
  ```

  This command will automatically try to install all dependencies required by the packages and will in addition attempt to load the package into your R session's namespace. Most likely, however, there will be errors accompanying the installation - these packages that cannot be installed will need to be installed from github and bioconductor. To install the Rtsne.multicore package, use:
  ``` r
  devtools::install_github("RGLab/Rtsne.multicore") 
  ```
  The rest of the packages will most likely need to be installed from Bioconductor using biocLite - use these commands in an R session to install a package from Bioconductor:
  ```r
  source("http://bioconductor.org/biocLite.R")
  biocLite(<name of the package>)
  ```
  How to Run FastProjectR
  -----------------------

  Once all necessary package dependencies are installed, to load FastProjectR into an R session, you can use
  ```r
  devtools::load_all()
  ```
  Note that if you are not running R from within the package directory, you must provide the path link as an argument to this function. This function will simulate loading in a package, as would happen if a package were being loaded in which was already installed on your computer. To run FastProjectR, a minimum of three files must be provided: an expression matrix (.txt), a list of housekeeping genes (.txt) , and list of signatures (.gmt or .txt). Refer to the packages wiki for more details on the input files. With these three files, a FastProjectR object can be created:
  ```r 
  fp <- FastProject(<expression file>, <housekeeping genes>, <list of signatures>)
  ```
  Note that the signature files must be provided in a list structure (even if you are only inputting one file) -- this is to enable the user to provide multiple signature files without having to collate them all first. There are many other arguments that can be provided to the FastProjectR object upon initialization, all of which can be found in the wiki or documentation. With the FastProject object created, the analysis can be run with 
  ``` r
  fpo <- Analyze(fp)
  ```

  This anlaysis will return a FastProjectOutput object which stores all relevant information for the visualization of the single cell data; this object can also be saved in an RDS file for future analysis. 

  Analyzing FastProjectR Results
  ------------------------------
  With the returned FastProjectOutput object, a dynamic report can be generated with the following command:
  ```r 
  viewResults(fpo)
  ```
  Where **fpo** is the FastProjectOutput object generated from analysis. This will launch a local server on which the report will be shown; for best results, use Google Chrome. Alternatively, the FastProjectOutput object can be saved first as an RDS object and then analyzed via the online report:
  ```r 
  saveFPOutAndViewResults(fpo)
  ```

  Sample Output
  -------------
  ![FastProjectR Output Sample Image](/SampleOutput.png?raw=true)

  Additional Info
  ---------------
  All additional documentation is in the [FastProjectR Wiki](https://github.com/YosefLab/FastProjectR/wiki)
