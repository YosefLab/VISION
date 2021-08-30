# GET: /SessionInfo
Returns high-level configuration info on this VISION session. Mostly used
to enable/disable features in the visualizer.

- 'name': A title for this session (shown in the top title bar)
- 'meta\_sigs': Names of meta-data signatures.  Sent here as a convenience as
numeric meta-data signatures and gene signatures get combined in the display
and often lookups are needed to distinguish.
- 'ncells': Total number of cells
- 'pooled': Whether or not there is micro-pooling
- 'has\_sigs': Whether or not there are gene signatures
- 'has\_proteins': Whether or not there are CITE-seq proteins
- 'has\_lca': Whether or not LC Annotator analysis is available
- 'has\_tree': Whether or not Trajectories are available

Response structure:
```
    {
        'name': 'Session Name',
        'meta_sigs': [list of string],
        'ncells': int,
        'pooled': bool,
        'has_sigs': bool,
        'has_proteins': bool,
        'has_lca': bool,
        'has_tree': bool,
    }
```

# GET: /Signature/Scores/<sig_name>
Return the per-cell signature scores for signature: `<sig_name>`
Data comes from obj@SigScores matrix

Response structure:

```
    {
        'cells': [list of str,  cell ids...],
        'values': [list of Number, per cell signature scores],
    }
```

# GET: /Signature/Meta/<sig_name>
Return the per-cell meta-data values for variable: `<sig_name>`
Data comes from obj@metaData dataframe

Response structure:

```
    {
        'cells': [list of str,  cell ids...],
        'values': [list of Number or string, per-cell meta-data values],
    }
```

# GET: /Signature/Info/<sig_name>
Returns information associated with a signature
Data comes from object@sigData

Response structure:

```
    {
        'sigDict': { # signature-gene weights
            'GeneName1': 1, 
            'GeneName2': -1,
            ...
        },
        'name': 'Signature Name',
        'source': 'Signature source', # typically filename where signature was loaded from
        'metaData': 'Signature description...', # second column in .gmt file
        'geneImportance': { # signature-gene importance values
            'GeneName1': .1,
            'GeneName2': .01, 
            ...
        }
    }
```

# GET: /Signature/Expression/<sig_name>/<cluster_var>
Returns expression values for `<sig_name>` grouped by `<cluster_var>`.
This is used to generate the gene x cluster heatmap in the bottom left
Values are computed on the fly at the server

Response structure
```
{
    'data': list of list of Number  (2d matrix, outer index is row index)
    'sample_labels': column labels for the matrix in 'data'
    'gene_labels': row labels for the matrix in 'data'
}
```

# GET: /Proteins/\<protein\_name\>/Values
Returns the per-cell protein expression values for `<protein_name>`
Data comes from obj@proteinData matrix

Response structure:
```
    {
        'cells': [list of str,  cell ids...],
        'values': [list of Number, per cell protein expression],
    }
```

# GET: /FilterGroup/SigClusters/Normal
For each signature, describes the assignment of signatures to signature-clusters
This is used to group signatures in the upper-left section

Response structure:
```
{
    "Signature Name" (str): <cluster_number> (int),
    ...
}
```

# GET: /FilterGroup/SigClusters/Proteins
For each protein, describes the assignment of proteins to protein-clusters
Proteins aren't actually clustered in the interface.  Every protein is assigned to cluster #1

Response structure:
```
{
    "Protein Name" (str): 1,
    ...
}
```

# GET: /FilterGroup/SigClusters/Meta
For each meta-data variable, describes the assignment of each variable to meta-data clusters
This is used when micropools are created to group factor variables with their expanded
percent-value-per-micropool representations

Response structure:
```
{
    "Metadata variable Name" (str): <cluster_number> (int),
    ...
}
```

# GET: /Projections/<proj_name>/coordinates/<proj_col>
Returns the coordinates for projection `<proj_name>` in column `<proj_col>`
Since the scatter plot can show more than just projects, this has been expanded to also
return coordinates from the latent space or other meta-data variables.  `<proj_name>` can
come from names in object@Projections, or can be "Meta Data" to refer to object@metaData or
"Latent Space" to refer to object@LatentSpace.  In either case `<proj_name>` should refer
to an associated column.

Response structure:

List of list of `[coordinate, 'cell id']`

```
[
    [coordinate (float): 'cell id' (str)],
    [coordinate (float): 'cell id' (str)],
    [coordinate (float): 'cell id' (str)],
    ...
[
```

# GET: /Projections/list
Returns information on projections and other variables that can be plotted in the main
scatter-plot panel.  Used to create the dropdown menus.

Response structure:

Dictionary where keys are the name of projections/meta data/latent space and the values are
a list of the associated columns

```
{
    'projection_name_1': list of column names in the projection,
    'projection_name_2': list of column names in the projection,
    'Meta Data': list of numeric meta-data variables,
    'Latent Space': list of column names in the latent space
}
```

# GET: /PearsonCorr/Normal
Gets information for the LC Annotator - the correlation of each signature with each latent component

Response structure:
```
{
    "zscores": list of list of correlation coefficients (N signatures x N components),
    "pvals": list of list of p-values (N signatures x N components),
    "proj_labels": list of latent space column names in the same order as columns of 'zscores'/'pvals',
    "sig_labels": list of signature names in the same order as columns of 'zscores'/'pvals',
}
```

# GET: /PearsonCorr/Proteins
Gets information for the LC Annotator - the correlation of each protein with each latent component

Response structure:
```
{
    "zscores": list of list of correlation coefficients (N proteins x N components),
    "pvals": list of list of p-values (N proteins x N components),
    "proj_labels": list of latent space column names in the same order as columns of 'zscores'/'pvals',
    "sig_labels": list of proteins names in the same order as columns of 'zscores'/'pvals',
}
```

# GET: /PearsonCorr/Meta
Gets information for the LC Annotator - the correlation of each numeric meta-data variable with each latent component

Response structure:
```
{
    "zscores": list of list of correlation coefficients (N variables x N components),
    "pvals": list of list of p-values (N variables x N components),
    "proj_labels": list of latent space column names in the same order as columns of 'zscores'/'pvals',
    "sig_labels": list of variables names in the same order as columns of 'zscores'/'pvals',
}
```

# POST: /DE
Requests a differential expression analysis from the server.

Determines the set of cells used for the numerator (i.e. _num or _n) and for the denominator (i.e. _denom or _d)

Examples:
- To select cells in cluster 5 for the numerator set `'type_n': 'meta'`, `'subtype_n': 'cluster'`, `'group_num': 'Cluster 5'`
- To select cells in a saved manual selection for the numerator set `'type_n': 'saved_selection'`, `'group_num': 'name_of_selection'`
- To select a specific set of cell ids for the numerator set `'type_n': 'current'`, `'group_num': list of cell ids`

Posted data structure:
```
{
    'type_n': either 'current', 'meta', or 'saved_selection'
    'type_d':  either 'remainder', 'meta', or 'saved_selection'
    'subtype_n': which meta-data variable to group cells by if type_n is 'meta'
    'subtype_d': which meta-data variable to group cells by if type_d is 'meta'
    'group_num': list of cell ids (if type_n is 'current') OR name of saved selection (if type_n is 'saved_selection')
                 OR value of categorical meta-data variable (if type_n is 'meta')
    'group_denom': same behavior as 'group_num'
    'min_cells': (int) pre-filter genes expressed in less than this many cells
    'subsample_groups': (bool), whether or not to subsample each comparison group first
    'subsample_N': (int) size of the groups after subsampling.  Ignored if 'subsample_groups' is false
}
```

Response structure:

DE Results table with {'column name' -> list of values}

```
{
    'Feature': list of feature (gene or protein) names
    'Type': list of either 'Gene' or 'Protein' 
    'logFC': list of log-fold-change values
    'stat': list of AUC values
    'pval': list of FDR-corrected p-values
}
```


# GET: /Clusters/MetaLevels
The display of the top-left area for signatures/proteins/meta-data depends on the selected
grouping variable in the dropdown.  Commands scoped to a particular grouping variable are
nested under "/Clusters".

This call returns a list of grouping variables and the levels/categories associated with them.

Response structure:
```
    {
        'name of variable': [list of str - values this variable can take],
        'name of variable': [list of str - values this variable can take],
        'Cell Type': ['CD4+', 'CD8+', 'NK'],   # For example
    }
```

# GET: /Clusters/<cluster_variable>/SigProjMatrix/Normal
Returns matrices consisting of test statistic and p-value for 1 vs. all differential signature
tests and for local autocorrelation.

First column is "Score" and values are for local autocorrelation. Statistic is the 1 - Geary's C.

Other columns are for 1 vs. all differential signature tests for each value in the selected grouping
variable in the dropdown.  Test statistic are the AUC values from a ranksums test.  P-values also
from this test.

Response structure:
```
    {
        'sig_labels': [list of str]  # Rows of matrices in zscores and pvals
        'proj_labels': [list of str] # Columns of the matrices in zscores and pvals
        'zscores': [list of list of number]  # signatures x cell group (+ 'Score') matrix of test statistics
        'pvals': [list of list of number] # signatures x cell group (+ 'Score') matrix of p-values
    }
```

# GET: /Clusters/<cluster_variable>/SigProjMatrix/Meta
Same as `/Clusters/<cluster_variable>/SigProjMatrix/Normal`, but for meta-data variables

Response structure:
```
    {
        'sig_labels': [list of str]  # Rows of matrices in zscores and pvals
        'proj_labels': [list of str] # Columns of the matrices in zscores and pvals
        'zscores': [list of list of number]  # variables x cell group (+ 'Score') matrix of test statistics
        'pvals': [list of list of number] # variables x cell group (+ 'Score') matrix of p-values
    }
```
# GET: /Clusters/<cluster_variable>/ProteinMatrix
Same as `/Clusters/<cluster_variable>/SigProjMatrix/Normal`, but for CITE proteins

Response structure:
```
    {
        'sig_labels': [list of str]  # Rows of matrices in zscores and pvals
        'proj_labels': [list of str] # Columns of the matrices in zscores and pvals
        'zscores': [list of list of number]  # proteins x cell group (+ 'Score') matrix of test statistics
        'pvals': [list of list of number] # proteins x cell group (+ 'Score') matrix of p-values
    }
```

# GET: /Clusters/list
Not used anymore

# GET: /Clusters/<cluster_variable>/Cells
Not used anymore


# GET: /Expression/Genes/List
Retrieves a list of genes for which expression values can be plotted.
Used to populate the dropdown in the 'Genes' tab.

Response structure:
```
[list of str - gene names]
```
# GET: /Expression/Gene/<gene_name>
Gets the per-cell expression of the gene specified by `<gene_name>`.
Expression values returned are on a log2 scale.

Response structure:
Items in 'values' correspond to the cells in 'cells' in order
```
    {
        'cells': [list of str - cell identifiers],
        'values': [list of number - gene expression values],
    }
```

# GET: /Cell/<cell_id>/Meta
Retrieves meta-data information for an individual cell specified by `<cell_id>`

Response structure:
```
    {
        'variable name': value (str or number),
        'variable name': value (str or number),
        'variable name': value (str or number),
        ...
    }
```

# POST: /Cells/Meta
Retrieves meta-data information for a group of cells

Posted data structure
```
[list of cell id]
```

Response structure:
```
    {
        'numeric': {
            '<variable name>': {
                'Min': <min value>,
                'Median': <median value>,
                'Max': <max value>,
            },
            '<variable name>': {
                'Min': <min value>,
                'Median': <median value>,
                'Max': <max value>,
            },
            # Repeat for additional numeric meta-data variables
        },
        'factor': {
            '<variable name>': {
                '<factor level 1>': number (proportion of cells assigned to this level,
                '<factor level 2>': number (proportion of cells assigned to this level,
                '<factor level 3>': number (proportion of cells assigned to this level,
                ...
            }
            # Repeat for additional factor meta-data variables
        }
    }
```

# GET: /Cells/Selections
Used when loading a selection of cells.

Response structure:
```
[list of str - names of saved selections]
```

# GET: /Cells/Selections/<selection_id>
Retrieves the cells associated with the selection, `<selection_id>`

Response structure:
```
[list of str - cell ids in the indicated cell selection]
```
# POST: /Cells/Selections/<selection_id>
Saves a cell selection with name `<selection_id>`

```
[list of str - cell ids to be assigned to this cell selection]
```

# GET: /Tree/Projections/list
Get the list of different tree layouts that are available
Used to populate the 'Projection' dropdown under the 'Trajectories' view

R version: returns the names of object@TrajectoryProjections

Response structure:
```
[list of str - names of tree projections]
```

# GET: /Tree/Projections/<proj_name>/coordinates
Returns the coordinates for the tree layout `<proj_name>` in addition
to the position of individual cells

Response structure:
```
    [
        [
            [<x val>, <y val>, <cell id>],  # coordinates for individual cells
            [<x val>, <y val>, <cell id>],
            [<x val>, <y val>, <cell id>],
            ...
        ],
        [
            [<x val>, <y val>],  # coordinates for internal tree nodes
            [<x val>, <y val>], 
            [<x val>, <y val>], 
            ...
        ],
        [list of list of number]  # 0/1 adjacency matrix for internal tree nodes
    ]
```


# GET: /Tree/SigProjMatrix/Normal
Returns matrices consisting of test statistic and p-value trajectory autocorrelation
Format is similar to `/Clusters/<cluster_variable>/SigProjMatrix/Normal`, but without
the 1 vs all test results.

'zscores' and 'pvals' matrices have just a single column

Response structure:
```
    {
        'sig_labels': [list of str]  # Rows of matrices in zscores and pvals
        'proj_labels': ["Score"],
        'zscores': [list of list of number]  # signatures x 1 matrix of test statistics
        'pvals': [list of list of number] # signatures x 1 matrix of p-values
    }
```

# GET: /Tree/SigProjMatrix/Meta
Same as `/Tree/SigProjMatrix/Normal` but for meta-data variables

'zscores' and 'pvals' matrices have just a single column

Response structure:
```
    {
        'sig_labels': [list of str]  # Rows of matrices in zscores and pvals
        'proj_labels': ["Score"],
        'zscores': [list of list of number]  # variables x 1 matrix of test statistics
        'pvals': [list of list of number] # variables x 1 matrix of p-values
    }
```

# GET: /Tree/ProteinMatrix
Same as `/Tree/SigProjMatrix/Normal` but for CITE protein expression values

'zscores' and 'pvals' matrices have just a single column

Response structure:
```
    {
        'sig_labels': [list of str]  # Rows of matrices in zscores and pvals
        'proj_labels': ["Score"],
        'zscores': [list of list of number]  # proteins x 1 matrix of test statistics
        'pvals': [list of list of number] # proteins x 1 matrix of p-values
    }
```
