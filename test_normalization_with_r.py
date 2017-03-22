from __future__ import absolute_import, print_function, division;
import numpy as np;
from FastProject import FileIO;
from FastProject import Signatures;
from FastProject import Pipelines;
from FastProject import HtmlViewer;
from FastProject import Transforms;
from FastProject.DataTypes import ExpressionData;
from FastProject.Global import FP_Output, get_housekeeping_dir;
from FastProject import CLI;
from FastProject import NormalizationMethods
import FastProject

file_dict = {"signatures": ["../../data/geneSetLibrary.gmt"],
           "data_file": "../../data/expression_matrix.txt",
           "housekeeping": "../../data/Gene Name Housekeeping.txt",
           "precomputed": ["../../data/precomputed_signatures.txt"],
           "projections": None,
           "weights": None}

files = CLI.loadFilesFromDisk(file_dict)

expressionMatrix = files[0]

no_norm = NormalizationMethods.no_normalization(expressionMatrix.base)
col_norm = NormalizationMethods.col_normalization(expressionMatrix.base)
row_norm = NormalizationMethods.row_normalization(expressionMatrix.base)
print(row_norm)
print(row_norm.shape)
row_col_norm = NormalizationMethods.row_and_col_normalization(expressionMatrix.base)
col_rank_norm = NormalizationMethods.col_rank_normalization(expressionMatrix.base)
