from __future__ import absolute_import, print_function, division;
import numpy as np;
from FastProject import FileIO;
from FastProject import Signatures;
from FastProject import Pipelines;
from FastProject import HtmlViewer;
from FastProject import Transforms;
from FastProject import SigScoreMethods;
from FastProject import NormalizationMethods
from FastProject.DataTypes import ExpressionData;
from FastProject.Global import FP_Output, get_housekeeping_dir;
from FastProject import CLI;
import FastProject

file_dict = {"signatures": ["../../data/geneSetLibrary.gmt"],
           "data_file": "../../data/expression_matrix.txt",
           "housekeeping": "../../data/Gene Name Housekeeping.txt",
           "precomputed": ["../../data/precomputed_signatures.txt"],
           "projections": None,
           "weights": None}

files = CLI.loadFilesFromDisk(file_dict)

expressionMatrix = files[0]
signatures = files[1]
housekeeping_genes = files[3]

sig_norm_method = NormalizationMethods.row_normalization;
sig_data = expressionMatrix.get_normalized_copy(sig_norm_method);

for sig in signatures:
    try:
        ss = SigScoreMethods.naive_eval_signature(sig_data, sig, None, 0)
        break
    except ValueError:
        pass 
