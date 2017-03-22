from __future__ import print_function, division;
import os;
import logging;
from ..Pipelines import Analysis;
from ..CLI import loadFilesFromDisk, createOutputDirectories;
from ..FileIO import saveResultstoDisk;

this_directory = os.path.dirname(os.path.abspath(__file__));
TEST_FILE_DIR = os.path.join(this_directory, "TestFiles");


def get_default_args():
    args = {
        "data_file": os.path.join(TEST_FILE_DIR, "smallData.txt"),
        "signatures": [os.path.join(TEST_FILE_DIR, "sigsSmall.gmt")],
        "precomputed": [os.path.join(TEST_FILE_DIR, "precomputed_sigs.txt")],
        "housekeeping": "",
        "output": "",
        "nofilter": False,
        "nomodel": False,
        "pca_filter": False,
        "qc": False,
        "subsample_size": None,
        "min_signature_genes": 5,
        "projections": [],
        "weights": None,
        "threshold": None,
        "all_sigs": False,
        "debug": False,
        "lean": False,
        "sig_norm_method": "znorm_rows",
        "sig_score_method": "weighted_avg"}

    return args;


def test(func):

    def test_func():
        args = func()
        try:
            run_test(args);
        except:
            with open('errors.log', 'a') as fout:
                fout.write("Errors Found During Test: " + func.__name__ + "\n");
                print("Errors Found During Test: " + func.__name__ + "\n");
                fout.write("  ARG DUMP:\n");
                for key in args:
                        fout.write("    " + key + ": " + str(args[key]) + "\n");

                import traceback
                fout.write("  EXCEPTION:")
                for line in traceback.format_exc().splitlines():
                    fout.write("    " + line + "\n");
        else:
            with open('errors.log', 'a') as fout:
                fout.write("Test Passed: " + func.__name__ + "\n");
                print("Test Passed: " + func.__name__ + "\n");

    return test_func;


def run_test(args):

    dir_name = createOutputDirectories(args);  # Needs to be created first so logging can write here

    (expressionMatrix, signatures, precomputed_signatures,
     housekeeping_genes, input_projections,
     input_weights) = loadFilesFromDisk(args);

    models, qc_info = Analysis(expressionMatrix, signatures, precomputed_signatures,
        housekeeping_genes, input_projections, input_weights, args);

    saveResultstoDisk(models, signatures, qc_info, dir_name);

    # Cleanup
    import shutil
    import time

    # Close logger
    logger = logging.getLogger("FastProject")
    for handler in logger.handlers:
        handler.close()
        logger.removeHandler(handler)

    for x in range(10):  # Solution to Dropbox locking files.
        try:
            shutil.rmtree(dir_name);
        except:
            pass;
        else:
            break;
        time.sleep(0.1);


@test
def test_simple():
    args = get_default_args();

    return args;


@test
def test_hk():
    args = get_default_args();

    args["housekeeping"] = os.path.join(TEST_FILE_DIR, "housekeeping.txt");

    return args;


@test
def test_text_sigs():
    args = get_default_args();

    args["signatures"] = [os.path.join(TEST_FILE_DIR, "sigsA.txt")];

    return args;


@test
def test_text_sigsB():
    args = get_default_args();

    args["signatures"] = [os.path.join(TEST_FILE_DIR, "sigsB.txt")];

    return args;


@test
def test_multi_sigs():
    args = get_default_args();

    args["signatures"] = [os.path.join(TEST_FILE_DIR, "sigsA.txt"),
                       os.path.join(TEST_FILE_DIR, "sigsSmall.gmt")];

    return args;


@test
def test_no_precomputed():
    args = get_default_args();

    args["precomputed"] = [];

    return args;


@test
def test_no_sigs():
    args = get_default_args();

    args["signatures"] = [];

    return args;


@test
def test_qc():
    args = get_default_args();

    args["qc"] = True

    return args;


@test
def test_nofilter():
    args = get_default_args();

    args["nofilter"] = True;

    return args;


@test
def test_nomodel():
    args = get_default_args();

    args["nomodel"] = True;

    return args;


@test
def test_pca_filter():
    args = get_default_args();

    args["pca_filter"] = True

    return args;


@test
def test_weights():
    args = get_default_args();

    args["weights"] = os.path.join(TEST_FILE_DIR, "weights.txt");

    return args;


@test
def test_projections():
    args = get_default_args();

    args["projections"] = [os.path.join(TEST_FILE_DIR, "input_projection_good.txt")];

    return args;


@test
def test_many_projections():
    args = get_default_args();

    projection_files = ["input_projection_good.txt",
        "input_projection_RFormat_ok.txt",
        "input_projection_header_ok.txt",
        "input_projection_header_ok2.txt"];

    args["projections"] = [os.path.join(TEST_FILE_DIR, x) for x in projection_files];

    return args;


@test
def test_threshold():
    args = get_default_args();

    args["threshold"] = 10;

    return args;


@test
def test_all_sigs():
    args = get_default_args();

    args["all_sigs"] = True;

    return args;


@test
def test_debug():
    args = get_default_args();

    args["debug"] = True;

    return args;


@test
def test_lean():
    args = get_default_args();

    args["lean"] = True;

    return args;


@test
def test_sig_norm_method_none():
    args = get_default_args();

    args["sig_norm_method"] = "none";

    return args;


@test
def test_sig_norm_method_znorm_columns():
    args = get_default_args();

    args["sig_norm_method"] = "znorm_columns";

    return args;


@test
def test_sig_norm_method_znorm_rows():
    args = get_default_args();

    args["sig_norm_method"] = "znorm_rows";

    return args;


@test
def test_sig_norm_method_znorm_rows_then_columns():
    args = get_default_args();

    args["sig_norm_method"] = "znorm_rows_then_columns";

    return args;


@test
def test_sig_norm_method_rank_norm_columns():
    args = get_default_args();

    args["sig_norm_method"] = "rank_norm_columns";

    return args;


@test
def test_sig_score_method_naive():
    args = get_default_args();

    args["sig_score_method"] = "naive";

    return args;


@test
def test_sig_score_method_weighted_avg():
    args = get_default_args();

    args["sig_score_method"] = "weighted_avg";

    return args;


@test
def test_sig_score_method_imputed():
    args = get_default_args();

    args["sig_score_method"] = "imputed";

    return args;


@test
def test_sig_score_method_only_nonzero():
    args = get_default_args();

    args["sig_score_method"] = "only_nonzero";

    return args;


@test
def test_output():
    args = get_default_args();

    args["output"] = "outfolder";

    return args;


def run_all():

    test_simple();
    test_hk();
    test_text_sigs();
    test_text_sigsB();
    test_multi_sigs();
    test_no_precomputed();
    test_no_sigs();
    test_qc();
    test_nofilter();
    test_nomodel();
    test_pca_filter();
    test_weights();
    test_projections();
    test_many_projections();
    test_threshold();
    test_all_sigs();
    test_debug();
    test_lean();
    test_sig_norm_method_none();
    test_sig_norm_method_znorm_columns();
    test_sig_norm_method_znorm_rows();
    test_sig_norm_method_znorm_rows_then_columns();
    test_sig_norm_method_rank_norm_columns();
    test_sig_score_method_naive();
    test_sig_score_method_weighted_avg();
    test_sig_score_method_imputed();
    test_sig_score_method_only_nonzero();
    test_output();


