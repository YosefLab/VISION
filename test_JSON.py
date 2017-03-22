"""
Test conversions of FastProject datatypes to/from JSON
"""
import numpy as np
import unittest
from ..DataTypes import ExpressionData
from ..Signatures import Signature, SignatureScores
from .. import jsonIO

SAMPLE_EXPRESSION_DATA_JSON = """
{"col_labels": ["sample1", "sample2", "sample3"],
 "row_labels": ["feature1", "feature2", "feature3"],
 "data": [[1.2302, 3.43, 2.34],
            [3.23, 5.93, 1.03],
            [4.12, 9.23, 1.00]]
}
"""

SAMPLE_SIGNATURE_JSON = """
[{
    "sig_dict": {"gene1": 1, "gene2": -1, "gene3": 1},
    "signed": true,
    "source": "mysignaturefile.gmt",
    "name": "sample signature"
}]"""

SAMPLE_PRECOMPUTED_SIGNATURE_JSON = """
{"sample precomputed": {
        "scores": [3.42, 5.04, 5.00],
        "name": "sample precomputed",
        "sample_labels": ["sample1", "sample2", "sample3"],
        "isFactor": false,
        "isPrecomputed": true,
        "numGenes": 0
    }
}"""

SAMPLE_PRECOMPUTED_SIGNATURE2_JSON = """
{"sample precomputed": {
    "scores": ["level1", "level2", "level1"],
    "name": "sample precomputed",
    "sample_labels": ["sample1", "sample2", "sample3"],
    "isFactor": true,
    "isPrecomputed": true,
    "numGenes": 0
    }
}"""


class TestJsonIO(unittest.TestCase):

    def testExpressionData_roundTrip(self):
        """
        Test roundtrip conversion
        """

        data = np.random.randint(0, 100, size=(3, 3)) / 100.0;
        row_labels = ['a', 'b', 'c'];
        col_labels = ['d', 'e', 'f'];

        testData = ExpressionData(data,
                row_labels=row_labels,
                col_labels=col_labels);

        json_str = jsonIO.expressionData_to_JSON(testData);

        testData2 = jsonIO.JSON_to_ExpressionData(json_str);

        data2 = testData2.base;
        row_labels2 = testData2.row_labels;
        col_labels2 = testData2.col_labels;

        for rl1, rl2 in zip(row_labels, row_labels2):
            self.assertEqual(rl1, rl2)

        for cl1, cl2 in zip(col_labels, col_labels2):
            self.assertEqual(cl1, cl2)

        self.assertEqual(data.shape[0], data2.shape[0]);
        self.assertEqual(data.shape[1], data2.shape[1]);

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                self.assertEqual(data[i, j], data2[i, j]);

    def testJSON_to_Signature_roundTrip(self):
        """
        Test roundtrip conversion
        """

        sig_dict = {'gene1': 1, 'gene2': -1, 'gene3': 1};
        signed = False;
        source = 'myfile.gmt';
        name = 'test signature';

        testSig = Signature(sig_dict, signed, source, name);

        testSigs = [testSig];

        testJson = jsonIO.signatures_to_JSON(testSigs);

        resultSigs = jsonIO.JSON_to_Signatures(testJson);

        for a, b in zip(testSigs, resultSigs):
            self.assertEqual(a.signed, b.signed)
            self.assertEqual(a.source, b.source)
            self.assertEqual(a.name, b.name)

            self.assertEqual(a.sig_dict, b.sig_dict);

    def testJSON_to_SignatureScore_roundTrip(self):
        """
        Test roundtrip conversion
        """

        scores = [3.42, 10, 2.3]
        name = 'sigName'
        sample_labels = ['sample1', 'sample2', 'sample3']
        isFactor = False
        isPrecomputed = True
        numGenes = 4

        testSigScores = SignatureScores(scores, name, sample_labels,
                isFactor, isPrecomputed, numGenes);

        testSigScoresDict = {name: testSigScores};

        testJson = jsonIO.precomputed_signatures_to_JSON(testSigScoresDict);

        resultSigScoresDict = jsonIO.JSON_to_SignatureScores(testJson);

        for key in testSigScoresDict:
            orig = testSigScoresDict[key];
            result = resultSigScoresDict[key]
            self.assertEqual(orig.scores, result.scores);
            self.assertEqual(orig.name, result.name);
            self.assertEqual(orig.sample_labels, result.sample_labels);
            self.assertEqual(orig.isFactor, result.isFactor);
            self.assertEqual(orig.isPrecomputed, result.isPrecomputed);
            self.assertEqual(orig.numGenes, result.numGenes);

    def testJSON_to_ExpressionData(self):

        json_str = SAMPLE_EXPRESSION_DATA_JSON;

        result = jsonIO.JSON_to_ExpressionData(json_str);

        self.assertIsInstance(result, ExpressionData)

    def testJSON_to_Signature(self):

        json_str = SAMPLE_SIGNATURE_JSON;

        result = jsonIO.JSON_to_Signatures(json_str);

        self.assertIsInstance(result, list);
        self.assertIsInstance(result[0], Signature);

    def testJSON_to_SignatureScore(self):

        # Test first example str
        json_str = SAMPLE_PRECOMPUTED_SIGNATURE_JSON;

        result = jsonIO.JSON_to_SignatureScores(json_str);

        self.assertIsInstance(result, dict);
        self.assertIsInstance(result.values()[0], SignatureScores);

        # Test second example string (for a factor precomputed sig)
        json_str = SAMPLE_PRECOMPUTED_SIGNATURE2_JSON;

        result = jsonIO.JSON_to_SignatureScores(json_str);

        self.assertIsInstance(result, dict);
        self.assertIsInstance(result.values()[0], SignatureScores);
