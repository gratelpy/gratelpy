import unittest
import sys

from test_reversible_substrate_inhibition import TestReversibleSubstrateInhibition
from test_stoich_ranks import TestStoichRanks
from test_critical_fragments import TestCriticalFragments

def runtests():
    suite = unittest.TestSuite()

    # add TestCase objects
    suite.addTest(TestReversibleSubstrateInhibition('test_stoich_rank'))
    suite.addTest(TestReversibleSubstrateInhibition('test_alpha'))
    suite.addTest(TestReversibleSubstrateInhibition('test_beta'))
    suite.addTest(TestReversibleSubstrateInhibition('test_species'))
    suite.addTest(TestReversibleSubstrateInhibition('test_constants'))
    suite.addTest(TestReversibleSubstrateInhibition('test_fragments'))
    suite.addTest(TestReversibleSubstrateInhibition('test_subgraphs'))

    suite.addTest(TestStoichRanks('test_reversible_substrate'))
    suite.addTest(TestStoichRanks('test_cdc42_yeast'))
    suite.addTest(TestStoichRanks('test_glycolysis_gluconeogenesis'))
    suite.addTest(TestStoichRanks('test_single_layer_mapk'))
    suite.addTest(TestStoichRanks('test_double_layer_mapk'))

    suite.addTest(TestCriticalFragments('test_reversible_substrate'))
    suite.addTest(TestCriticalFragments('test_cdc42_yeast'))
    suite.addTest(TestCriticalFragments('test_glycolysis_gluconeogenesis_rank_2'))
    suite.addTest(TestCriticalFragments('test_glycolysis_gluconeogenesis_rank_3'))
    suite.addTest(TestCriticalFragments('test_glycolysis_gluconeogenesis_rank_4'))
    suite.addTest(TestCriticalFragments('test_glycolysis_gluconeogenesis_rank_5'))
    suite.addTest(TestCriticalFragments('test_single_layer_mapk'))
    
    unittest.TextTestRunner().run(suite)
    return suite
