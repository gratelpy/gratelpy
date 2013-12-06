import unittest
import sys

from test_reversible_substrate_inhibition import TestReversibleSubstrateInhibition

def runtests():
    suite = unittest.TestSuite()
    
    # add TestCase objects
    suite.addTest(TestReversibleSubstrateInhibition('test_alpha'))
    suite.addTest(TestReversibleSubstrateInhibition('test_beta'))
    suite.addTest(TestReversibleSubstrateInhibition('test_species'))
    suite.addTest(TestReversibleSubstrateInhibition('test_constants'))
    suite.addTest(TestReversibleSubstrateInhibition('test_fragments'))
    suite.addTest(TestReversibleSubstrateInhibition('test_subgraphs'))

    unittest.TextTestRunner().run(suite)
    return suite
