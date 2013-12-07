import unittest

from gratelpy import get_mechanism
from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.stoich import get_graph_stoich

class TestStoichRanks(unittest.TestCase):
    def check_rank(self, mechanism, no_spec, expected_rank):
        alpha, beta = get_network_from_mechanism(mechanism, no_spec)[:2]
        stoich_rank = get_graph_stoich(alpha, beta)[2]
        self.assertEqual(stoich_rank, expected_rank)

    def test_reversible_substrate(self):
        mechanism = get_mechanism('reversible_substrate_inhibition.txt')
        no_spec = 4
        expected_rank = 3
        self.check_rank(mechanism, no_spec, expected_rank)

    def test_cdc42_yeast(self):
        mechanism = get_mechanism('cdc42_yeast.txt')
        no_spec = 8
        expected_rank = 5
        self.check_rank(mechanism, no_spec, expected_rank)

    def test_glycolysis_gluconeogenesis(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        expected_rank = 5
        self.check_rank(mechanism, no_spec, expected_rank)

    def test_single_layer_mapk(self):
        mechanism = get_mechanism('single_layer_mapk_mechanism.txt')
        no_spec = 9
        expected_rank = 6
        self.check_rank(mechanism, no_spec, expected_rank)

    def test_double_layer_mapk(self):
        mechanism = get_mechanism('double_layer_mapk_mechanism.txt')
        no_spec = 12
        expected_rank = 9
        self.check_rank(mechanism, no_spec, expected_rank)
