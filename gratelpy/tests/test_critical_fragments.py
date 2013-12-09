import unittest

try:
    from collections import Counter
except ImportError:
    from gratelpy.counter_python26 import Counter

from gratelpy import get_mechanism, analyze_one_proc

ks_index = -1
frag_i = 0
sg_i = 2

class TestCriticalFragments(unittest.TestCase):
    def get_critical_count(self, results):
        return sum([1 for r in results if r[ks_index] < 0.])
    def check_multiplicity(self, results):
        for result in results:
            species = Counter(result[frag_i][0])
            reactions = Counter(result[frag_i][1])
            subgraphs = result[sg_i]
            for subgraph in subgraphs:
                sg_s = Counter([el[0] for el in subgraph])
                sg_r = Counter([el[1] for el in subgraph])
                self.assertEqual(species, sg_s)
                self.assertEqual(reactions, sg_r)

    def test_reversible_substrate(self):
        mechanism = get_mechanism('reversible_substrate_inhibition.txt')
        no_spec = 4
        expected_critical = 1
        results = analyze_one_proc(mechanism, no_spec)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)

    def test_cdc42_yeast(self):
        mechanism = get_mechanism('cdc42_yeast.txt')
        no_spec = 8
        rank = 5
        expected_critical = 35
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)

    def test_glycolysis_gluconeogenesis_rank_2(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 2
        expected_critical = 1
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)

    def test_glycolysis_gluconeogenesis_rank_3(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 3
        expected_critical = 8
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)

    def test_glycolysis_gluconeogenesis_rank_4(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 4
        expected_critical = 12
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)

    def test_glycolysis_gluconeogenesis_rank_5(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 5
        expected_critical = 5
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)

    def test_single_layer_mapk(self):
        mechanism = get_mechanism('single_layer_mapk_mechanism.txt')
        no_spec = 9
        rank = 6
        expected_critical = 9
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)
