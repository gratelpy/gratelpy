import unittest

try:
    from collections import Counter
except ImportError:
    from gratelpy.counter_python26 import Counter

from gratelpy import get_mechanism, analyze_one_proc
from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.stoich import get_graph_stoich
from gratelpy.utils import (result_get_fragment,
                            result_get_sc,
                            result_get_sg,
                            fragment_get_species,
                            fragment_get_reactions,
                            subgraph_get_species,
                            subgraph_get_reactions,
                            species_get_index,
                            reaction_get_index)
from gratelpy.subgraphs import get_subgraph_motifs

ks_index = -1
frag_i = 0
sc_i = 1
sg_i = 2

spec_i = 0
rxn_i = 1


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

    def check_duplicate_fragments(self, results):
        fragments = []
        for result in results:
            species = tuple(result[frag_i][0])
            reactions = tuple(result[frag_i][1])
            fragments.append((species, reactions))
        fragments = Counter(fragments)
        for f in fragments:
            self.assertEqual(fragments[f], 1)

    def check_fragment_notation(self, results, name, no_species, des_order):
        alpha, beta, _, _, _, _ = get_network_from_mechanism(name,
                                                             no_species)
        G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

        for result in results:
            sc = result_get_sc(result)
            fragment = result_get_fragment(result)
            species = fragment_get_species(fragment)
            self.assertEqual(len(species), des_order)
            reactions = fragment_get_reactions(fragment)

            fragment_e = {}  # edges encoded in fragment syntax
            for s_i, s in enumerate(species):
                fragment_e[s] = (s, reactions[s_i])

            # edges stored in subgraph components dictionary
            for s in species:
                # exactly one edge per species
                self.assertEqual(len(sc[s]['edges']), 1)
                # syntax-encoded edge equals sc-stored edge
                self.assertEqual(sc[s]['edges'][0], fragment_e[s])
                # edges saved in networkx.DiGraph object 'G'
                self.assertTrue(sc[s]['edges'][0] in G.edges())
                # edge (s1, r2) means alpha[0][1] > 0
                s_i = species_get_index(s)
                reaction = sc[s]['edges'][0][1]
                r_i = reaction_get_index(reaction)
                self.assertTrue(alpha[s_i][r_i] > 0)

    def check_edges_only_subgraphs(self, results):
        for result in results:
            sg = result_get_sg(result)
            self.assertTrue(sum(all([len(el) == 2 for el in a_sg])
                                for a_sg in sg) == 1)

    def check_sg_fragment_conform(self, results):
        for result in results:
            fragment = result_get_fragment(result)
            species = Counter(fragment_get_species(fragment))
            reactions = Counter(fragment_get_reactions(fragment))
            for sg in result_get_sg(result):
                self.assertEqual(species, Counter(subgraph_get_species(sg)))
                self.assertEqual(reactions,
                                 Counter(subgraph_get_reactions(sg)))

    def check_sg_sgmotif_conform(self, results):
        for result in results:
            sc = result_get_sc(result)
            sg_motifs = get_subgraph_motifs(sc)
            sgm_keys = [set(key) for key in sg_motifs.keys()]
            subgraphs = result_get_sg(result)
            self.assertEqual(len(sgm_keys), len(subgraphs))
            for sg in result_get_sg(result):
                sg_motif_match = [set(sg) == key for key in sgm_keys]
                self.assertEqual(sum(sg_motif_match), 1)

    def run_all(self, results,
                no_critical, expected_critical,
                mechanism, no_spec, rank):
        self.assertEqual(no_critical, expected_critical)
        self.check_multiplicity(results)
        self.check_duplicate_fragments(results)
        self.check_fragment_notation(results, mechanism, no_spec, rank)
        self.check_edges_only_subgraphs(results)
        self.check_sg_fragment_conform(results)
        self.check_sg_sgmotif_conform(results)

    def test_reversible_substrate(self):
        mechanism = get_mechanism('reversible_substrate_inhibition.txt')
        no_spec = 4
        rank = 3
        expected_critical = 1
        results = analyze_one_proc(mechanism, no_spec)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)

    def test_cdc42_yeast(self):
        mechanism = get_mechanism('cdc42_yeast.txt')
        no_spec = 8
        rank = 5
        expected_critical = 35
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)

    def test_glycolysis_gluconeogenesis_rank_2(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 2
        expected_critical = 1
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)

    def test_glycolysis_gluconeogenesis_rank_3(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 3
        expected_critical = 8
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)

    def test_glycolysis_gluconeogenesis_rank_4(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 4
        expected_critical = 12
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)

    def test_glycolysis_gluconeogenesis_rank_5(self):
        mechanism = get_mechanism('glycolysis_mechanism.txt')
        no_spec = 7
        rank = 5
        expected_critical = 5
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)

    def test_single_layer_mapk(self):
        mechanism = get_mechanism('single_layer_mapk_mechanism.txt')
        no_spec = 9
        rank = 6
        expected_critical = 9
        results = analyze_one_proc(mechanism, no_spec, rank)
        no_critical = self.get_critical_count(results)
        self.run_all(results, no_critical, expected_critical,
                     mechanism, no_spec, rank)
