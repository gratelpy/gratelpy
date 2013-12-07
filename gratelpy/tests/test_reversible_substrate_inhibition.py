import unittest
import cPickle as pickle

import os
path = os.path.split(os.path.realpath(__file__))[0]

from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.stoich import get_graph_stoich
from gratelpy.fragments import get_sensible_fragments
from gratelpy.subgraphs import get_all_valid_subgraphs
from gratelpy.pyin import resource_path

alpha_reference = [[1, 0, 1, 0, 1, 0],
                   [0, 0, 1, 0, 0, 0],
                   [0, 0, 0, 1, 1, 0],
                   [0, 0, 0, 0, 0, 1]]

beta_reference = [[0, 1, 0, 0, 0, 1],
                  [0, 0, 0, 1, 0, 0],
                  [0, 0, 1, 0, 0, 1],
                  [0, 0, 0, 0, 1, 0]]

species_reference = {'[A1]': 0, '[A2]': 1, '[A3]': 2, '[A4]': 3}
constants_reference = {'k1': 0, 'k2': 1, 'k3': 2, 'k4': 3, 'k5': 4, 'k6': 5}

fragments_reference = [(('s3', 's2', 's1'), ('w5', 'w3', 'w5')), 
                       (('s3', 's2', 's1'), ('w5', 'w3', 'w3')), 
                       (('s3', 's2', 's1'), ('w5', 'w3', 'w1')), \
                           (('s3', 's2', 's1'), ('w4', 'w3', 'w5')), 
                       (('s3', 's2', 's1'), ('w4', 'w3', 'w3')), \
                           (('s3', 's2', 's1'), ('w4', 'w3', 'w1')), 
                       (('s3', 's2', 's4'), ('w5', 'w3', 'w6')), \
                           (('s3', 's2', 's4'), ('w4', 'w3', 'w6')), 
                       (('s3', 's1', 's4'), ('w5', 'w5', 'w6')), \
                           (('s3', 's1', 's4'), ('w5', 'w3', 'w6')), 
                       (('s3', 's1', 's4'), ('w5', 'w1', 'w6')), \
                           (('s3', 's1', 's4'), ('w4', 'w5', 'w6')), 
                       (('s3', 's1', 's4'), ('w4', 'w3', 'w6')), \
                           (('s3', 's1', 's4'), ('w4', 'w1', 'w6')), 
                       (('s2', 's1', 's4'), ('w3', 'w5', 'w6')), \
                           (('s2', 's1', 's4'), ('w3', 'w3', 'w6')), 
                       (('s2', 's1', 's4'), ('w3', 'w1', 'w6'))]

subgraphs_reference = pickle.load(open(
        resource_path(path, 'test_reversible_substrate_inhibition.dat'), 'r'))
sg_fragments_reference = [sg[0] for sg in subgraphs_reference]
frag_i = 0
sc_i = 1
sg_i = 2
ks_i = 3

subgraphs_reference_dict = {}
for sg in subgraphs_reference:
    subgraphs_reference_dict[sg[frag_i]] = {'sc': sg[sc_i], 
                                           'sg': sg[sg_i], 
                                           'ks': sg[ks_i]}
sg_dict = subgraphs_reference_dict

class TestReversibleSubstrateInhibition(unittest.TestCase):

    _alpha, _beta, _species, _constants, _, _ =\
        get_network_from_mechanism(
        resource_path(os.path.join(path, '..', 'mechanisms'),
                      'reversible_substrate_inhibition.txt'), 
        4)
        
    _G, _stoich, _stoich_rank = \
        get_graph_stoich(_alpha, _beta)

    _vf = get_sensible_fragments(_G, _stoich_rank)
    
    _sg = []
    for frag in _vf:
        _sg.append(get_all_valid_subgraphs(_G, _stoich_rank, frag))
        
    def test_stoich_rank(self):
        self.assertEqual(self._stoich_rank, 3)

    def test_alpha(self):
        for row, row_ref in zip(self._alpha, alpha_reference):
            for el, el_ref in zip(row, row_ref):
                self.assertEqual(el, el_ref)

    def test_beta(self):
        for row, row_ref in zip(self._beta, beta_reference):
            for el, el_ref in zip(row, row_ref):
                self.assertEqual(el, el_ref)

    def test_species(self):
        self.assertEqual(self._species, species_reference)

    def test_constants(self):
        self.assertEqual(self._constants, constants_reference)

    def test_fragments(self):
        self.assertEqual(len(self._vf), len(fragments_reference))
        for frag in self._vf:
            self.assertTrue(frag in fragments_reference)

    def test_subgraphs(self):
        self.assertEqual(len(self._sg), len(subgraphs_reference))
        for sg in self._sg:
            self.assertEqual(sg[sc_i], sg_dict[sg[frag_i]]['sc'])
            self.assertEqual(sg[sg_i], sg_dict[sg[frag_i]]['sg'])
            self.assertEqual(sg[ks_i], sg_dict[sg[frag_i]]['ks'])
            
if __name__ == '__main__':
    unittest.main()
