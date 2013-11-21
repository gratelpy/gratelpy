VERSION = (0,1,1)

from multiprocessing import Pool
from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.stoich import get_graph_stoich
from gratelpy.fragments import get_sensible_fragments
from gratelpy.subgraphs import get_all_valid_subgraphs
from functools import partial
import os

def get_version(*args, **kwargs):
    return str(VERSION[0])+'.'+str(VERSION[1])+'.'+str(VERSION[2])

def analyze(name, no_species, no_proc=1):

    pool = Pool(no_proc)

    alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, \
        dict_constants_reverse = get_network_from_mechanism(name, no_species)
			
    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
    valid_fragments = get_sensible_fragments(G, stoich_rank)

    get_subgraphs = partial(get_all_valid_subgraphs, G, stoich_rank)

    results = []
    r = pool.map_async(get_subgraphs, valid_fragments, callback=results.append)
    r.wait()

    return results[0]

def get_mechanism(name):
    _me = path = os.path.split(os.path.realpath(__file__))[0]
    return os.path.join(_me, 'mechanisms', name)
