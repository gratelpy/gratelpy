"""
GraTeLPy
"""

__version__ = '0.2.0.1'

import os
import itertools

def get_version(*args, **kwargs):
    return __version__

def analyze_one_proc(name, no_species, rank=None):
    from gratelpy.parse_mechanism import get_network_from_mechanism
    from gratelpy.stoich import get_graph_stoich
    from gratelpy.fragments import get_sensible_fragments
    from gratelpy.subgraphs import get_all_valid_subgraphs

    alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, \
        dict_constants_reverse = get_network_from_mechanism(name, no_species)
			
    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
    if rank:
        stoich_rank = rank

    valid_fragments = get_sensible_fragments(G, stoich_rank)

    results = []
    for f in valid_fragments:
        results += [get_all_valid_subgraphs(G, stoich_rank, f)]
    
    return results

def analyze(name, no_species, no_proc=1):

    from multiprocessing import Pool
    from gratelpy.parse_mechanism import get_network_from_mechanism
    from gratelpy.stoich import get_graph_stoich
    from gratelpy.fragments import get_sensible_fragments
    from functools import partial

    pool = Pool(no_proc)

    alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, \
        dict_constants_reverse = get_network_from_mechanism(name, no_species)
			
    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
    valid_fragments = get_sensible_fragments(G, stoich_rank)

    

    results = []
    for res in pool.imap(get_subgraphs,
                         itertools.izip(itertools.repeat(G),
                                        itertools.repeat(stoich_rank),
                                        valid_fragments),
                         chunksize=round(float(len(valid_fragments))/float(no_proc))):
        results += [res]

    # The above looks a bit esoteric because it is a workaround of a bug
    # that existed in Python 2.6 and that was fixed only in Python 2.7.
    # To support 2.6 I had to use the above workaround -- including
    # making 'get_subgraphs' a proper method of this package (see below).
    # http://techguyinmidtown.com/2009/01/23/hack-for-functoolspartial-and-multiprocessing/
    # http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
    #
    # In Python 2.7 the following works:
    #get_subgraphs = partial(get_all_valid_subgraphs, G, stoich_rank)
    #r = pool.map_async(get_subgraphs, valid_fragments, callback=results.append)
    #r.wait()

    return results[0]

def get_mechanism(name):
    _me = path = os.path.split(os.path.realpath(__file__))[0]
    return os.path.join(_me, 'mechanisms', name)

def get_subgraphs(arg):
    from gratelpy.subgraphs import get_all_valid_subgraphs
    return get_all_valid_subgraphs(arg[0], arg[1], arg[2])
