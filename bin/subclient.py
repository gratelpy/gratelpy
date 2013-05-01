#!/usr/bin/env python
"""subclient.py

Usage: subclient.py mechanism_file num_complexes [fragment_server_hostname=localhost] [port number on fragments server = 50000]
"""
import sys
import networkx as nx

import Queue
from multiprocessing.managers import SyncManager

# resource package to limit memory used by subclient
import resource

import numpy as np

from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.subgraphs import get_all_valid_subgraphs
from gratelpy.stoich import get_graph_stoich

def process_fragment_q(fragment_q, valid_frag_q, G, stoich_rank):
    while True:
        try:
            f = fragment_q.get_nowait()
            print 'Got', f
            valid_sgs = get_all_valid_subgraphs(G, stoich_rank, f)
            valid_frag_q.put(valid_sgs)

        except Queue.Empty:
            print 'All done'
            return

def main():
    try:
        mechanism_file = sys.argv[1]
        num_complexes = int(sys.argv[2])
    except IndexError:
        print __doc__
        sys.exit(2)

    try:
        frag_server = sys.argv[3]
    except IndexError:
        frag_server = 'localhost'

    try:
	    fragments_server_port = int(sys.argv[4])
    except IndexError:
	    fragments_server_port = 50000
	    

    class QueueManager(SyncManager): pass
    
    QueueManager.register('get_fragment_q')
    QueueManager.register('get_valid_frag_q')
    
    m = QueueManager(address=(frag_server, fragments_server_port), authkey='smallg')
    
    m.connect()
    
    fragment_q = m.get_fragment_q()
    valid_frag_q = m.get_valid_frag_q()
    
    alpha, beta, _, _, _, _ = get_network_from_mechanism(mechanism_file, num_complexes)

    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
    
    process_fragment_q(fragment_q, valid_frag_q, G, stoich_rank)

if __name__ == '__main__':
    # limit virtual memory (address space) to 2 GB
    # see man getrlimit for identifier RLIMIT_AS
    resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648)) 
    
    main()
