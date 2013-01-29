#!/usr/bin/env python
"""subclient.py

Usage: subclient.py mechanism_file num_complexes fragment_server_hostname
"""
import socket
hostname = socket.gethostname()

import sys
if hostname == 'pinguinzinho':
	path_to_networkx_dev = '/home/waltherg/Dropbox/GratelPy/networkx-dev'
elif hostname == 'pinguim':
	path_to_networkx_dev = '/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/networkx-dev/'
else:
	raise Exception('hostname unknown')

sys.path.insert(1, path_to_networkx_dev)
import networkx as nx

import Queue
from multiprocessing.managers import SyncManager

# import and use resource stuff to limit amount of memory each sublcient is allowed
import resource

# # DEBUG
# from pympler import tracker
import gc
# # DEBUG END

import numpy as np

from parse_mechanism import get_network_from_mechanism
from subgraphs import get_all_valid_subgraphs
from stoich import get_graph_stoich

# #DEBUG
# memory_tracker = None

def process_fragment_q(fragment_q, valid_frag_q, G, stoich_rank):
    while True:
        try:
            f = fragment_q.get_nowait()
            print 'Got', f
            valid_sgs = get_all_valid_subgraphs(G, stoich_rank, f)
            valid_frag_q.put(valid_sgs)
            #DEBUG
            #memory_tracker.print_diff()

        except Queue.Empty:
            print 'All done'
            return

def main():
    try:
        mechanism_file = sys.argv[1]
        num_complexes = int(sys.argv[2])
        frag_server = sys.argv[3]
    except IndexError:
        print __doc__
        sys.exit(2)

    class QueueManager(SyncManager): pass
    
    QueueManager.register('get_fragment_q')
    QueueManager.register('get_valid_frag_q')
    
    m = QueueManager(address=(frag_server, 50000), authkey='smallg')
    
    m.connect()
    
    fragment_q = m.get_fragment_q()
    valid_frag_q = m.get_valid_frag_q()
    
    alpha, beta, _, _, _, _ = get_network_from_mechanism(mechanism_file, num_complexes)
#    alpha = np.array(
#         [[1,0,0,0,0,0,0,0,0,0,0,0],
#          [1,0,0,1,0,0,0,0,0,0,0,0],
#          [0,1,1,0,0,0,0,0,0,0,0,0],
#          [0,0,0,1,0,0,0,0,0,1,0,0],
#          [0,0,0,0,1,1,0,0,0,0,0,0],
#          [0,0,0,0,0,0,1,0,0,0,0,0],
#          [0,0,0,0,0,0,1,0,0,1,0,0],
#          [0,0,0,0,0,0,0,1,1,0,0,0],
#          [0,0,0,0,0,0,0,0,0,0,1,1]]
#         )
#    
#    beta = np.array(
#         [[0,1,0,0,0,0,0,0,0,0,0,1],
#          [0,1,1,0,1,1,0,0,0,0,0,0],
#          [1,0,0,0,0,0,0,0,0,0,0,0],
#          [0,0,1,0,1,0,0,0,1,0,0,0],
#          [0,0,0,1,0,0,0,0,0,0,0,0],
#          [0,0,0,0,0,1,0,1,0,0,0,0],
#          [0,0,0,0,0,0,0,1,1,0,1,1],
#          [0,0,0,0,0,0,1,0,0,0,0,0],
#          [0,0,0,0,0,0,0,0,0,1,0,0]]
#         )
#

    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
#    stoich_rank = 6
    
    process_fragment_q(fragment_q, valid_frag_q, G, stoich_rank)

if __name__ == '__main__':
    #DEBUG
    #memory_tracker = tracker.SummaryTracker()
    
    # limit virtual memory (address space) to 2 GB
    # see man getrlimit for identifier RLIMIT_AS
    resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648)) 
    
    #DEBUG
    # if not gc.isenabled():
    #     gc.enable()
    # gc.set_debug(gc.DEBUG_LEAK) 
    main()
