#!/usr/bin/env python
"""subclient.py

Usage: subclient.py mechanism_file num_complexes
"""

import sys
import Queue
from multiprocessing.managers import SyncManager

# import and use resource stuff to limit amount of memory each sublcient is allowed
import resource

# # DEBUG
# from pympler import tracker
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
    except IndexError:
        print __doc__
        sys.exit(2)

    class QueueManager(SyncManager): pass
    
    QueueManager.register('get_fragment_q')
    QueueManager.register('get_valid_frag_q')
    
    m = QueueManager(address=('pinguinzinho', 50000), authkey='smallg')
    
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
    resource.setrlimit(resource.RLIMIT_DATA,(2147483648,2147483648)) # limit heap size to 2 GB
    main()
