import unittest

import os
path = os.path.split(os.path.realpath(__file__))[0]

from multiprocessing import freeze_support
from multiprocessing.managers import SyncManager
import Queue
from gratelpy.gensg import make_server_manager
from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.stoich import get_graph_stoich
from gratelpy.fragments import get_sensible_fragments
from gratelpy.pyin import resource_path

def main():
    alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, \
        dict_constants_reverse = get_network_from_mechanism(
                resource_path(os.path.join(path, '..', 'mechanisms'),
                              'reversible_substrate_inhibition.txt'), 
                4)
			
    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
    valid_fragments = get_sensible_fragments(G, stoich_rank)

    # we will emulate the behavior of both gratelpy_fragment_server
    # and gratelpy_subclient here

    # server sets up fragment queue
    server_manager = make_server_manager(50000, 'smallg')
    server_fragment_q = server_manager.get_fragment_q()
    server_valid_frag_q = server_manager.get_valid_frag_q()
    for f in valid_fragments:
        print('Putting fragment %s in queue' % str(f))
        server_fragment_q.put(f)
    print('Fragment queue filled.')

    # client fetches all fragments one-by-one and processes them
    class QueueManager(SyncManager): pass    
    QueueManager.register('get_fragment_q')
    QueueManager.register('get_valid_frag_q')
    
    client_manager = QueueManager(address=('localhost', 50000), 
                                  authkey='smallg')
    print('Client connecting ...')
    client_manager.connect()
    print('Client connected.')

    client_fragment_q = client_manager.get_fragment_q()
    client_valid_frag_q = client_manager.get_valid_frag_q()

    print('Client starts fetching fragments off the queue ...')
    while True:
        try:
            f = client_fragment_q.get_nowait()
            print('Client fetched fragment %s.' % str(f))
            client_valid_subgraphs = get_all_valid_subgraphs(G, 
                                                             stoich_rank, 
                                                             f)
            client_valid_frag_q.put(client_valid_subgraphs)
        except Queue.Empty:
            break

    # now let the server go through the list of data deposited by the client
    while True:
        try:
            server_valid_subgraphs = server_valid_frag_q.get()
            print server_valid_subgraphs
        except:
            break

if __name__ == '__main__':
    # Required on Windows 
    # http://docs.python.org/2/library/multiprocessing.html#multiprocessing.freeze_support
    freeze_support()
    main()
