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

import resource
# limit virtual memory (address space) to 2 GB
# see man getrlimit for identifier RLIMIT_AS
resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648)) 
#resource.setrlimit(resource.RLIMIT_AS,(4294967296,4294967296)) 

from parse_mechanism import get_network_from_mechanism
from stoich import get_graph_stoich
from subgraphs import get_subgraph_components, get_sensible_subgraphs
from fragments import get_sensible_fragments

import cPickle as pickle
import errno

mechanism_file = '../mechanisms/double_layer_mapk_mechanism.txt'
num_complexes = 12
alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism(mechanism_file, num_complexes)

G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

frag = (('s9', 's8', 's2', 's7', 's6', 's5', 's4', 's12', 's11'), ('w16', 'w10', 'w7', 'w8', 'w13', 'w6', 'w4', 'w17', 'w15'))

fn_vfrags = '../double_layer_mapk.vfrags'
try:
        valid_fragments = pickle.load(open(fn_vfrags))
except IOError, e:
    if e.errno == errno.ENOENT:
        valid_fragments = get_sensible_fragments(G, stoich_rank)
        pickle.dump(valid_fragments, open(fn_vfrags, 'wb'))
    else: raise

sc = get_subgraph_components(G, frag)
sensible_sgs = get_sensible_subgraphs(sc)

sg_ctr = 0
print 'sensible subgraphs: '
for sg in sensible_sgs:
	print sg
	sg_ctr = sg_ctr + 1

print sg_ctr
