#!/usr/bin/env python
"""benchmark.py

Usage: benchmark.py bechnmark_file your_file"""

import cPickle as pickle
import sys

# data format
frag_i = 0
sc_i = 1
sg_i = 2
ks_i = 3

# global variables available to interactive python sessions
benchmark_file_name = None
your_file_name = None
benchmark_fs = None
your_fs = None

b_fragments = None
y_fragments = None
b_critical = None
y_critical = None
b_critical_extras = None
y_critical_extras = None
b_sc = None
y_sc = None
b_sg = None
y_sg = None
b_ks = None
y_ks = None
b_dup = None

def data_summary(data, data_filename):
    print '----------------------'
    print data_filename
    print 'number of data sets:',str(len(data))
    print 'number of critical fragments:',str(len([s[ks_i] for s in data if s[ks_i]<0]))
    print '----------------------'

def duplicate_fragments(b_data, b_filename, y_data, y_filename):
    global b_fragments
    global y_fragments
    global b_dup

    b_duplicates = False
    y_duplicates = False

    if len(b_fragments) < len(b_data):
        print str(len(b_data)),'fragments in data set',b_filename,'-',str(len(b_fragments)),'unique fragments found'
        b_duplicates = True
    elif len(b_fragments) > len(b_data):
        raise
    
    if len(y_fragments) < len(y_data):
        print str(len(y_data)),'fragments in data set',y_filename,'-',str(len(y_fragments)),'unique fragments found'
        y_duplicates = True
    elif len(y_fragments) > len(y_data):
        raise

    if b_duplicates:
        b_dup = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): [] for s in b_data}
        for s in b_data:
            b_dup[(frozenset(s[frag_i][0]), frozenset(s[frag_i][1]))].append(s)
    
    
def data_compare(b_data, b_filename, y_data, y_filename):
    global b_fragments
    global y_fragments
    global b_critical
    global y_critical
    global b_critical_extras
    global y_critical_extras
    global b_sc
    global y_sc
    global b_sg
    global y_sg
    global b_ks
    global y_ks

    b_fragments = frozenset([(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])) for s in b_data])
    y_fragments = frozenset([(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])) for s in y_data])

    b_sc = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): s[sc_i] for s in b_data}
    y_sc = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): s[sc_i] for s in y_data}

    b_sg = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): frozenset([frozenset([path for path in sg]) for sg in s[sg_i]]) for s in b_data}
    y_sg = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): frozenset([frozenset([path for path in sg]) for sg in s[sg_i]]) for s in y_data}

    b_ks = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): s[ks_i] for s in b_data}
    y_ks = {(frozenset(s[frag_i][0]), frozenset(s[frag_i][1])): s[ks_i] for s in y_data}

    print 'set of fragments in',b_filename,'no. of fragments:',len(b_fragments)
    print 'set of fragments in',y_filename,'no. of fragments:',len(y_fragments)
    print 'size of intersection of these sets of fragments:',len(b_fragments.intersection(y_fragments))

    print 'fragments in',b_filename,'but not in',y_filename,':',len(b_fragments.difference(y_fragments))
    print 'fragments in',y_filename,'but not in',b_filename,':',len(y_fragments.difference(b_fragments))

    b_critical = frozenset([((frozenset(s[frag_i][0]), frozenset(s[frag_i][1])), frozenset([frozenset([path for path in sg]) for sg in s[sg_i]]), s[ks_i]) for s in b_data if s[ks_i]<0.0])

    y_critical = frozenset([((frozenset(s[frag_i][0]), frozenset(s[frag_i][1])), frozenset([frozenset([path for path in sg]) for sg in s[sg_i]]), s[ks_i]) for s in y_data if s[ks_i]<0.0])

    y_critical_extras = list(y_critical.difference(b_critical))
    b_critical_extras = list(b_critical.difference(y_critical))

    for yf in y_critical_extras:
        print 'sg in b:',len(b_sg[yf[0]]),', sg in y:',len(y_sg[yf[0]])
    
def main():
    # declare global variables
    global benchmark_file_name
    global your_file_name
    global benchmark_fs
    global your_fs

    try:
        benchmark_file_name = sys.argv[1]
        your_file_name = sys.argv[2]
    except IndexError:
        print __doc__
        sys.exit(2)

    try:
        print 'loading benchmark data from', benchmark_file_name 
        benchmark_fs = pickle.load(open(benchmark_file_name))
    except IOError, e:
        if e.errno == errno.ENOENT: # no such file or directory
            raise Exception('unable to open benchmark file \'',benchmark_file_name,'\'')
        else:
            raise Exception('unknown error')

    try:
        print 'loading your data from', your_file_name
        your_fs = pickle.load(open(your_file_name))
    except IOError, e:
        if e.errno == errno.ENOENT: # no such file or directory
            raise Exception('unable to open your file \'',your_file_name,'\'')
        else:
            raise Exception('unknown error')

    # data summary
    data_summary(benchmark_fs, benchmark_file_name)
    data_summary(your_fs, your_file_name)

    # compare data sets
    data_compare(benchmark_fs, benchmark_file_name, your_fs, your_file_name)
        
if __name__ == '__main__':
    main()

# old_cf_name = '/usr/users/cbu/waltherg/JIC/Projects/my_stuff/networkx/LPA/saved_states/double_layer_mapk_mechanism.fs'
# new_cf_name = 'single_layer_mapk_mechanism.vsg'

# #old_cf_name = '/usr/users/cbu/waltherg/JIC/Projects/my_stuff/networkx/LPA/saved_states/single_layer_mapk_mechanism.fs'
# #new_cf_name = 'single_layer_mapk_mechanism.vsg'

# old_cf = pickle.load(open(old_cf_name))
# old_cf = [(ocf, False) for ocf in old_cf] # mark each old critical fragment with a 'visited' switch
# #old_fragments = pickle.load(open(old_fragments_name))
# new_cf = pickle.load(open(new_cf_name))

# cf_matched_ctr = 0

# for ncf in new_cf:
#     in_old_cf = False
#     in_old_cf_more_than_once = False
#     for ocf_tuple in old_cf:
#         ocf = ocf_tuple[0]
#         ocf_visited = ocf_tuple[1]
#         if set(ncf[0][0]) == set(ocf[0][0]) and set(ncf[0][1]) == set(ocf[0][1]):
#             if ocf_visited == False:
#                 # substance and reaction node sets of new and old critical fragments agree
#                 # current old critical fragment hasn't been visited before, so current ocf is seen for first time as a match
#                 # let's now check if ncf and ocf also share exactly the same subgraphs
#                 sg_ocf_sets = [set(sg_ocf) for sg_ocf in ocf[-2]]
#                 if all([set(sg_ncf) in sg_ocf_sets for sg_ncf in ncf[-2]]):
#                     in_old_cf = True
#                 ocf_visited = True
#             else:
#                 in_old_cf_more_than_once = True

#     if (not in_old_cf) or in_old_cf_more_than_once:
#         print ncf
#         break
#     else:
#         cf_matched_ctr = cf_matched_ctr + 1
#         #print 'critical fragment matched:'+str(ncf[0])
