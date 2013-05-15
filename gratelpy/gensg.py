import sys
import time
import Queue
import functools
from multiprocessing import Pool
from multiprocessing.managers import SyncManager

from subgraphs import get_all_valid_subgraphs

def make_server_manager(port, authkey):
    fragment_q = Queue.Queue()
    valid_frag_q = Queue.Queue()

    class JobQueueManager(SyncManager): pass

    JobQueueManager.register('get_fragment_q', callable=lambda: fragment_q)
    JobQueueManager.register('get_valid_frag_q', callable=lambda: valid_frag_q)

    manager = JobQueueManager(address=('', port), authkey=authkey)

    manager.start()

    return manager


def gen_valid_subgraphs(G, valid_fragments, stoich_rank):
    all_subgraphs = [get_all_valid_subgraphs(G, stoich_rank, f) for f in valid_fragments]
    valid_subgraphs = [sg for sg in all_subgraphs if sg is not None]

    return valid_subgraphs

def gen_valid_subgraphs_mp(G, valid_fragments, stoich_rank):
    pool = Pool()
    chunksize = 10

    get_val_sub = functools.partial(get_all_valid_subgraphs, G, stoich_rank)
    all_subgraphs = pool.imap(get_val_sub, valid_fragments, chunksize)
  
    #all_subgraphs = [get_val_sub(f) for f in valid_fragments]
    valid_subgraphs = [sg for sg in all_subgraphs if sg is not None]

    return valid_subgraphs

def gen_valid_subgraphs_mps(G, valid_fragments, stoich_rank):
    print 'gensg: gen_valid_subgraphs_mps received',str(len(valid_fragments)),'fragments. the first up to 20 fragments are:'
    for ctr in range(min(len(valid_fragments),20)):
        print valid_fragments[ctr]

    server_manager = make_server_manager(50000, 'smallg')
    
    fragment_q = server_manager.get_fragment_q()
    valid_frag_q = server_manager.get_valid_frag_q()
    for f in valid_fragments:
        fragment_q.put(f)

    valid_subgraphs = []

    n_results = 0
    target_results = len(valid_fragments)
    start = time.time()
    last_t = start
    intervals = [1] * 10
    while n_results < target_results:
        vsg = valid_frag_q.get()
        len_vsg = 0
        if vsg is not None:
            len_vsg = len(vsg)
        if vsg is not None:
            # here we expect that vsg includes the score K_S of the corresponding fragment at index -1
            if type(vsg[-1]) == type(float()):
                if vsg[-1] < 0:
                    valid_subgraphs.append(vsg) # only collect critical fragments
            else:
                raise

        n_results += 1
        perc_done = int(float(100 * n_results) / target_results)
        interval = time.time() - last_t
        last_t = time.time()
        intervals.append(interval)
        del(intervals[0])
        iavg = sum(intervals) / len(intervals)
        freq = 1 / iavg
        est_time = iavg * (target_results - n_results)
        print '%d (%d subgraphs) / %d (%d%%), %f, %ds' % (n_results, len_vsg, target_results, 
            perc_done, freq, est_time)

    return valid_subgraphs
