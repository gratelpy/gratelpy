import cPickle as pickle
import sys

def main():
    try:
        fs_file_1 = sys.argv[1]
        fs_file_2 = sys.argv[2]
    except IndexError, e:
        print __doc__
        sys.exit(2)

    fs_1 = pickle.load(open(fs_file_1))
    fs_2 = pickle.load(open(fs_file_2))

    crit_1 = set()
    crit_2 = set()
    
    for s in fs_1:
        if s[-1]<0:
            crit_1.add(s[0])
    for s in fs_2:
        if s[-1]<0:
            crit_2.add(s[0])

    if crit_1 == crit_2:
        print 'sets are the same'

if __name__ == '__main__':
    main()
