import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re

def get_longest_common_lb_sequence(seq_a, seq_b):
    seq_length = len(seq_a)
    assert (len(seq_a) == len(seq_b))
    common_lb = seq_a & seq_b
    #common_lb = np.insert(common_lb, 0, 1)
    common_lb[0]  = 1 # np.concatenate([[1], common_lb])
    common_lb[-1] = 1

    common_lb_iteration = np.argwhere(common_lb == 1).flatten()

    return np.array(list(zip(common_lb_iteration, common_lb_iteration[1:])))

def tryint(s):
    try:
        return int(s)
    except ValueError:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

def compute_performance_diff(a, b):
    return (b / a - 1.0) * 100.0

def get_params(fname):
    with open(fname) as f:
        I, W, t, cum, dec, tim, C, max, avg = f.readlines()
        W =   np.asarray([float(x) for x in   W.strip().split(' ')])
        cum2= np.asarray([float(x) for x in cum.strip().split(' ')])
        dec2= np.asarray([int(x)   for x in dec.strip().split(' ')])
        time2=np.asarray([float(x) for x in tim.strip().split(' ')])
        max=  np.asarray([float(x) for x in max.strip().split(' ')])
        avg=  np.asarray([float(x) for x in avg.strip().split(' ')])

    return {'I': int(I), 'W': W, 't': float(t), 'cumli': cum2, 'decision': dec2, 'time': time2, 'C': float(C), 'max': max, 'avg': avg}

dir    = "results/"
fmenon = dir + "menon-solution.txt"
fmenon1 = dir + "menon-solution-1.txt"

babfiles  = list( [dir + fname for fname in os.listdir(dir) if "optimal-solution" in fname] )
#filter( lambda f: "optimal-solution" in f, os.listdir('results/') ) )

sort_nicely(babfiles)

configs = {
    'menon': {'fname': "menon-solution.txt", 'data': []},
    'bastien': {'fname': "bastien-solution.txt", 'data': []},
    #'procassini': {'fname': "proca-solution.txt", 'data': []},
    'static': {'fname': "static-solution.txt", 'data': []},
    #'freq100': {'fname': "freq-100-solution.txt", 'data': []},
    #'freq50': {'fname': "freq-50-solution.txt", 'data': []},
    #'freq25': {'fname': "freq-25-solution.txt", 'data': []},
}

for k, cfg in configs.items():
    configs[k]['data'] = get_params(dir + cfg['fname'])

fig, ax = plt.subplots(2, 1)

#ax[3].plot(max, c='C2', ls='-', label='max U>C')
#ax[3].plot(avg, c='r', ls='-', label='avg curve')

#I, W, tpro, cumpro, decpro, timepro, C, max, avg = get_params(fpro)
#I, W, tstatic, cumstatic, decstatic, timestatic, C, max, avg = get_params(fstatic)
#I, W, t100, cum100, dec100, time100, C, max, avg = get_params(ffreq100)

ax[0].set_title(sys.argv[1])

#ax[0].plot(np.array(W) /100.0, label='AVG Workload')


for i, fbab in enumerate(babfiles):

    babcfg = get_params(fbab)
    #ax[0].plot(babcfg['max']-babcfg['avg'], label='Optimum')
    ax[0].plot(babcfg['time'], label='Optimum')
    ax[1].plot(babcfg['cumli'], ls='-', label='Optimum')
    #ax[3].plot(babcfg['max'], label='opt')
    for name, cfg in configs.items():
        print("Bab", i, "with", babcfg['time'][-1], "is", compute_performance_diff(babcfg['time'][-1], float(cfg['data']['time'][-1])), "% faster than", name, "with", float(cfg['data']['time'][-1]) )

for name, cfg in configs.items():
    if name != 'static':
        #ax[0].plot(cfg['data']['max'] - cfg['data']['avg'], label=name)
        ax[0].plot(cfg['data']['time'], label=name)
        ax[1].plot(cfg['data']['cumli'], label=name)

for canvas in ax:
    canvas.legend()

#ax[0].set_ylabel('Imbalance Time')
ax[0].set_ylabel('Simulated Parallel Time')

ax[1].plot([configs['static']['data']['C']]*configs['static']['data']['I'], label='C')

ax[1].set_xlabel('Iteration')
ax[1].set_ylabel('Cumulative Imbalance Time')

plt.tight_layout()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2], dpi=300)
else:
    plt.show()

