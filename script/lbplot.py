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

fig, ax = plt.subplots(4, 1, figsize=(8.27, 11.69))

#ax[3].plot(max, c='C2', ls='-', label='max U>C')
#ax[3].plot(avg, c='r', ls='-', label='avg curve')

#I, W, tpro, cumpro, decpro, timepro, C, max, avg = get_params(fpro)
#I, W, tstatic, cumstatic, decstatic, timestatic, C, max, avg = get_params(fstatic)
#I, W, t100, cum100, dec100, time100, C, max, avg = get_params(ffreq100)

ax[0].set_title(sys.argv[1])

ax[0].plot(configs['static']['data']['W'], label='Workload')
#ax[0].plot(np.array(W) /100.0, label='AVG Workload')
ax[0].set_ylabel('Application Workload')
ax[0].legend()

for i, fbab in enumerate(babfiles):
    babcfg = get_params(fbab)
    ax[1].plot(babcfg['time'], label='Optimum')
    ax[2].plot(babcfg['cumli'], ls='-', label='Optimum')
    ax[3].plot(babcfg['max'], label='opt')
    for name, cfg in configs.items():
        print("Bab", i, "with", babcfg['time'][-1], "is", compute_performance_diff(babcfg['time'][-1], float(cfg['data']['time'][-1])), "% faster than", name, "with", float(cfg['data']['time'][-1]) )

    #print("Bab", i, "is", compute_performance_diff(timebab[-1], timemen1[-1]),"% faster than U>C -1")
    #print("Bab", i, "is", compute_performance_diff(timebab[-1], timepro[-1]), "% faster than procassini")
    #print("Bab", i, "is", compute_performance_diff(timebab[-1], time100[-1]), "% faster than rebalancing every 100 it")

#ax[3].legend()
#print(np.sum(np.asarray(cummen)) + np.sum(decmen) * C)
#common_sub_seqs = get_longest_common_lb_sequence(decbab, decmen)
#sa, sb = [], []
#ca, cb = [], []
#ta, tb = [], []
#
#for p, c in common_sub_seqs:
#    sa.append(decbab[p:c])
#    sb.append(decmen[p:c])
#    ca.append(cumbab[p:c])
#    cb.append(cummen[p:c])
#    ta.append(timebab[p:c])
#    tb.append(timemen[p:c])

#ax[1].plot(timepro, label='procassini')
for name, cfg in configs.items():
    if name != 'static':
        ax[1].plot(cfg['data']['time'], label=name)
        ax[2].plot(cfg['data']['cumli'], label=name)
        ax[3].plot(cfg['data']['max'], label=name)


#ax[2].plot(configs['menon']['data']['cumli'],label='menon')
#ax[2].plot(configs['menon1']['data']['cumli'],label='menon++')

#ax[1].plot(-compute_performance_diff(timebab, timemen1), label='U + $\Delta_{i+1}$ > C')

ax[1].set_ylabel('Relative performance')
ax[1].legend()

ax[2].plot([configs['static']['data']['C']]*configs['static']['data']['I'], label='C')

ax[2].set_xlabel('iteration')
ax[2].set_ylabel('Cumulative imbalance time')
ax[2].legend()
ax[3].legend()
plt.tight_layout()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()

"""
fig, ax = plt.subplots(len(sa), 1, sharex=True)
if len(sa) > 1:
    for seq_a, seq_b, li_a, li_b, time_a, time_b, s, i in (list(zip(sa, sb, ca, cb, ta,tb, common_sub_seqs, range(len(sa))))):
        it = range(*s)
        ax[i].plot(it, li_a)
        ax[i].plot(it, li_b)
        ax[i].plot(it, [C]*len(it), label='C')
else:
    for seq_a, seq_b, li_a, li_b, time_a, time_b, s, i in (list(zip(sa, sb, ca, cb, ta,tb, common_sub_seqs, range(len(sa))))):
        it = range(*s)
        ax.plot(it, li_a)
        ax.plot(it, li_b)
        ax.plot(it, [C]*len(it), label='C')
plt.title('cumulative LI')
plt.show()
plt.close()

fig, ax = plt.subplots(len(sa), 1, figsize=(10,50), sharex=True)
if len(sa) > 1:
    for seq_a, seq_b, li_a, li_b, time_a, time_b, s, i in (list(zip(sa, sb, ca, cb, ta,tb, common_sub_seqs, range(len(sa))))):
        it = range(*s)
        ax[i].plot(it, time_a, label='bab')
        ax[i].plot(it, time_b, label='bab')
        ax[i].legend()
else:
    for seq_a, seq_b, li_a, li_b, time_a, time_b, s, i in (list(zip(sa, sb, ca, cb, ta,tb, common_sub_seqs, range(len(sa))))):
        it = range(*s)
        ax.plot(it, time_a)
        ax.plot(it, time_b)

plt.title('time')
plt.show()
plt.close()
"""
