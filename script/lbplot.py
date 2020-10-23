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

    return int(I), W, float(t), cum2, dec2, time2, float(C), max, avg

dir    = "results/"
fmenon = dir + "menon-solution.txt"
fmenon1 = dir + "menon-solution-1.txt"

babfiles  = list( [dir + fname for fname in os.listdir(dir) if "optimal-solution" in fname] )
#filter( lambda f: "optimal-solution" in f, os.listdir('results/') ) )

sort_nicely(babfiles)

fpro   = dir + "proca-solution.txt"
ffreq100=dir + "freq-100-solution.txt"

fig, ax = plt.subplots(4,1, figsize=(8.27, 11.69))

I, W, tmen, cummen, decmen, timemen, C, max, avg = get_params(fmenon)
I, W, tmen1, cummen1, decmen1, timemen1, C, max, avg = get_params(fmenon1)
ax[3].plot(max, c='C2', ls='-', label='max U>C')
ax[3].plot(avg, c='r', ls='-', label='avg curve')
I, W, tpro, cumpro, decpro, timepro, C, max, avg = get_params(fpro)
I, W, t100, cum100, dec100, time100, C, max, avg = get_params(ffreq100)

# if len(sys.argv) > 2:
#     data = np.asarray([len(x[1:]) for x in np.split(decmen, np.where(decmen == 1)[0].astype(int))])
#     bins = np.arange(0, data.max() + 1.5) - 0.5
#     plt.hist(data, bins=bins)
#     plt.title('Distribution of distance between two load balancings')
#     plt.savefig('men-histogram-tau-'+sys.argv[2])
#     plt.close()
#
# if len(sys.argv) > 2:
#     I, W, tbab, cumbab, decbab, timebab, C = get_params(babfiles[0])
#     data = np.asarray([len(x[1:]) for x in np.split(decbab, np.where(decbab == 1)[0].astype(int))])
#     bins = np.arange(0, data.max() + 1.5) - 0.5
#     plt.hist(data, bins=bins)
#     plt.title('Distribution of distance between two load balancings')
#     plt.savefig('bab-histogram-tau-'+sys.argv[2])
#     plt.close()

ax[0].set_title(sys.argv[1])

ax[0].plot(W, label='Workload')
#ax[0].plot(np.array(W) /100.0, label='AVG Workload')
ax[0].set_ylabel('Application Workload')
ax[0].legend()

for i, fbab in enumerate(babfiles):
    I, W, tbab, cumbab, decbab, timebab, C, max, avg = get_params(fbab)
    ax[1].plot(-compute_performance_diff(timebab, timebab), label='Branch and Bound %d' % i)
    ax[2].plot(cumbab, ls='-', label='Branch and Bound %d' % i)
    ax[3].plot(max, c='C0', ls='-', label='BaB max %d' % i)
    #ax[3].plot(avg, c='r', ls='-', label='avg curve')

    print("Bab", i, "is", compute_performance_diff(timebab[-1], timemen[-1]), "% faster than U>C")
    print("Bab", i, "is", compute_performance_diff(timebab[-1], timemen1[-1]),"% faster than U>C -1")
    print("Bab", i, "is", compute_performance_diff(timebab[-1], timepro[-1]), "% faster than procassini")
    print("Bab", i, "is", compute_performance_diff(timebab[-1], time100[-1]), "% faster than rebalancing every 100 it")

    # print("  Absolute difference is: ", (timemen[-1] - timebab[-1]), "(Bab =", timebab[-1], ")")

    # print(np.sum(np.asarray(cumbab)) + np.sum(decbab) * C)
ax[3].legend()
print(np.sum(np.asarray(cummen)) + np.sum(decmen) * C)
common_sub_seqs = get_longest_common_lb_sequence(decbab, decmen)
sa, sb = [], []
ca, cb = [], []
ta, tb = [], []

for p, c in common_sub_seqs:
    sa.append(decbab[p:c])
    sb.append(decmen[p:c])
    ca.append(cumbab[p:c])
    cb.append(cummen[p:c])
    ta.append(timebab[p:c])
    tb.append(timemen[p:c])

#ax[1].plot(timepro, label='procassini')
ax[1].plot(-compute_performance_diff(timebab, timemen), label='U > C')
ax[1].plot(-compute_performance_diff(timebab, timemen1), label='U + $\Delta_{i+1}$ > C')

ax[1].set_ylabel('Relative performance')
ax[1].legend()

#ax[2].plot(cumpro,  ls='dotted', label='procassini ($t\' < t$)' )
ax[2].plot(cummen,  ls='dotted',  label='U>C')
ax[2].plot(cummen1, ls='dashdot', label='U + $\Delta_{i+1}$ > C')

ax[2].plot([C]*I, label='C')

ax[2].set_xlabel('iteration')
ax[2].set_ylabel('Cumulative imbalance time')
ax[2].legend()
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
