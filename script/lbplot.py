import matplotlib.pyplot as plt
import numpy as np
import sys
import os

from matplotlib import rc
rc('text', usetex=True)


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
    import re
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


def build_figures(title, optimums, configs):
    fig, ax = plt.subplots(2, 1)
    ax[0].set_title(title)
    # plot optimum's data
    for i, fbab in enumerate(optimums):
        babcfg = get_params(fbab)
        ax[0].plot(babcfg['time'], label='Optimum')
        ax[1].plot(babcfg['cumli'], ls='-', label='Optimum')
        for name, cfg in configs.items():
            print("Bab", i, "with", babcfg['time'][-1], "is", compute_performance_diff(babcfg['time'][-1], float(cfg['data']['time'][-1])), "% faster than", name, "with", float(cfg['data']['time'][-1]) )
    # plot the criteria' data
    for name, cfg in configs.items():
        if name != 'static':
            ax[0].plot(cfg['data']['time'], label=name)
            ax[1].plot(cfg['data']['cumli'], label=name)
    # plot C line
    ax[1].plot([configs['static']['data']['C']]*configs['static']['data']['I'], label='C')
    # put legends in all plots
    for canvas in ax:
        canvas.legend()
    # write labels
    ax[0].set_ylabel('Simulated Parallel Time')
    ax[1].set_xlabel('Iteration')
    ax[1].set_ylabel('Cumulative Imbalance Time')
    return fig, ax


folder  = os.getcwd().split('/')[-1]
if 'LBOPT' in folder:
    dir = 'results/'
elif 'script' in folder:
    dir = '../results/'
else:
    raise Exception('unknown folder')

babfiles  = list( [os.path.join(dir, fname) for fname in os.listdir(dir) if "optimal-solution" in fname] )

sort_nicely(babfiles)

configs = {
    'menon': {'fname': "menon-solution.txt", 'data': []},
    'bastien': {'fname': "bastien-solution.txt", 'data': []},
    'static': {'fname': "static-solution.txt", 'data': []},

    #'procassini': {'fname': "proca-solution.txt", 'data': []},
    #'freq100': {'fname': "freq-100-solution.txt", 'data': []},
    #'freq50': {'fname': "freq-50-solution.txt", 'data': []},
    #'freq25': {'fname': "freq-25-solution.txt", 'data': []},
}

for k, cfg in configs.items():
    configs[k]['data'] = get_params(os.path.join(dir, cfg['fname']))

fig, ax = build_figures(sys.argv[1], babfiles, configs)

plt.tight_layout()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2], dpi=300)
else:
    plt.show()

