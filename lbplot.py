import matplotlib.pyplot as plt
import numpy as np
import sys
def get_params(fname):
    with open(fname) as f:
        I, W, t, cum, dec, tim,C = f.readlines()
        W =   [float(x) for x in W.strip().split(' ')]
        cum2= [float(x) for x in cum.strip().split(' ')]
        dec2= [int(x)   for x in dec.strip().split(' ')]
        time2=[float(x) for x in tim.strip().split(' ')]
    return int(I), W, float(t), cum2, dec2, time2, float(C)


fmenon = "menon-solution.txt"
fbab   = "optimal-solution-0.txt"

I, W, tmen, cummen, decmen, timemen, C = get_params(fmenon)
I, W, tbab, cumbab, decbab, timebab, C = get_params(fbab)

fig, ax = plt.subplots(3,1, figsize=(20, 20))
ax[0].set_title(sys.argv[1])

ax[0].plot(W, label='Workload')
ax[0].set_ylabel('Application Workload')
ax[0].legend()

ax[1].plot(timebab, label='Branch and Bound')
ax[1].plot(timemen, label='U > C')
ax[1].set_ylabel('CPU Time')


ax[1].legend()

ax[2].plot(cumbab, ls='--', label='Branch and Bound')
ax[2].plot(cummen, ls='dotted', label='U>C')
ax[2].plot([C]*I, label='C')
ax[2].set_xlabel('iteration')
ax[2].set_ylabel('Cumulative imbalance time')
ax[2].legend()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()
