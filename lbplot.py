import matplotlib.pyplot as plt
import sys
with open('best-cpu-time-%d.txt'%int(sys.argv[1]), 'r') as f:
    lines = [float(l) for l in f.readlines()]
plt.plot(lines)
plt.show()
