#!/usr/bin/zsh

# 52 * 10^7 / 1Gf
P=10649600
T=500

W=`python -c "print((52*10**-2 * $P))"`
C=`python -c "print($W / $P*10**2)"`

make -j4

# static workload
./bin/lb -W $W -p $P -i $T -C $C -w 0 -u 0 && python script/lbplot.py "$\iota(t) = 0.1$" results/constant.pdf > results/constant-perfs.txt

./bin/lb -W $W -p $P -i $T -C $C -w 0 -u 1 && python script/lbplot.py "$\iota(t) = \frac{1}{0.4*t+1}$" results/sublinear.pdf > results/sublinear-perfs.txt

./bin/lb -W $W -p $P -i $T -C $C -w 0 -u 2 && python script/lbplot.py "$\iota(t) = 0.02*t$" results/linear.pdf > results/linear-perfs.txt

./bin/lb -W $W -p $P -i $T -C $C -w 0 -u 3 && python script/lbplot.py "$\iota(t) = -(0.1*(t \% 17)) + 8$" results/autocorrect.pdf > results/autocorrect-perfs.txt

# dynamic workload
./bin/lb -W $W -p $P -i $T -C $C -w 1 -u 1 && python script/lbplot.py "$\iota(t) = \frac{1}{0.4*t+1}$" results/sublineardynamic.pdf > results/sublinear-perfs-dynamic.txt
