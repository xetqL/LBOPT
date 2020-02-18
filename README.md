# lb_exhaustive_search
Lightning fast code for computing load balancing scenario from application parameters

# install
``git clone https://github.com/xetqL/lb_exhaustive_search.git --recursive``

``cmake -DCMAKE_BUILD_TYPE=Release .``

``make``

# run
```console
./lb_exhaustive_search -h

Usage: lb_exhaustive_search  [-hVv] -W=<DOUBLE> [-d <INT>] -p=<INT> -i=<INT> -C=<DOUBLE> [-s <INT>]

  Required options:
  -W, --W0=DOUBLE           Initial workload(default: 0)
  -p, --processor=INT       Number of processors(default: 0)
  -i, --iteration=INT       Number of iteration to simulate(default: 2776629248)
  -C, --lbcost=DOUBLE       Load balancing cost(default: 0)

  Optional options:
  -h, --help                print this help and exit
  -V, --version             print version and exit
  -d, --deltaW=INT          Select the workload increase rate function: (0) dw=1.0, (1) dw=5.0/(1.0+0.2*i, (2) dw=1.0 + 0.2*i, (3) dw=sin(0.2*i)(default: 0)
  -s, --solution=INT        Number of output solution from branch and bound(default: 1)
  -v, --verbose             (1) output iteration by iteration, (2) produce output 
```  
# example
```console  
./lb_exhaustive_search -W 10 -p 10 -C 20 -d 1 -s 3 -i 10 -v                                                                                    255 ↵
=================================================================================================================================
Results from search using branch and bound: 
            Solution 0:              iteration  10 -> { 0 0 0 0 1 0 0 0 0 0 } = 90.11871
            Solution 1:              iteration  10 -> { 0 0 0 0 0 1 0 0 0 0 } = 90.46755
            Solution 2:              iteration  10 -> { 0 0 0 1 0 0 1 0 0 0 } = 92.65171

Results using mathematical formulation: 
sqrt[2*C/dW*(1-1.0/P)]: Tau=3.81941; iteration  10 -> { 0 0 0 1 0 0 0 1 0 0 } = 93.51931
sqrt[2*C/dW*(1+1.0/P)]: Tau=3.45479; iteration  10 -> { 0 0 1 0 0 0 1 0 0 1 } = 109.02764
sqrt(2*C/dW):           Tau=3.62341; iteration  10 -> { 0 0 0 1 0 0 1 0 0 0 } = 92.65171
=================================================================================================================================
Solution (0) -------------------------------------------------
iteration   1 -> { 0 } = 1.00000
iteration   2 -> { 0 0 } = 6.16667
iteration   3 -> { 0 0 0 } = 14.30952
iteration   4 -> { 0 0 0 0 } = 24.68452
iteration   5 -> { 0 0 0 0 1 } = 47.27083
iteration   6 -> { 0 0 0 0 1 0 } = 52.35714
iteration   7 -> { 0 0 0 0 1 0 0 } = 59.48891
iteration   8 -> { 0 0 0 0 1 0 0 0 } = 68.32522
iteration   9 -> { 0 0 0 0 1 0 0 0 0 } = 78.60383
iteration  10 -> { 0 0 0 0 1 0 0 0 0 0 } = 90.11871
Solution (1) -------------------------------------------------
iteration   1 -> { 0 } = 1.00000
iteration   2 -> { 0 0 } = 6.16667
iteration   3 -> { 0 0 0 } = 14.30952
iteration   4 -> { 0 0 0 0 } = 24.68452
iteration   5 -> { 0 0 0 0 0 } = 36.79563
iteration   6 -> { 0 0 0 0 0 1 } = 59.65972
iteration   7 -> { 0 0 0 0 0 1 0 } = 64.79654
iteration   8 -> { 0 0 0 0 0 1 0 0 } = 71.82729
iteration   9 -> { 0 0 0 0 0 1 0 0 0 } = 80.46061
iteration  10 -> { 0 0 0 0 0 1 0 0 0 0 } = 90.46755
Solution (2) -------------------------------------------------
iteration   1 -> { 0 } = 1.00000
iteration   2 -> { 0 0 } = 6.16667
iteration   3 -> { 0 0 0 } = 14.30952
iteration   4 -> { 0 0 0 1 } = 36.58333
iteration   5 -> { 0 0 0 1 0 } = 41.63492
iteration   6 -> { 0 0 0 1 0 0 } = 48.90873
iteration   7 -> { 0 0 0 1 0 0 1 } = 72.02282
iteration   8 -> { 0 0 0 1 0 0 1 0 } = 77.22024
iteration   9 -> { 0 0 0 1 0 0 1 0 0 } = 84.18048
iteration  10 -> { 0 0 0 1 0 0 1 0 0 0 } = 92.65171
```

  
