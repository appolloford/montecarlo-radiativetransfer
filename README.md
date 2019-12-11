# montecarlo-radiativetransfer
Simple Monte Carlo Radiative Transfer code (RT course assignment)

### Compile
g++ --std=c++11 -o mcrt.exe mcrt.cpp

The structure changed after assignment 3.
For assignment 1 and 2, use : g++ --std=c++11 -o mcrt.exe test.cpp

### Usage / Benchmark Test

+ Assignment 1: ./mcrt.exe --singleruntest (--inhomogeneous)
+ Assignment 2: ./mcrt.exe (--albedo=\<albedo\>) (--reflect=\<surface albedo\>) (--for-back-ratio=\<ratio\>) (--inhomogeneous)
    + ./mcrt.exe --albedo=0.8 --reflect=0.3 --for-back-ratio=0.7
+ Assignment 3: ./mcrt.exe (--Tau_z==(\<\tau^*\>) (--albedo=\<albedo\>) (--reflect=\<surface albedo\>) (--nPhotons=\<number\>)
    + ./mcrt.exe --Tau_z=1.0 --albedo=0.5 --reflect=0.8 --nPhotons=10000
+ Assignment 4: ./mcrt.exe (--Tau_z==(\<\tau^*\>) (--albedo=\<albedo\>) (--reflect=\<surface albedo\>) (--nPhotons=\<number\>) --testHG (--HGAsymParam=\<param\>)
    + ./mcrt.exe --Tau_z=5.0 --albedo=0.8 --reflect=0.7 --nPhotons=10000 --testHG

### Plot
