# README #

This folder contains Stochastic 1D-Burgers, QG and SWE model coupled with DO code. It solves the QG problem with
same parameters as in __qg__ repository( see __JBIMAU/DO-methods/QG/matlab__ folder). QG code is
not parallelized. Burgers equation has been verified to give correct same
results in serial as well as in parallel. SWE model also run in parallel but
results are not verified.

## Build

You must have already built and installed the following libraries somewhere

* Trilinos (OpenMPI built)
* Hymls (OpenMPI)
* BOOST or TRNG library (for random number generation)

Set `TRILINOS_HOME`, `HYMLS_HOME` environment variables. If BOOST library is used then we are already expecting it to be installed system-wide and included in default system's paths. If TRNG library is used then set  `TRNG_HOME` environment variable.

Then  create a directory called __build__ in the root folder  then type the following in terminal:

````
cd build
CXX=mpicxx cmake ..
````

 To build stochastic BURGER, type the follwing

```
make BRGRDO
```

similarly, for stochastic QG and SWE type the following respectively

```
make QGDO
make SWEDO
```

to build all three together just type 

```
make
```

 ## Running the code

folder __paramsXML__ contains three directories containing the essential xml files for the respective problems. All three directories contain __params.xml__ (states the values related to mean solver) and __StochasticParams.xml__ ( DO solver parameter).  

## To Do
* Change `Topography` parameter in __params.xml__ for QG to `Bottom Topograpgy coefficient` or somethig similar. 
