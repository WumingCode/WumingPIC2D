# Wuming PIC2D
Two-dimentional, special relativistic, electromagnetic particle-in-cell simulation code for general puposes in space and astrophysical plasmas.

## Features
* Buneman-Boris particle pusher
* Implicit FDTD scheme for the Maxwell equation (Hoshino, ApJ, 2013)
* Esirkepov's charge conservation scheme for the current deposit with the 2nd-order shape function (Esirkepov, CPC, 2001)
* Hybrid parallelization with MPI and OpenMP
* SIMD optimization
* JSON-based metadata with MPI-IO raw data output
* Python scripts for HDF5 format convertor and quicklook

## Requirements
* Fortran compiler
  - We prepare setting files for Makefile for different compilers, including Intel Fortran, GCC-Fortran, Fujitsu compiler
  - Because of the requirement of JSON-Fortran library used in this code (https://github.com/jacobwilliams/json-fortran), GCC-Fortran's version should be greater than 4.9.

* MPI library
  - We have tested with MPICH and Open MPI
  - The code works with the vender's MPI library on Fujitsu's supercomputer systems
   
* Python
  - The code generates raw binary data and JSON based metadata files. A Python script in each physics directory converts them to HDF5 files
  - A sample python script is prepared for quick look of the results

## Installation
```bash
$ git clone git@github.com:WumingCode/WumingPIC2D.git
```

## Code structure
``` 
WumingPIC2D
├── Makefile
├── README.md
├── common
│   └── common files of PIC algorithms
├── common.mk
├── compiler-fujitsu.mk
├── compiler-gcc.mk
├── compiler-intel.mk
├── compiler.mk
├── include
│   └── directory for module files
├── lib
│   └── directory for common library
├── proj
│   ├── shock
│   │   └── collsion-less shock simulation setup files
│   └── weibel
│       └── Weibel instability simulation setup files
└── utils
    └── utility files for MPI-IO and JSON output
```

## Preparation
1. Choose either of `comiler-*.mk` files depending on your compiler environment and copy it to `compiler.mk`

   ``` bash
   $ cp compiler-gcc.mk compiler.mk
   ```

2. Compile common library

   ``` bash
   $ make
   ```

3. Compile in a physics problem directory

   ``` bash
   $ cd proj/weibel
   $ make
   ```

## Physics calculation
Please refer to `README.md` in each physics problem in `proj/*`

