# Wuming PIC2D
Two-dimentional, special relativistic, electromagnetic particle-in-cell simulation code for general puposes in space and astrophysical plasmas.

## Features
* Solves the Vlasov-Maxwell equations by the particle-in-cell method
* Buneman-Boris method for the equation of motion of particles.
* Implicit FDTD scheme for the Maxwell equation (Hoshino, ApJ, 2013)
* Esirkepov's charge conservation scheme for the current deposit with the 2nd-order shape function (Esirkepov, CPC, 2001)
* Written in Fortran 90/95
* Hybrid parallelization by MPI and OpenMP
* SIMD optimization and efficient cache usage
* MPI-IO raw data output with JSON-based metadata
* Python scripts for HDF5 format convertor and quicklook

## Requirements
* Fortran compiler
  - We prepare setting files for Makefile for different compilers, including Intel Fortran, GCC-Fortran, Fujitsu compiler
  - Because of the requirement of JSON-Fortran library used in this code (https://github.com/jacobwilliams/json-fortran), **GCC-Fortran's version must be greater than 4.9**.

* MPI library
  - We have tested with MPICH and Open MPI
  - The code works with the vender's MPI library on Fujitsu's supercomputer systems
   
* Python [OPTIONAL]
  - The code generates raw binary data and JSON-based metadata files. A Python script in each physics directory converts them to HDF5 files as a post-process
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
├── include
│   └── directory for module files
├── lib
│   └── directory for common library
├── proj
│   ├── shock
│   │   └── collsion-less shock simulation setup files and scripts for post process
│   └── weibel
│       └── Weibel instability simulation setup files and scripts for post process
└── utils
    └── utility files for MPI-IO and JSON output
```

## Preparation
1. Move to the installed directory
   ``` bash
   $ cd ./WumingPIC2D
   ```
   
2. Copy one of `comiler-*.mk` files depending on your compiler environment to `compiler.mk`

   ``` bash
   $ cp compiler-gcc.mk compiler.mk
   ```

3. Make a common library

   ``` bash
   $ make
   ```

4. Make an excutable file in a physics problem directory

   ``` bash
   $ cd proj/weibel
   $ make
   ```
## Physics calculation
Please refer to `README.md` in each physics problem in `proj/*`

* [Collision-less shock](proj/shock/README.md)
* [Weibel instability](proj/weibel/README.md)

## Credits
WumingPIC2D code uses 
* [JSON-Fortran](https://github.com/jacobwilliams/json-fortran) API for reading/writing JSON files from Fortran.
* [Amano's MPI-IO, JSON, HDF5 utitlity files](https://github.com/amanotk)
* Anymore?

## LICENSE
TO BE ADDED (any ideas?)

