# Collisionless Shock

## How to run
In short, compile and edit a configuration file and run it!

Compiling by

```console
$ make
```

will create an executable `main.out`. This will read parameters from a configuration file in JSON format. You may copy a sample configuration file `config_sample.json` to `config.json`:

```console
$ cp config_sample.json config.json
```

and edit it as you like. By default, running the code via

```console
$ mpiexec -np 4 ./main.out
```

will try to read `config.json`.  
If you want, you may explicitly specify the filename with a command line argument:

```console
$ mpiexec -np 4 ./main.out somethingelse.json
```

in which case `sometingelse.json` will be read.

## Configuration file
See the detailed descriptions provided below.

- `config`
  - `datadir`  
    Data directory to which all the output will be saved.
  - `max_elapsed`  
    Maximum elapsed time. The code will stop automatically when the elapsed  
    time goes beyond this limit. A snapshot will be saved for restart.
  - `max_it`  
    Maximum number of iteration step.
  - `intvl_ptcl`  
    Interval of time step for entire particle data output.
  - `intvl_mom`  
    Interval of time step for for moment data output.
  - `intvl_orb`  
    Interval of time step for tracer particle data output.
  - `intvl_expand`  
    Interval of time step for expanding the box size in x direction.
  - `restart`  
    Snapshot to be read. If this is specified, the code will start from  
    this state. Otherwise, it will start from the initial condition.  
    Note that this will be overwritten by the code when it finishes.
- `parameter` : Parameters used for setting initial conditions.
  - `num_process`  
     Number of MPI processes. This must be the same as the number specified  
     in execusion with `mpiexec` or `mpirun`.
  - `n_ppc`  
     Number of particle per cell in the upstream region.
  - `n_x`  
     Number of grid in x direction.
  - `n_y`  
     Number of grid in y direction.
  - `n_x_ini`  
     Initial number of grid in x direction. The simulation box will expand  
     every `intvl_expand` step.
  - `u_injection`  
     Injection four-velocity from the boundary.
  - `mass_ratio`  
     Ion to electron mass ratio: m_i/m_e.
  - `sigma_e`  
     Ratio of electron cyclotron-to-plasma frequency squared.
  - `omega_pe`  
     Proper plasma frequency in the upstream.
  - `v_thi`  
     Thermal velocity of ions in the upsteram.
  - `v_the`  
     Thermal velocity of electrons in the upsteram.
  - `theta_bn`  
     Magnetic obliquity in the upstream or the shock angle in degrees.
  - `phi_bn`  
     Magnetic field orientation in the upstream with respect to the  
     simulation plane in degrees. The in-plane and out-of-plane magnetic  
     field configurations correspond to 0 and 90, respectively.
  - `l_damp_ini`  
     Length scale for the initial velocity profile in unit of grid size.  
     This is to smooth out initial disturbances and nothing physical.


## How it works

The code will produce a lot of files in a specified directory. The output should
appear in pairs of `*.json` and `*.raw` files. A `.json` file in JSON format
describes meta data and how the actual simulation data are stored in a `.raw`
file, which contains raw data in binary format.

The JSON files can be processed to generate HDF5 format files for data analysis
via a script `json2hdf5.py`. For instance in the working directory,

```console
$ python json2hdf5.py ./*.json
```

will process all JSON files in the current directory and generate HDF5 format
files for each JSON format file. You can use the generated HDF5 files for data
analysis.

When the code running time goes beyond a given elapsed time limit, it will save
a snapshot data and stop running. By specifying a snapshot in the configuration
file, you may restart the run.

By default, the configuration file will be overwritten by the code so that you
can restart the run next time with the exactly same command.

For instance, if you run the code via

```console
$ mpiexec -np 4 ./main.out
```

and the elapsed time limit is reached, it will overwrite the configuration file
`config.json` to properly set `restart` option.
In the next time, you may run again via

```console
$ mpiexec -np 4 ./main.out
```

then the previous snapshot data will be read automatically.

