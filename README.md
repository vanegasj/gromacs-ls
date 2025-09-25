# **GROMACS-LS: Local Stress and Elasticity Calculations from MD Simulations**


GROMACS-LS is a heavily modified version of GROMACS (<https://www.gromacs.org/>) designed to compute the local stress and elasticity tensors in 1-, 2-, or 3-dimensions from equilibrium molecular dynamics (MD) simulations. Most of the local stress/elasticity calculations in GROMACS-LS are performed by the MDStress library (<https://github.com/vanegasj/mdstress-library>). 

## **1. New development version**

---

We have new versions of GROMACS-LS and MDStress library that are ready for testing. Some of these new features include:

- Compute the local elasticity tensor in addition to the local stress tensor using the stress-stress fluctuation formula (see the papers at the bottom of the page). This feature is limited at the moment to force-fields with up to three-body potentials such as coarse-grained models. Four-body dihedral potentials will be implemented in the future.
- CPU Parallelization using threads (does not work with MPI, so can only be used in a single node/machine). The code is not optimized for GPUs or SIMD.
- Ability to compute the stress and elasticity while running a simulation on the fly, so there's no need to save the trajectory first and then post-process the trajectory. GROMACS-LS now also creates checkpoint files to restart local stress and elasticity calculations.
- Inclusion of the header-only eigen library (https://gitlab.com/libeigen/eigen) to speed-up calculations and remove the dependence on the LAPACK library.

### 1.1 Installation
External requirements:

- FFTW3 (double precission, libfftw3.so)
- MDStress library (libmdstress.so)
- CMake
- Python 3 with Numpy and Scipy to use the `tensortools` analysis tool

GROMACS-LS needs to be compiled in double precission and therefore the FFTW3 library also needs to be compiled and installed in double precission and as a shared library. For FFTW3, it is easiest to let GROMACS-LS automatically download, compile, and install the library for you. The GROMACS-LS package available from the Downloads menu is bundled with the MDStress library, so it it not necessary to install the library separately. 

If you install FFTW3 with a linux distribution such as Ubuntu or Fedora, you will also need to install the development (header) packages. If you have installed the FFTW3 library in a non-standard location (i.e. other than /usr/lib or /usr/local/lib), then before running cmake you should export the variables FFTW3_ROOT_DIR and CMAKE_PREFIX_PATH, e.g.

```
export FFTW3_ROOT_DIR=/path/to/fftw3
```

Now the installation of the custom GROMACS-LS package. Download the latest development package from [Releases](https://github.com/vanegasj/gromacs-ls/releases)

Configure, make and install:
```
tar -xvf gromacs-ls*.tar.gz
cd gromacs-ls*
mkdir build
cmake ../ -DCMAKE_INSTALL_PREFIX=/path/to/gromacs-ls
```
Add `-DGMX_BUILD_OWN_FFTW=ON` to the cmake line above if you want automatic build and configuration of the FFTW3 library.
```
make
make install
```
By default, the GROMACS-LS binary is called `gmx_LS` to distinguish it from other GROMACS installations present in your system.

### 1.2 Basic usage for computing the stress and elasticity on the fly
For local elasticity calculations, we recommend that one first simulate the system under constant temperature and pressure (*NPT*) for >100 ns and then run a constant volume and temperature (*NVT*) simulation based on a frame that most closely resembles the average system dimensions from the *NPT* run. You may get an incorrect result if you run the local elasticity analysis under *NPT* conditions! Please also read our papers at the end of the page for other important options in the gromacs mdp file.
To compute the bulk or local stress for a simulation, generate a tpr file using the familiar step of
```
gmx_LS grompp -c file.pdb -f grompp.mdp -o nvt.tpr
```
followed by an mdrun
```
gmx_LS mdrun -deffnm nvt -ols nvt
```
When the run is completed, `gmx_LS mdrun` will produce six output files in addition to the typical trajectory output:

1. `nvt_stress.dat0` - stress tensor.
2. `nvt_eltot.dat0` - total 'intrinsic' elasticity tensor, $c_{ijkl}$.
3. `nvt_eltothooke.dat0` - total 'extrinsic' elasticity tensor, $\tilde{c}_{ijkl}$, what one would measure within the context of Hooke's law. The two tensors are equal under conditions of zero stress, but may be significantly different for systems under high pressure or non-uniform systems such as lipid membranes where there is a large internal local pre-stress even if the ambient stress/pressure is small.
4. `nvt_elfluct.dat0` - fluctuation contribution to the elasticity tensor
5. `nvt_elkin.dat0` - kinetic contribution to the elasticity tensor
6. `nvt_elborn.dat0` - Born contribution to the elasticity tensor

See our paper at the end of the page for more details on these quantities.

Relevant options/flags that can be passed to `gmx_LS mdrun`:
```
 -ols     [<.dat>]           (localstress.dat)
 -localsgrid <real>         (0.1)
           Spacing for local stress grid (default = 0.1 nm)
 -localsiw <real>           (0.0001)
           Apply an impulsive correction for plain cutoff potentials (VDW or
           elec) using a delta function with a finite width (default = 0.0001
           nm)
 -lsgridx <int>             (0)
           Set the local stress grid size in the x direction (default use
           box[XX][XX]/localsgrid)
 -lsgridy <int>             (0)
           Set the local stress grid size in the y direction (default use
           box[YY][YY]/localsgrid)
 -lsgridz <int>             (0)
           Set the local stress grid size in the z direction (default use
           box[ZZ][ZZ]/localsgrid)
 -lscont <string>           (all)
           Select which contribution to write to output (default = all): all,
           vdw, coul, angles, bonds, dihp, dihi, dihrb, lincs, settle, shake,
           cmap, vel, none
 -lsfd   <string>           (ccfd)
           Select the type of force decomposition to be used: ccfd (covariant
           central force decomposition, default), ncfd (non-covariant central
           force decomposition), or gld (Goetz-Lipowsky decomposition)
 -[no]lsdispcor             (yes)
           Include contribution from dispersion correction.
 -[no]lspbc                 (no)
           Correct periodic boundary conditions in mdstress library. Typically
           not needed as gromacs passes the PBC-corrected positions and distances.
 -lsmindihang <real>        (0)
           Don't include dihedral local stress contributions if the sin(|phi|)
           is less than this factor. Use this flag if there is a dihedral
           potential (e.g. CHARMM36 lipid FF) that has been parametrized with
           a min/max that is not 0/Pi and the stress profiles show large noise
           that does not converge with additional frames. A -lsmindihang value
           of 0.0005 is typically sufficient to fix this problem.
 -lsskip <int>              (1)
           Only compute the local stress every nth frame
```

- Compute the total stress and elastic coefficients for a bulk system on the fly:
```
gmx_LS mdrun -deffnm nvt -ntmpi 64 -dd 4 4 4 -lsgridx 1 -lsgridy 1 -lsgridz 1 -lsskip 50 -ols nvt -localsiw 0.00001
```

- Compute the local stress and elasticity profiles on the fly along the *z* dimension with a grid spacing of 0.1 nm:
```
gmx_LS mdrun -deffnm nvt -ntmpi 64 -dd 4 4 4 -lsgridx 1 -lsgridy 1 -localsgrid 0.1 -lsskip 50 -ols nvt -localsiw 0.00001
```

- To restart a crashed calculation, relaunch with the same command and add the `-cpi` flag as one would do with conventional GROMACS.

The local stress and elasticity `.dat0` files produced by GROMACS-LS are stored in a simple binary format in double precision:
```
datatype     = sizeof(int)*1
box          = sizeof(double)*9
gridx        = sizeof(int)*1
gridy        = sizeof(int)*1
gridz        = sizeof(int)*1
tensor       = sizeof(double)*dsize*gridx*gridy*gridz
```
where `dsize` has values of 1 for density fields, 3 for vector fields, 9 for stress fields, and 36 for elasticity fields. These binary `.dat0` files can be easily converted and manipulated using the included `tensortools` python utility, e.g.:
```
tensortools -f nvt_stress.dat0 -o nvt_stress.txt
```

### 1.3 Known issues/limitations
- The thread-MPI CPU-only versions of GROMACS 2016 have a bug where after some time the `gmx mdrun` executable goes into a zombie state and stops producing any output even though the program appears to be running. This affects our GROMACS-LS code under some hardware combinations and we have been unable to pinpoint the source of the bug as the program doesn't actually crash or produce any errors or debugging output (If you know how to fix this bug let us know!!). If you notice that the program has stopped producing any new checkpoint files after 7-8 hrs and other trajectory files remain unchanged, then you are affected by this bug. Luckily, there's a simple work-around, which is to limit the time for the calculation and use the checkpointing feature to restart the calculation. Here's an example of bash code that executes `gmx_LS` in 6 hr intervals until it reaches the end of the simulation or runs 50 attempts:
```
runs=0
while [ ! -f final.pdb ] && [ -f nvt.tpr ] && [ $runs -le 50 ]
do
    gmx_LS mdrun -deffnm nvt -c final.pdb -cpi -maxh 6
    ((runs++))
done
```
- Elasticity calculations are limited at the moment to force-fields with only three-body potentials such as coarse-grained models (e.g., MARTINI). Four-body dihedral potentials will be implemented in the future.
- Calculation of long-range electrostatic interactions with PME are not included in the local stress/elasticity calculations.

## **2. Previous stable version based on GROMACS 2016.3**

---

The stable version of GROMACS-LS is used as a post-processing tool to analyze an existing trajectory and output the time-averaged local stress. Trajectories need to include both positions and velocities at every frame to be analyzed.

### 2.1 Installation
External requirements:

- FFTW3 (double precission, libfftw3.so)
- LAPACK (liblapack.so)
- MDStress library (libmdstress.so)
- CMake
- Python 3 with Numpy and Scipy to use the `tensortools` analysis tool


GROMACS-LS needs to be compiled in double precission and therefore the FFTW3 and LAPACK libraries also need to be compiled and installed in double precission and as shared libraries. For the FFTW3, it is easiest to let GROMACS-LS automatically download, compile, and install the library for you. The GROMACS-LS package available from the Downloads menu is bundled with the MDStress library, so it it not necessary to install the library separately. 

If you don't have LAPACK installed, download a recent verson from here: http://www.netlib.org/lapack (v 3.5.0 works well). Compiling LAPACK with the autoconf tools can be hassle, but recent versions can be easily compiled and installed with CMake:

```
tar -zxvf lapack-3.5.0.tgz
cd lapack-3.5.0
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/path/to/lapack ../
make
make install
```

If you install FFTW3 and LAPACK with a linux distribution such as Ubuntu or Fedora, you will also need to install the development (header) packages. If you have installed the FFTW3 or LAPACK libraries in a non-standard location (i.e. other than /usr/lib or /usr/local/lib), then before running cmake you should export the variables FFTW3_ROOT_DIR and CMAKE_PREFIX_PATH, e.g.

```
export FFTW3_ROOT_DIR=/path/to/fftw3
export CMAKE_PREFIX_PATH=/path/to/lapack
```

Now the installation of the custom GROMACS-LS package. Download the latest stable package from [Releases](https://github.com/vanegasj/gromacs-ls/releases)

```
wget gromacs-ls-2016.3-Dec-28-2019.tar.gz

```

Configure, make and install:
```
tar -xvf gromacs-ls*.tar.gz
cd gromacs-ls*
mkdir build
cmake ../ -DCMAKE_INSTALL_PREFIX=/path/to/gromacs-ls
```
Add `-DGMX_BUILD_OWN_FFTW=ON` to the cmake line above if you want automatic build and configuration of the FFTW3 library.
```
make
make install
```
By default, the GROMACS-LS binary is called `gmx_LS` to distinguish it from other GROMACS installations present in your system.

### 2.2 Basic usage for post-processing of a trajectory

For calculating the stress tensor in an existing MD simulation, we need a trajectory file (`.trr`) that contains both positions and velocities at the same points in time. In cases where system components drift over time, you may want to center your molecule/group of interest at every frame ([see the FAQ](FAQ.md)). Note that because of finite size effects, re-centering may change the resulting stress profile, so only do this if you are certain that your system is drifting.

The stress tensor is obtained by "rerunning" the trajectory with the `gmx_LS mdrun -rerun`. Since the `-rerun` option of `gmx_LS mdrun` outputs new `.log`, `.edr`, and `.trr` files, we suggest you first create a new folder and
analyze the trajectory within it:
```
mkdir stress
cd stress
gmx_LS mdrun -s ../topol.tpr -rerun ../traj.trr -lsgridx 1 -lsgridy 1 -localsgrid 0.1 -ols localstress
```

This will produce a file called `localstress.dat0` which can be processed with the included `tensortools` python utility:
```
tensortools -f localstress.dat0 -o localstress.txt
```

To save space you can recreate the `.tpr` file by setting all the output
`nst*` variables to 0 in the `.mdp` file:

```
nstxout   = 0
nstvout   = 0
nstlog    = 0
nstenergy = 0
nstxtcout = 0
```

If your simulation was performed with PME or another long-range electrostatic method, you will need to create a new `.tpr` file with the `coulombtype` set to `Cut-off` or `Reaction-field` and a cutoff radius
of at least 2.0-2.2 nm (see the Ref. 2 for more details). Note that this does not change your original trajectory but the way forces are calculated when analyzing each frame. Do not change the VdW cut-off as this will give incorrect results! The `gmx_LS grompp` utility has been modified to allow the use of different cutoff radii for electrostatic and VdW interactions.

### 2.3 Performance and parallelization

This version of the code only runs serially, but each frame in a trajectory can be analyzed completely independent of every other frame. Therefore, one method to parallelize the analysis is to separate the trajectory into smaller chunks with equal number of frames (using `gmx_LS trjconv -split`), and analyze each chunk separately. You can then average all of the resulting stress files together using the `tensortools` utility. The advantage of this method is that it always scales linearly although it takes a bit of work to set up. CAUTION! When doing this make sure to use the same `.tpr` file and same grid spacing to analyze each chunk of the trajectory, otherwise the grid sizes may differ and the average will fail.

Also, the cost of the local stress calculation grows significantly as the grid dimensions grow. Therefore, do not run a full 3D calculation unless it is completely necessary. If you only want a 1D stress profile along the z-dimension, you can use the flags `-lsgridx 1 -lsgridy 1`.

## **3. References/citations**

---

If you use GROMACS-LS or MDStress library in your research, please read and cite the following:

### 3.1 Local stress calculations:
```
@article{vanegas_importance_2014,
	title = {Importance of force decomposition for local stress calculations in biomembrane molecular simulations},
	volume = {10},
	number = {2},
	url = {http://pubs.acs.org/doi/abs/10.1021/ct4008926},
	doi = {10.1021/ct4008926},
	journal = {J. Chem. Theory Comput.},
	author = {Vanegas, Juan M and Torres-Sánchez, Alejandro and Arroyo, Marino},
	month = feb,
	year = {2014},
	pages = {691--702}
}

@article{torres-sanchez_examining_2015,
	title = {Examining the mechanical equilibrium of microscopic stresses in molecular simulations},
	volume = {114},
	number = {25},
	journal = {Phys. Rev. Lett.},
	author = {Torres-Sánchez, Alejandro and Vanegas, Juan M and Arroyo, Marino},
	url = {http://link.aps.org/doi/10.1103/PhysRevLett.114.258102},
	doi = {10.1103/PhysRevLett.114.258102},
	month = jun,
	year = {2015},
	pages = {258102}
}

@article{torres-sanchez_geometric_2016,
	title = {Geometric derivation of the microscopic stress: {A} covariant central force decomposition},
	volume = {93},
	number = {C},
	url = {http://dx.doi.org/10.1016/j.jmps.2016.03.006},
	doi = {10.1016/j.jmps.2016.03.006},
	journal = {J. Mech. Phys. Solids},
	author = {Torres-Sánchez, Alejandro and Vanegas, Juan M and Arroyo, Marino},
	month = aug,
	year = {2016},
	pages = {224--239}
}
```

### 3.2 Local elasticity calculations:
```
@article{lewis_elast_2025a,
	title = {Microscopic elasticity from MD part I: Bulk solid and fluid systems},
	volume = {TBD},
	number = {TBD},
	journal = {TBD},
	author = {Lewis, Andrew L and Himberg, Benjamin and Torres-Sánchez, Alejandro and Vanegas, Juan M},
	month = jan,
	year = {2025},
	pages = {TBD}
}
@article{lewis_elast_2025b,
	title = {Microscopic elasticity from MD part II: Liquid interfaces and lipid membranes},
	volume = {TBD},
	number = {TBD},
	journal = {TBD},
	author = {Lewis, Andrew L and Himberg, Benjamin and Torres-Sánchez, Alejandro and Vanegas, Juan M},
	month = jan,
	year = {2025},
	pages = {TBD}
}
```
