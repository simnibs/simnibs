# SimNIBS

THIS REPOSITORY CONTAINS SIMNIBS 3.0 ONLY.

SIMNIBS 3.0 IS STILL IN ALPHA VERSION.

The main goal of SimNIBS is to calculate electric fields caused by Transcranial Electrical Stimulation (TES) and Transcranial Magnetic Stimulation (TMS).
The pipeline is divided in three parts:
1. Automatic segmentation of MRI images and meshing to create individualized head models
2. Calculation of electric fields through the Finite Element Method (FEM)
3. Post-processing of results for further analysis.

## Build Status
| Linux   | Windows    | OSX |
|---------|------------|-----|
| [![Build Status](https://dev.azure.com/simnibs/simnibs/_apis/build/status/Linux?branchName=master)](https://dev.azure.com/simnibs/simnibs/_build/latest?definitionId=4&branchName=master) | [![Build Status](https://dev.azure.com/simnibs/simnibs/_apis/build/status/Windows?branchName=master)](https://dev.azure.com/simnibs/simnibs/_build/latest?definitionId=5&branchName=master) |     |

## Getting Started
 
SimNIBS runs on 64bit Windows, Linux and OSX machines.
Please visit [the SimNIBS website](http://www.simnibs.org) for instructions into how to download and install SimNIBS.


## Authors
Please see [the SimNIBS website](http://www.simnibs.org) for a complete list of contributors.

## 3rd Party Files
We have included code or binaries from the following project to this repository:
* [Gmsh](www.gmsh.info)
* [meshfix](https://github.com/MarcoAttene/MeshFix-V2.1)
* [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
* [CAT12](http://www.neuro.uni-jena.de/cat/)
* [PETSc](https://www.mcs.anl.gov/petsc/)
* [HYPRE](https://www.mcs.anl.gov/petsc://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods)
* [MPICH](https://www.mpich.org/)
* [MSMPI](https://github.com/Microsoft/Microsoft-MPI)
* [CYGWIN](https://www.cygwin.com/)
* [pygpc](https://github.com/konstantinweise/pygpc)

For a full list of files and licenses, please see the [3RD-PARTY.md](3RD-PARTY.md) file
