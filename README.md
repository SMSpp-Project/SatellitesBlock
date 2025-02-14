# SatellitesBlock

Implementation of different `Blocks` for the Satellite Constellation 
Design Problem (SCDP) for the minimization of the number of satellites
in the constellation or of the sum of the maximum revisit times per target.

- `ConstellationBlock` is a collection of `SatelliteBlocks`, defining the
single satellites of the constellation, linked by the observation constraints,
imposing that each target should be observed by at least one satellite within
the revisit time. These `Blocks` are necessary to minimize their number of
satellites active in the constellation and able to observe the targets within
the revisit time (see `test_Constellation`).

- `MultiTargetBlock` is a collection of `SingleTargetBlocks`, defining the
single target to be observed, linked by the configuration constraints,
imposing that exactly one orbital configuration is selected for each satellite
in the constellation. These `Blocks` are necessary to minimize the maximum 
revisit time per target (see `test_Targets`)

### Bibliography

- L. Mencarelli. An MILP approach to minimize the maximum revisit time in 
the satellite constellation design problem. Working paper, 2025.
- L. Mencarelli. On the Lagrangian relaxation for the 
satellite constellation design problem. Working paper, 2025.

## Getting started

These instructions will let you build `SatellitesBlock` on your system.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `SatellitesBlock` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the
[tester to minimize the number of satellites and the maximum revisit time ]
(https://gitlab.com/smspp/tests/-/tree/develop/SatellitesBlock).
The makefiles in turn recursively include all the required other makefiles,
hence one should only need to edit the "main makefile" for compilation type
(C++ compiler and its options) and it all should be good to go. In case some
of the external libraries are not at their default location, it should only be
necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/satellitesblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

- **Luca Mencarelli**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
