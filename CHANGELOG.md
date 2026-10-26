# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the makefiles build the C++20 module `smspp.satellitesblock` when
  `SMSppMOD` is not empty (e.g., `make SMSppMOD=1`), with clang++ 16 or
  later [see `lib/makefile-mod` of SMS++]

- the C++20 module `smspp.satellitesblock` (`src/SatellitesBlock.cppm`),
  which exports every name that the headers declare at namespace scope and
  re-exports `smspp`, so that a source file can write `import
  smspp.satellitesblock;` in place of the `#include` of the headers; it is
  built with the CMake option `SMSpp_BUILD_MODULE`, as the target
  `SMS++::SatellitesBlock_module`

### Changed

- the pipeline builds the module also with its makefiles, after CMake, and
  runs on a change of the makefiles too [see `ci/makefile` of the umbrella]

- `SMSpp_BUILD_MODULE` no longer stops the configuration when the C++20
  module cannot be built, i.e., with a CMake older than 3.28, a generator
  other than Ninja or Visual Studio, a compiler other than Clang 16 (with
  `clang-scan-deps`) or MSVC 19.34, or a module it imports that is not
  there: the option is turned off with a warning that says why, and the
  build goes on with the headers [see `cmake/CheckModule.cmake` of SMS++]

### Fixed

- the makefiles of the tests and of the tools build also with g++: they
  take the option that stops at the first error and the one of the
  debugging information, which clang++ and g++ spell differently
  (`-ferror-limit=1` and `-glldb` against `-fmax-errors=1` and `-ggdb`),
  from `SMSppSWERR` and `SMSppSWDBG` of `lib/makefile-OS` of SMS++


## [0.1.0] - 2026-10-09

### Added

- ConstellationBlock, SatelliteBlock, SingleTargetBlock, MultiTargetBlock
  (and its MultiTargetBlockv2 variant) for the continuous Satellite
  Constellation Design Problem, with the SatelliteSolver that solves a
  SatelliteBlock by inspection

- DiscreteConstellationBlock and DiscreteSatelliteBlock for the version of
  the problem where the observability threshold of each satellite takes one
  of finitely many levels, with the DiscreteSatelliteSolver

- the two instances of the module in data/txt

- the unit test of the module on the core alone: the Block are in the
  factory, generating their abstract representation issues no Modification,
  and a second load() gives the same representation

- module structure generated from ModuleTemplate: version derived from the
  git tag, versioned shared library, relative install RPATH, the symbol that
  keeps the module linked, CI building the commit of its pipeline

### Fixed

- no Modification is issued while constructing the abstract representation
  (eNoMod/eNoBlck in the generate_* methods): this made the Blocks unusable
  as sub-Block of a Lagrangian decomposition

- a second load() rebuilds the Block from scratch: the sub-Block are
  deleted rather than leaked, and the flags saying which part of the
  abstract representation is there are reset, so that it is generated again

- load() reports a malformed instance with an exception (no altitude in
  [400, 1400] km, fewer than three where three are needed, latitude out of
  [-90, 90], fewer than two discretization points) instead of reading out
  of range or printing an error and going on, and the library prints
  nothing on the standard output

- the void observation constraints of DiscreteConstellationBlock name each
  Variable once, as a LinearFunction requires

- SatelliteSolver and DiscreteSatelliteSolver release their mutex on every
  return of compute(), read the coefficients of the Objective by Variable
  rather than by position, count its constant term, and have a solution
  only after a successful compute()

- is_feasible() of SatelliteBlock, SingleTargetBlock, MultiTargetBlockv2 and
  DiscreteSatelliteBlock checks the Variable and the constraints instead of
  always answering false, and the Solution of the Blocks read and write the
  Variable instead of being empty

- MultiTargetBlockv2 deleted twice the LinearFunction of two observation
  constraints that were never added to the model

[Unreleased]: https://gitlab.com/smspp/satellitesblock/-/compare/0.1.0...develop
[0.1.0]: https://gitlab.com/smspp/satellitesblock/-/tags/0.1.0
