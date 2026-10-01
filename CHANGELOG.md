# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
