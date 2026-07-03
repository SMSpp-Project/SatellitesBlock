# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- no Modification is issued while constructing the abstract representation
  (eNoMod/eNoBlck in the generate_* methods): this made the Blocks unusable
  as sub-Block of a Lagrangian decomposition

### Added

- Lagrangian decomposition test (LagrangianDualSolver + BundleSolver +
  *MILPSolver) on a small multi-target SCDP instance

## [0.1.0] - 2026-07-02

### Added

- ConstellationBlock, SatelliteBlock, SingleTargetBlock, MultiTargetBlock
  (and its MultiTargetBlockv2 variant) for the continuous Satellite
  Constellation Design Problem, with the SatelliteSolver heuristic

- module structure generated from ModuleTemplate

### Changed

- sources realigned to the SMS++ style
