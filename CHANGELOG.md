## [1.2.0] - 2026-09-18
### Added
- Added support for collapsing connecting lines (CL), including iterative and small CL collapsing.
- Added support for generating connecting points (CP) from undershoot detection.
- Added support for generating CP from CL using an improved geometry computation method.
- Added support for managing multiple connecting features between more than two countries.
- Added support for detecting and managing connections to a third country.
- Added support for managing multiple fictitious values in connecting line processing.
- Added support for managing existing connections close to international borders.
- Added a command-line parameter to specify the target database name.
- Added the `CL_MIN_LENGTH` parameter for watercourse links.

### Changed
- Improved the handling of connecting lines inside areas, including fictitious edges and incident edges.
- Improved the merging of connecting points, including CPs originating from the same CL and different merging distances.
- Changed the CP and vertex border snapping distances.
- Changed the output table suffix so that it no longer includes country codes and is based only on the user-defined suffix.
- Improved the processing of clean faces and antenna features.
- Consolidated connecting lines inside areas.
- Added in-memory storage of selected tables to reduce the number of SQL queries and improve processing performance.
- Updated the application for the new SOCLE version.
- Updated the Docker build and runtime environment.
- Adapted the application to the IGN-MUT deployment environment.

### Fixed
- Fixed cases where the projection of a connecting line onto an edge results in a point.
- Fixed issues related to SRID handling for CP and CL tables.
- Fixed issues in CP generation and merging.
- Fixed an infinite loop in connecting line processing inside areas.
- Fixed issues occurring when processing already matched data.


## [1.1.0] - 2025-06-17
### Added
- [configuration] configuration extended to new countries
- [all] implementation of NotDestroyedTools facilities not to iterate on destructed objects
- [step] added step CorrectCountryConnectivity at the very begining of the process to remedy source data connectivity issues

### Changed
- [parameters] refactoring/renaming
- [documentation] various corrections
- [CFeatConnectionOp] added method to treat the networks of the two countries being matched at the same time. Previous edges displacement were processed country by country at step ConnectionConnectingLines generating disconnection between networks and CL (connecting lines). Indeed, CL generation that occure in previous steps results in an intrication of the two countries networks. As a consequence the whole data must be treated in the step ConnectionConnectingLines as a single network.


### Fixed
- [CLInAreaGenerationOp] correction not to merged edges from same country
- [EdgeCleaningOp] added deletion instruction that was missing in method _cleanTinyEdges
- [FillFictitiousFieldOp] management of not boolean values for the field 'fictitious'


## [1.0.0] - 2025-03-24
### Added
- Initial release of the project

### Changed
- NTR

### Fixed
- NTR