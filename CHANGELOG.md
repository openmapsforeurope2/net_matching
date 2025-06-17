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


## [1.0.0] - 2025-03-24
### Added
- Initial release of the project

### Changed
- NTR

### Fixed
- NTR