# Changelog

## [v0.0.3] - 2026-06-05

### Added

* Adds a Dockerfile for the project with samtools pre-installed
* Includes (and thus supersedes) feat: add onyx project arg #3 and fix: importlib to use joinpath #4
* Adds a --database arg so you can re-use a taxaplease database when running in nextflow
* Will accept local BAM files as well as climb ids as input
  * How this works is by checking if each input is_file() - if it is we assume its a local BAM. This makes it easier to use in nextflow as we can stage in BAMs from S3 to processes but otherwise climb ids are still downloaded.


### Changed

* Changes CI so docker file is built off of main and publishes image to GHCR
* Now check for samtools in PATH before falling back to previous conda behaviour
* Removed some duplicate taxaPlease and onyx inits

### Fixed

* Fixes a dodgy lookup, but lookup.txt has non taxid entries in the taxid field #5 still remains unfixed broadly

## [v0.0.2] - 2026-05-08

### Added

None

### Changed

* Orthomyxoviridae segment numbers automatically mapped to segment names.
* Dropdown entries are sorted alphabetically
* Script now prints more status info while running

### Fixed

Tolerate empty figlists when grouping data.


## [v0.0.1] - 2026-05-08

### Added

Initial commit

### Changed

None

### Fixed

None
