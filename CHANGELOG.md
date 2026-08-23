# Changelog

## 0.2.1

First release since 0.1.6. Version 0.2.0 was tagged but never published.

### Added

- Guided enumeration: steer the search with a scoring function through a genetic
  algorithm or Bayesian optimizer. Install with `pip install 'mol-healer[opt]'`.
  See [Guided Enumeration](docs/guided-enumeration.md).
- `BBID` result columns carrying the catalog identifier of each building block,
  read from the `id`, `catalog_id`, or `molport id` SDF property.
- `Molport_Fast_Delivery` named building block source.
- Configurable logging: `configure_logging()` and the `HEALER_LOG_LEVEL`
  environment variable. Progress bars are now controlled separately through
  `show_progress`.
- `max_retro_nodes`, a node budget for the retrosynthesis tree.
- Web interface: building block libraries show their size, jobs report the stage
  they are in, and finished runs report molecule count and elapsed time.
- Documentation under [docs/](docs/README.md) and a demo notebook.

### Changed

- `get_results()` keeps one row per distinct route, so the same product can
  appear more than once with different building blocks or reactions. Scripts
  that assumed unique products will see more rows.
- `verbose` is deprecated in favour of `show_progress`.

### Fixed

- Building block preprocessing recorded every reaction for every block.
  **Libraries processed with an earlier version must be re-run through
  `preprocess-bb`** — until then, reaction matching is wrong.
- Importing `healer` now applies the RDKit patch, so the documented quick start
  works without an extra import.
- `enumerate()` can be called repeatedly on `SiteHEALER` and `FragmentHEALER`.
- `SiteHEALER` no longer mutates the query molecule passed to it, and no longer
  shares one rules dictionary across instances.
- An unparsable reaction template is reported and skipped instead of stopping
  the whole library from loading.
