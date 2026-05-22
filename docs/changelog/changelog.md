# Release notes and Changelog

More details are provided along with version tags on [https://github.com/ahellers/GravTools/tags](https://github.com/ahellers/GravTools/tags).

## 0.3.3 (2025-01-08)
- New dialog for drift determination based on the reduced observations as shown in the observations plot. One drift 
polynomial is fitted per station. This feature is useful for calculating the drift correction coefficient for Scintrex 
CG5 and CG6 meters based on stationary observations. 

## 0.3.4 (2025-02-07)
- Fixed bug in the data export dialog.

## 0.3.5 (2025-09-24)
- Fixed bug: The survey name was recognized as instrument S/N when loading CG6 solo observation files.

## 0.3.6 (2026-01-02)
- Loading a new json gravimeter file now triggers the calculation of updated gravity reductions.
- Future pandas warning concerning timezone aware datetime64 types fixed.
- Bugfix: Problem when loading stations with integer-type statio names from CSV files fixed (resulted in duplicated entries in the station list).
- Handling of station filters in the GUI (stations tab) improved. 

## 0.3.7 (2026-03-27)
- Fixed syntaxWarning that was caused by regex expressions in the CG5 observation file loading method.
- The setup IDs are now created based on the observation reference time and the survey name (combined hash). This allows multiple observations at the same time in different surveys within one campaign. Now it is possible to analyze campaigns with two ore more gravimeters being used in parallel. 
- Fixed a bug that caused an error when reading negative dhf values from the instrument height column in LynxLG CG6 observation files.

## 0.3.8 (2026-05-20)
- Fixed compatibility with pandas >= 3.0: replaced in-place `.values[:] = 0` writes with direct Series assignment to comply with Copy-on-Write semantics (survey.py); fixed `ValueError: assignment destination is read-only` in Longman (1959) tidal correction caused by `.values` returning a read-only array from a pandas Index (longman1959.py).
- Fixed compatibility with NumPy >= 2.0: replaced deprecated `np.NAN` with `np.nan` (mlr_bev_legacy.py) and replaced broken `.values.astype(np.int64) / 10**9` datetime conversion with a resolution-independent `to_unix_seconds()` helper (survey.py, correction_time_series.py).
- Packaging migrated from setup.cfg to pyproject.toml (PEP 621); explicit minimum version requirements added for all dependencies.
- Code quality refactoring: file renames, unused imports/variables removed, bare `except` clauses narrowed, optional geopandas import centralised, f-string modernisation, context managers for file I/O.
- Added an API reference, changelog and a landing page for documentation build with mkdocs (html). The documentation is accessible via the about menu in the GUI.
