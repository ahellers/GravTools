# API Reference

This reference covers the public Python API of GravTools.

## Package structure

```
gravtools/
├── models/                       # Core domain model
│   ├── campaign.py               # Campaign — top-level container
│   ├── survey.py                 # Survey — one field session
│   ├── station.py                # Station — known-point database
│   ├── gravimeter.py             # Gravimeters — instrument metadata
│   ├── lsm_diff.py               # LSM: differential observations
│   ├── lsm_nondiff.py            # LSM: non-differential observations
│   ├── vg_lsm.py                 # LSM: vertical gravity gradient
│   ├── mlr_bev_legacy.py         # BEV legacy MLR processing
│   └── atmosphere_correction.py  # Atmospheric pressure corrections
├── tides/                        # Tidal corrections
│   ├── longman1959.py            # Longman (1959) tidal model
│   ├── correction_time_series.py # TSF-based external corrections
│   └── tide_data_tfs.py          # TSoft TSF file parser
├── CG5_utils/                    # Scintrex CG-5 file parser
├── CG6_utils/                    # Scintrex CG-6 file parsers (LynxLG v1/v2, CG6 Solo)
├── settings.py                   # Enumerations, lookup dicts, GUI constants
└── const.py                      # Physical constants
```

## Data model hierarchy

```
Campaign
├── Survey (one per field session / instrument day)
│   └── Observation / Setup  (stored in obs_df, a pandas DataFrame)
├── Station                  (stat_df — known-point coordinates, datum flags)
├── Gravimeters              (instrument metadata & scale factors)
├── LSM run results          (list of LSMDiff / LSMNonDiff / VGLSM results)
└── CorrectionTimeSeries     (external tidal / loading corrections)
```

## Conventions

| Convention | Detail |
|---|---|
| Gravity unit | µGal throughout; DataFrame columns use the `_mugal` suffix |
| Height unit | metres; columns use the `_m` suffix |
| Angles | degrees in DataFrames (`_deg`); radians internally where noted |
| Timestamps | timezone-naive UTC `datetime` / `numpy.datetime64` |
| Persistence | Campaigns are serialised as pickle files (`.pkl`, protocol 4) |
| Optional dep. | `geopandas` is guarded by `try/except`; check `_has_geopandas` flags |

## Where to start

- **Loading data** → [`Campaign`](campaign.md) and [`Survey`](survey.md)
- **Station management** → [`Station`](station.md)
- **Least-squares adjustment** → [Parameter Estimation](lsm_diff.md)
- **Tidal & atmospheric corrections** → [Corrections](tides.md)
- **Configuration** → [`Settings`](settings.md)

