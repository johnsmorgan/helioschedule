# helioschedule

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

Schedule MWA observations close to the Sun

## Instructions
All of the metadata required to schedule the entire observing run (as well as the locations of generic data files) are all recorded in a `yaml`-format configuration file. See the examples.

A typical scheduling run consists of 4 discrete commands. Each command reads from a single file and produces a single output file. These files are all specified in the configuration file. There is a single python script in the `helioschedule` directory corresponding to each of these commands, which can be run as follows.

Run `get_local_noons [yaml_file]` - this determines the UTC time (and Local Sidereal Time) of local noon for each day. All scheduling is done relative to this LST.

Run `convert_coordinates [yaml_file]` to get position in ecliptic and equatorial coordinates for each pointing for each day

Run `schedule [yaml_file]` This is where the actual optimised scheduling happens.

Run `generate_observations [yaml_file]` to generate the commands that the MWA operations team can use to schedule the observations.

It is suggested to use a simple Makefile to orchestrate this process. See the examples.

## Project Organization

```
├── LICENSE            <- Open-source license if one is chosen
├── Makefile           <- Makefile with convenience commands like `make data` or `make train`
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── README.md      <- README
│   ├── azel.json      <- Mapping from delays to nominal azimuth and elevation.
│   ├── generate_beams <- Generate beams in HA/Decl. from `mwa_lookup_beam` lookup table
│
├── docs               <- A default mkdocs project; see mkdocs.org for details
│
├── examples           <- Set of example scheduling scenarios
│
├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
│
├── pyproject.toml     <- Project configuration file with package metadata for helioschedule
│                         and configuration for tools like black
│
├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
│                         generated with `pip freeze > requirements.txt`
│
├── setup.cfg          <- Configuration file for flake8
│
└── helioschedule                <- Source code for use in this project.
    │
    ├── __init__.py    <- Makes helioschedule a Python module
    │
    ├── get_local_noons.py             <- Pre-calculate local noons over the whole observing timerange
    ├── convert_coordinates.py         <- Convert target coordinates to RA and Decl.
    ├── schedule.py                    <- Carry out scheduling
    ├── generate_observation_script.py <- Convert to a script containing commands to schedule on MWA
```

--------

