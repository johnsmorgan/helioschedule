"""
Compute local solar noon in LST and UTC and write to file.
"""

import argparse
from datetime import datetime, timedelta, timezone
from yaml import safe_load
from skyfield import almanac
from skyfield.api import load, Topos
from astropy import units as u
from astropy.coordinates import Longitude

INTERVAL = 1  # days


def main():
    eph = load("de421.bsp")
    ts = load.timescale()

    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    sun = eph["sun"]

    location = Topos(
        latitude_degrees=conf["lat"], longitude_degrees=conf["lon"], elevation_m=conf["alt"]
    )

    t_start = ts.from_datetime(
        datetime.combine(
            conf["startDate"],
            datetime.min.time(),
            tzinfo=timezone(timedelta(hours=conf["timezone"])),
        )
        + timedelta(hours=12)
    )
    t_end = ts.from_datetime(
        datetime.combine(
            conf["stopDate"],
            datetime.min.time(),
            tzinfo=timezone(timedelta(hours=conf["timezone"])),
        )
        + timedelta(hours=12)
    )

    f = almanac.meridian_transits(eph, sun, location)
    times, y = almanac.find_discrete(t_start, t_end, f)

    with open(conf["files"]["noons"], "w") as f:
        print("local_noon_str,local_noon_lst", file=f)
        for t, time in enumerate(times):
            if y[t]:
                local_noon_lst = Longitude(
                    time.gmst * u.hourangle + conf["lon"] * u.deg, wrap_angle=180 * u.deg
                ).deg
                print(time.utc_iso()[:-1], "%.3f" % local_noon_lst, sep=",", file=f)


if __name__ == "__main__":
    main()
