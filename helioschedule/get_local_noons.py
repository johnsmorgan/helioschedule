"""
Compute local solar noon in LST and UTC and write to file.
"""

import argparse
from datetime import datetime, timedelta, timezone
from yaml import safe_load
from skyfield import almanac
from skyfield.api import load, N, E, wgs84
from astropy import units as u
from astropy.coordinates import Longitude


class SolarTransiter:
    def __init__(self, lat_deg_n, lon_deg_e, elevation_m):
        self.eph = load("de421.bsp")
        self.lat, self.lon, self.el = lat_deg_n, lon_deg_e, elevation_m

        self.sun = self.eph["Sun"]
        self.observer = self.eph["Earth"] + wgs84.latlon(self.lat * N, self.lon * E)
        self.ts = load.timescale()

    def get_transit(self, t_start, t_end):
        times = almanac.find_transits(self.observer, self.sun, t_start, t_end)
        return times


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    transiter = SolarTransiter(
        lat_deg_n=conf["lat"], lon_deg_e=conf["lon"], elevation_m=conf["alt"]
    )

    t_start = transiter.ts.from_datetime(
        datetime.combine(
            conf["startDate"],
            datetime.min.time(),
            tzinfo=timezone(timedelta(hours=conf["timezone"])),
        )
        + timedelta(hours=12)
    )
    t_end = transiter.ts.from_datetime(
        datetime.combine(
            conf["stopDate"],
            datetime.min.time(),
            tzinfo=timezone(timedelta(hours=conf["timezone"])),
        )
        + timedelta(hours=12)
    )

    times = transiter.get_transit(t_start, t_end)

    with open(conf["files"]["noons"], "w") as f:
        print("local_noon_str,local_noon_lst,local_noon_gps", file=f)
        for t, time in enumerate(times):
            local_noon_lst = Longitude(
                time.gmst * u.hourangle + conf["lon"] * u.deg, wrap_angle=180 * u.deg
            ).deg
            print(
                time.utc_iso()[:-1],
                "%.3f" % local_noon_lst,
                "%.1f" % time.to_astropy().gps,
                sep=",",
                file=f,
            )


if __name__ == "__main__":
    main()
