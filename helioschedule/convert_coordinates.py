"""
Calculate coordinates of each observation
"""

import argparse
from datetime import datetime
from numpy import radians, sin, cos, arcsin, arctan2, ones
from yaml import safe_load

from astropy import io
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, Longitude, get_sun, GeocentricTrueEcliptic, GCRS
from astropy.table import Table, Column
from astropy.time import Time


def destination(theta, d, phi1=0.0, lambda1=0.0):
    """
    from http://www.movable-type.co.uk/scripts/latlong.html

    Given a starting lat and long (phi1, lambda1), bearing theta and angular distance d

    Returns new lat and long
    """
    phi2 = arcsin(sin(phi1) * cos(d) + cos(phi1) * sin(d) * cos(theta))
    lambda2 = lambda1 + arctan2(sin(theta) * sin(d) * cos(phi1), cos(d) - sin(phi1) * sin(phi2))
    return phi2, lambda2


def parse_time(t_string):
    """
    Parse time without date (i.e. HH:MM:SS) and return as number of seconds
    since midnight, returning as degrees
    """
    t = Angle(t_string, unit=u.hour)
    return t.deg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    noons = io.ascii.read(conf["files"]["noons"])

    times = Time(noons["local_noon_str"])
    start_time = Time(datetime.combine(conf["startDate"], datetime.min.time()))
    stop_time = Time(datetime.combine(conf["stopDate"], datetime.min.time()))
    good_times = (times > start_time) & (times < stop_time)
    times = times[good_times]

    out_table = noons[good_times]
    out_table = Table(out_table, masked=True, copy=False)  # convert to masked table

    sun_3d = get_sun(times)
    sun_eq = SkyCoord(sun_3d.ra, sun_3d.dec)

    sun_ra = Column(data=sun_eq.ra, name="ra_sun")
    sun_dec = Column(data=sun_eq.dec, name="dec_sun")
    out_table.add_columns([sun_ra, sun_dec])

    sun_ecliptic = sun_eq.geocentrictrueecliptic
    for target in conf["priority"]:
        if conf["fields"][target]["system"] == "Heliocentric":
            l1, phi = conf["fields"][target]["coordinates"]
            target_offset = destination(radians(phi), radians(l1))
            target_ecliptic = GeocentricTrueEcliptic(
                sun_ecliptic.lon + u.rad * target_offset[1],
                sun_ecliptic.lat + u.rad * target_offset[0],
            )
            target_eq = target_ecliptic.transform_to(GCRS())
            lon = Column(data=target_ecliptic.lon.deg, name="lon_%s" % (target))
            lat = Column(data=target_ecliptic.lat, name="lat_%s" % (target))
            ra = Column(data=target_eq.ra.deg, name="ra_%s" % (target))
            ha = Column(
                data=Longitude(target_eq.ra - sun_eq.ra, wrap_angle=180 * u.deg).deg,
                name="ha_%s" % (target),
            )
            dec = Column(data=target_eq.dec, name="dec_%s" % (target))
            out_table.add_columns([lon, lat, ra, dec, ha])
        elif conf["fields"][target]["system"] == "Ecliptic":
            ra_coord, dec_coord = conf["fields"][target]["coordinates"]
            ra = Column(data=ra_coord * ones(times.shape) * u.deg, name="ra_%s" % (target))
            ha = Column(
                data=Longitude(
                    ra_coord * ones(times.shape) * u.deg - sun_eq.ra, wrap_angle=180 * u.deg
                ).deg,
                name="ha_%s" % (target),
            )
            dec = Column(data=dec_coord * ones(times.shape) * u.deg, name="dec_%s" % (target))
            out_table.add_columns([ra, dec, ha])
        if "skip" in conf["fields"][target]:
            assert "offset" in conf["fields"][target], "target %s has 'skip' but no 'offset'"
            out_table["ha_%s" % target].mask = True
            out_table["ha_%s" % target].mask[
                conf["fields"][target]["offset"] :: conf["fields"][target]["skip"]
            ] = False

    if "flags" in conf.keys():
        noon_deg = parse_time([t.isot[11:] for t in times])
        sidereal = float(1.0 * u.day / u.sday)
        for f, flag in enumerate(conf["flags"]):
            assert flag["type"] == "daily", "only daily type flags supported"
            print(f"start:{flag['start']} stop:{flag['stop']}")
            start = parse_time(flag["start"])
            stop = parse_time(flag["stop"])
            print(f"start:{start} deg stop:{stop} deg")
            out_table["start_flag_%d" % (f + 1)] = sidereal * (start - noon_deg)
            out_table["stop_flag_%d" % (f + 1)] = sidereal * (stop - noon_deg)

    out_table.write(conf["files"]["targets"], format="csv", overwrite=True)


if __name__ == "__main__":
    main()
