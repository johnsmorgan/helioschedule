"""
Calculate the beams, averaged over the bandwidth, utilising
pre-calculated beams from 
https://github.com/johnsmorgan/mwa_pb_lookup


Hour Angle (fastest moving axis)
Declination
Both are in degrees with 1 degree resolution
"""

import json
import argparse
import numpy as np
from scipy.interpolate import RectBivariateSpline

from h5py import File


# NB the following two are copied from lookup_beam.py
def trap(n):
    """
    return normalised trapezium
    """
    x = np.ones((n,), dtype=np.float32)
    x[0] = x[-1] = 0.5
    x /= sum(x)
    return x


def coarse_range(chans, coarse_str):
    int_chans = [int(c) for c in coarse_str.split("-")]
    edge_freq_hz = (int_chans[0] * 1280000 - 640000, int_chans[1] * 1280000 + 640000)

    lower = np.argwhere(chans == edge_freq_hz[0]).flatten()
    upper = np.argwhere(chans == edge_freq_hz[1]).flatten()
    if len(lower) == 0:
        raise IndexError("No match for lower coarse chan %d" % int_chans[0])
    if len(upper) == 0:
        raise IndexError("No match for upper coarse chan %d" % int_chans[1])
    return lower[0], upper[0] - lower[0] + 1


LAT = -26.703319  # degrees


azel = json.load(open("azel.json"))

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Input hdf5 file")
parser.add_argument(
    "coarse_chans", help="comma-separated list of coarse chan ranges: e.g. 121-132,157-180"
)
parser.add_argument("outfile", help="Output hdf5 file")
args = parser.parse_args()

coarse_chans = args.coarse_chans.split(",")
n_freq = len(coarse_chans)
print(coarse_chans)

# Generate grid of all HA and Declination coordinates
ha_scale = np.linspace(
    -179.5, 179.5, 360
)  # this has the advantage of avoiding the singularities at 0 and 180 degrees
dec_scale = np.linspace(-90, 90, 181)
has, decs = np.meshgrid(ha_scale, dec_scale)

# figure out corresponding azimuth and elevation
# for simplicity, use basic trig rather than trying to construct a time when RA0h is precisely on the meridian.
alt = np.degrees(
    np.arcsin(
        np.sin(np.radians(decs)) * np.sin(np.radians(LAT))
        + np.cos(np.radians(decs)) * np.cos(np.radians(LAT)) * np.cos(np.radians(has))
    )
)
# This doesn't quite work for the nadir (HA = 180deg), but who cares?
az = np.degrees(
    np.arccos(
        (np.sin(np.radians(decs)) - np.sin(np.radians(alt)) * np.sin(np.radians(LAT)))
        / (np.cos(np.radians(alt)) * np.cos(np.radians(LAT)))
    )
)
az = np.where(np.sin(np.radians(has)) < 0, az, 360 - az)

# next, read in each of the beams
# and interrogate for values at those coordinates
df_in = File(args.infile, "r")
sweetspots = df_in["sweetspot_number"][...]
n_sweetspots = len(sweetspots)

beam_data = np.zeros([n_freq, n_sweetspots] + list(az.shape))
beam_azs = np.zeros((n_freq, n_sweetspots))
beam_els = np.zeros((n_freq, n_sweetspots))
broadnesss = np.zeros((n_freq, n_sweetspots))

alt_in = df_in["alt_scale"][...]
az_in = df_in["az_scale"][...]

for c, chan_str in enumerate(coarse_chans):
    low_index, n_chan = coarse_range(df_in["chans"][...], chan_str)
    weights = trap(n_chan)
    for i in sweetspots:
        print("sweetspot %d/%d chan %d/%d" % (i + 1, n_sweetspots, c + 1, n_freq))
        beam_az = azel[str(i)][0]
        beam_el = azel[str(i)][1]
        print(beam_az, end=" ")
        print(beam_el)
        print(low_index, n_chan)
        print()

        # Use Sault weighting to combine polarisations
        d_chan = np.average(
            df_in["beams"][i, low_index : low_index + n_chan, :, :, :],
            axis=1,
            weights=df_in["beams"][i, low_index : low_index + n_chan, :, :, :] ** 2,
        )
        # Next combine frequencies
        d = np.average(d_chan, axis=0, weights=weights)

        d /= np.nanmax(d)

        # The above normalisation does not take into account that some beams are now wider than others;
        # The wider ones will let in more diffuse radiation, raising the tsys. We store this
        # parameter as the broadness; normalised to the zenith beam.
        broadnesss[c][i] = np.nansum(
            d * np.cos(np.radians(df_in["alt_scale"][...])).reshape((-1, 1))
        )

        d_interpolate = RectBivariateSpline(x=alt_in, y=az_in, z=np.nan_to_num(d))
        beam_data[c][i] = d_interpolate(alt, az, grid=False).reshape(az.shape)
        beam_data[c][i] = np.where(alt < 0, np.nan, beam_data[c][i])
        # print dmax_to_world(d)

# Normalise broadness by zenith pointing
broadnesss /= broadnesss[:, 0].reshape(-1, 1)

# We keep broadness as a separate parameter as we want to apply it when calculating the
# relative sensitivity to our target, but we will not apply it when
# calculating the attenuation at the location of the Sun
#
# A fuller treatment would take account of how the sensitivity depends on where the telescope
# is pointing on the sky. However this would complicated matters *considerably* compared
# to the approach we're taking.

azel_tuples = [(int(a[0]), int(a[1])) for a in np.swapaxes(np.vstack((beam_azs, beam_els)), 1, 0)]

# finally, store all of the relevant data in an hdf5 file
with File(args.outfile, "w") as df:
    df.create_dataset("beams", data=beam_data, compression="lzf", shuffle=True)

    df["beams"].dims[0].label = "freq"
    df.create_dataset("coarse_chans", data=[s.encode("ascii") for s in coarse_chans])
    df["beams"].dims.create_scale(df["coarse_chans"])
    df["beams"].dims[0].attach_scale(df["coarse_chans"])

    df["beams"].dims[1].label = "beam"
    df.create_dataset("sweetspot_number", data=sweetspots)
    df["beams"].dims.create_scale(df["sweetspot_number"])
    df["beams"].dims[1].attach_scale(df["sweetspot_number"])

    df["beams"].dims[2].label = "dec"
    df.create_dataset("dec_scale", data=dec_scale)
    df["beams"].dims.create_scale(df["dec_scale"])
    df["beams"].dims[2].attach_scale(df["dec_scale"])

    df["beams"].dims[3].label = "ha"
    df.create_dataset("ha_scale", data=ha_scale)
    df["beams"].dims.create_scale(df["ha_scale"])
    df["beams"].dims[3].attach_scale(df["ha_scale"])

    df.create_dataset("az", data=beam_azs)
    df.create_dataset("el", data=beam_els)
    df.create_dataset("broadness", data=broadnesss)
    df.create_dataset("delays", data=df_in["delays"][...])
