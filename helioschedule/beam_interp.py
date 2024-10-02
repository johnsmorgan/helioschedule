import numpy as np
from h5py import File
from scipy.interpolate import interp1d

SIDEREAL_FACTOR = 1.0027379093604878  # 86400/((1.0*units.sday).to(units.s).value)

arg_closest = lambda x, y: np.argmin(np.abs((x - y)))


def neighbours(arr, val):
    """
    return two closest values in arr to val
    assumes arr is sorted (lowest value first)
    """
    closest = arg_closest(arr, val)
    if arr[closest] > val:
        closest1 = closest - 1
        closest2 = closest
    else:
        closest1 = closest
        closest2 = closest + 1
    return closest1, closest2


def lin_interp(y1, y2, dx):
    """
    return linear interpolation of y1 and y2

    y1 and y2 are y(x1) and y(x2)

    dx controls the relative weighting of y1 and y2 in the interpolation
    """
    return y1 * (1 - dx) + y2 * dx


class Beams:
    def __init__(self, filename):
        self.df = File(filename, "r")
        # FIXME add some checking here that we are working
        # with the correct kind of file.

    def beam_str_to_idx(self, beam_str):
        """
        Convert a beam string (e.g. '121-132') to an index
        """
        beam_idx = np.argwhere(self.df["coarse_chans"][...] == beam_str.encode("ascii"))
        if beam_idx.shape[0] == 0:
            raise RuntimeError("course chan not in index!")
        return beam_idx[0][0]

    def interpolate_beam_2d(self, freq_idx, ha0_obstime, obstimes, ha_offset_deg, dec_deg):
        """
        Construct a 2D surface of beams consisting of all pointings for a series of times.
        - freq_idx - Frequency index (first axis of beam array)
        - obstime - times
        - ha0_obstime time (in units of obstimes)
        - ha_offset_deg - Hour angle offset in degrees from
        - Declination in degrees
        """

        dec_idx = neighbours(self.df["beams"].dims[2][0][...], dec_deg)
        beam_1deg = lin_interp(
            self.df["beams"][freq_idx, :, dec_idx[0], :],
            self.df["beams"][freq_idx, :, dec_idx[1], :],
            dec_deg - self.df["beams"].dims[2][0][dec_idx[0]],
        )
        beam_1deg = np.where(np.isnan(beam_1deg), 0, beam_1deg)
        # calculate gpstime corresponding to each HA for target
        if ha_offset_deg is None:
            # we are dealing with the Sun
            beam_x = ha0_obstime + 86400 * self.df["beams"].dims[3][0][...] / 360
        else:
            beam_1deg /= np.expand_dims(self.df["broadness"][freq_idx, ...], 1)
            ha_offset_s = 86400.0 * ha_offset_deg / 360.0 / SIDEREAL_FACTOR
            beam_x = ha0_obstime + SIDEREAL_FACTOR * (
                86400 * self.df["beams"].dims[3][0][...] / 360 + ha_offset_s
            )
        beam_interpolator = interp1d(x=beam_x, y=beam_1deg, kind="cubic", axis=1, fill_value='extrapolate')
        if obstimes.max() > beam_x.max():
            print(f"Warning! extrapolating beam by {obstimes.max()-beam_x.max()}s")
        return beam_interpolator(obstimes)
