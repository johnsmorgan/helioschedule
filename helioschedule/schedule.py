"""
Optimized scheduling of observations
"""

import argparse
import csv
import numpy as np
from h5py import File
from scipy.interpolate import interp1d
from yaml import safe_load

INTERVAL_DEGREES = 3.3424596978682920e-02  # interval in degrees for fine HA search
N = 75  # observation length in units of INTERVAL_DEGREES


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    targets = csv.DictReader(line for line in open(conf["files"]["targets"]))
    df = File(conf["files"]["beams"], "r")
    fine_scale = np.arange(
        np.min(df["beams"].dims[3][0]), np.max(df["beams"].dims[3][0]), INTERVAL_DEGREES
    )

    obs_ha = []

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

    for target in targets:
        beam_chan = None
        out_dict = {}
        day_has = []
        for key in ("local_noon_str", "local_noon_lst"):
            out_dict[key] = target[key]
        for c in conf["priority"]:
            if target["ha_%s" % c] == "":
                out_dict["ha_%s" % c] = np.nan
                out_dict["beam_%s" % c] = -1
                out_dict["sun_attenuation_%s" % c] = np.nan
                out_dict["target_sensitivity_%s" % c] = np.nan
                continue
            print(target["local_noon_str"][:10], c)
            if beam_chan is not None and conf["fields"][c]["beam_chan"] is beam_chan:
                # Only reconstruct sun_beam if we have switched frequency
                pass
            else:
                try:
                    beam_chan = conf["fields"][c]["beam_chan"]
                    f = np.argwhere(
                        df["coarse_chans"][...] == conf["fields"][c]["beam_chan"].encode("ascii")
                    )[0][0]
                except IndexError:
                    print("can't find %s in beam file" % conf["fields"][c]["beam_chan"])
                sun_dec = float(target["dec_sun"])
                sun_dec_idx = neighbours(df["beams"].dims[2][0][...], sun_dec)
                sun_beam = lin_interp(
                    df["beams"][f, :, sun_dec_idx[0], :],
                    df["beams"][f, :, sun_dec_idx[1], :],
                    sun_dec - df["beams"].dims[2][0][sun_dec_idx[0]],
                )
                ha_grid = None
                if "flags" in conf.keys():
                    flag_filter = np.zeros(df["beams"].dims[3][0].shape, dtype=bool)
                    for f, flag in enumerate(conf["flags"]):
                        start = float(target["start_flag_%d" % (f + 1)]) - (N * INTERVAL_DEGREES)
                        stop = float(target["stop_flag_%d" % (f + 1)]) + (N * INTERVAL_DEGREES)
                        flag_filter = flag_filter | (
                            (df["beams"].dims[3][0][...] > start)
                            & (df["beams"].dims[3][0][...] < stop)
                        )
                        print(flag_filter)
                        print(np.where(flag_filter, np.nan, 1.0)[None, :].shape)
                    sun_beam *= np.where(flag_filter, np.nan, 1.0)[None, :]
            target_dec = float(target["dec_%s" % (c)])
            target_dec_idx = neighbours(df["beams"].dims[2][0][...], target_dec)

            # only roll to integer degree
            target_ha_idx = int(round(float(target["ha_%s" % c])))
            # produce sun_beam and target beam
            # both are 2D arrays (pointing, ha)

            target_beam = lin_interp(
                df["beams"][f, :, target_dec_idx[0], :],
                df["beams"][f, :, target_dec_idx[1], :],
                target_dec - df["beams"].dims[2][0][target_dec_idx[0]],
            )

            target_beam = np.roll(target_beam, target_ha_idx, axis=1) / np.expand_dims(
                df["broadness"][f, ...], 1
            )
            sun_filter = sun_beam < 10 ** conf["solarAttenuationCutoff"]

            # applying sun_filter to target_beam will return a ravelled array.
            # we need to be able to identify the original location of our peak sensitivity within target_beam
            if ha_grid is None:
                ha_grid, beam_grid = np.meshgrid(
                    np.arange(target_beam.shape[1]), np.arange(target_beam.shape[0])
                )
                # also need a filter to stop subsequent observations from clashing
                target_filter = np.ones(target_beam.shape, bool)

            try:
                flat_idx = np.nanargmax(target_beam[sun_filter & target_filter])
            except ValueError:
                print(
                    "Warning, no observation meets criteria for %s target %s"
                    % (target["local_noon_str"], c)
                )
                out_dict["ha_%s" % c] = np.nan
                out_dict["beam_%s" % c] = -1
                out_dict["sun_attenuation_%s" % c] = np.nan
                out_dict["target_sensitivity_%s" % c] = np.nan
                # continue
                break

            ha_idx = ha_grid[sun_filter & target_filter][flat_idx]
            beam_idx = beam_grid[sun_filter & target_filter][flat_idx]
            ha = df["beams"].dims[3][0][ha_idx]
            beam = df["beams"].dims[1][0][beam_idx]
            # refine ha
            #
            tt_rounded = int(np.ceil(conf["timeTweakDegrees"]))
            fine_slice = slice(ha_idx - tt_rounded, ha_idx + tt_rounded + 1)

            fine_ha_interp = interp1d(
                df["beams"].dims[3][0][fine_slice],
                sun_beam[beam_idx, fine_slice],
                kind="quadratic",
            )
            fine_has = fine_scale[np.abs(fine_scale - ha) < conf["timeTweakDegrees"]]
            fine_sun_beam = fine_ha_interp(fine_has)

            # blank out existing with np.inf
            for ha_ in day_has:
                fine_sun_beam = np.where(
                    np.abs(fine_has - ha_) < N * INTERVAL_DEGREES, np.inf, fine_sun_beam
                )
            # if 'flags' in conf.keys():
            # for f, flag in enumerate(conf['flags']):
            # start = float(target['start_flag_%d' % (f+1)])
            # stop = float(target['stop_flag_%d' % (f+1)])
            # fine_sun_beam = np.where((fine_has<stop)&(fine_has>start), np.nan, fine_sun_beam)
            if np.all(fine_sun_beam == np.inf):
                print("all inf!")
            ha = fine_has[np.argmin(fine_sun_beam)]
            flag_range = neighbours(df["beams"].dims[3][0], ha)
            target_filter[:, flag_range[0] - 2 : flag_range[1] + 2] = False
            day_has.append(ha)

            out_dict["beam_%s" % c] = beam
            out_dict["ha_%s" % c] = ha
            out_dict["sun_attenuation_%s" % c] = np.log10(np.min(fine_sun_beam))
            out_dict["target_sensitivity_%s" % c] = target_beam[beam_idx, ha_idx]
        obs_ha.append(out_dict)
        print()

    with open(conf["files"]["observations"], "w") as csvfile:
        fieldnames = sorted(obs_ha[0].keys(), key=lambda k: k[::-1])
        fieldnames.insert(0, fieldnames.pop(-1))
        fieldnames.insert(0, fieldnames.pop(-1))
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for out_dict in obs_ha:
            writer.writerow(out_dict)


if __name__ == "__main__":
    main()
