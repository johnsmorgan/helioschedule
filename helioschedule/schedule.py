"""
Optimized scheduling of observations
"""

import argparse
import csv
import json
import numpy as np
from yaml import safe_load
from helioschedule.beam_interp import Beams

round_to_8 = lambda x: int(x - x % 8)


def get_min_max(solar_noon_gps, offset):
    mintime = round_to_8(solar_noon_gps - offset)
    maxtime = mintime + round_to_8(offset * 2)
    return mintime, maxtime


def get_all_steps(solar_noon_gps, mintime, maxtime, solar_offset):
    return np.linspace(mintime, maxtime, solar_offset * 2 // 8, endpoint=False)


def get_flags(flags, mintime, maxtime):
    """
    Select relevant flags
    flags - start/stop in gps seconds
    solar noon (gps sections)
    offset (seconds)
    """
    return [(o[0], o[1]) for o in flags if o[1] > mintime if o[0] < maxtime]


def flags_to_mask(flags, n_steps, mintime):
    """
    convert flags (in gpstime) to boolean mask
    """
    flag_mask = np.zeros(n_steps, dtype=bool)
    for f in flags:
        flag_mask[(f[0] - mintime) // 8 : (f[1] - mintime) // 8] = True
    return flag_mask


def num_unflagged(array, idx, forward=True):
    """
    find the length of run of False after/before idx
    inclusive of idx
    """
    if idx < 0:
        return 0
    if idx >= len(array):
        return 0
    if array[idx] is True:
        return 0
    if forward is True:
        arr = array[idx:]
    else:
        arr = array[:idx][::-1]
    if np.all(arr != True):
        return len(arr)
    return np.where(arr == True)[0][0]


def add_out_dict_fields(field, out_dict):
    out_dict["ha_idx_%s" % field] = np.nan
    out_dict["beam_%s" % field] = -1
    out_dict["sun_attenuation_%s" % field] = np.nan
    out_dict["starttime_%s" % field] = np.nan
    out_dict["target_sensitivity_%s" % field] = np.nan
    out_dict["unflagged_before_%s" % field] = np.nan
    out_dict["unflagged_after_%s" % field] = np.nan
    return out_dict



class Scheduler:
    def setup_day(self, ref_time_gps):
        """
        set up all members required for each day's observation:
        self.all_obstimes
        self.flag_mask
        """
        # Calculate and stop time based on solar noon
        mintime, maxtime = get_min_max(ref_time_gps, self.conf["solarOffset"])
        self.all_obstimes = get_all_steps(ref_time_gps, mintime, maxtime, self.conf["solarOffset"])

        obs_list = get_flags(self.flags, mintime, maxtime)
        self.flag_mask = flags_to_mask(obs_list, len(self.all_obstimes), mintime)

    def setup_obs(self, solar_noon_gps, ha, dec, dec_sun, beam_chan, regen_sun=True):
        """
        set up beams for a new pointing
        self.sun_beam
        self.target_beam

        FIXME: might be able to eliminate `c`.
        Also need to think about whether beam_chan should be an attribute
        """
        # generate sun beam, only if we are on a new day or a new observing frequency (chan)
        if regen_sun is True:
            self.sun_beam = self.beams.interpolate_beam_2d(
                self.beams.beam_str_to_idx(beam_chan),
                solar_noon_gps,
                self.all_obstimes,
                None,
                dec_sun,
            )
            self.sun_filter = self.sun_beam < 10 ** self.conf["solarAttenuationCutoff"]
        # Next, the target field grid
        self.target_beam = self.beams.interpolate_beam_2d(
            self.beams.beam_str_to_idx(beam_chan),
            solar_noon_gps,
            self.all_obstimes,
            ha,
            dec,
        )

    def schedule_day(
        self, solar_noon_gps, local_noon_str, has, decs, dec_sun, ref_time_gps=None, out_dict=None
    ):
        if out_dict is None:
            out_dict = {}
        if ref_time_gps is None:
            ref_time_gps = solar_noon_gps
        self.setup_day(solar_noon_gps)
        flag_mask_convolved = np.convolve(self.flag_mask, np.ones(75, dtype=bool), "same")
        # FIXME: hard-coded indices
        self.flag_mask[:37] = True
        self.flag_mask[-38:]= True
        target_mask = np.zeros_like(self.flag_mask)
        target_filter = ~target_mask.reshape(1, -1)
        flag_filter = ~flag_mask_convolved.reshape(1, -1)
        beam_chan = None
        for c in self.conf["priority"]:
            out_dict = add_out_dict_fields(c, out_dict)
            if has[c] is None:
                continue
            if (beam_chan is None) or (self.conf["fields"][c]["beam_chan"] != beam_chan):
                beam_chan = self.conf["fields"][c]["beam_chan"]
                regen_sun = True
            else:
                regen_sun = False
            self.setup_obs(solar_noon_gps, ha=has[c], dec=decs[c], dec_sun=dec_sun, beam_chan=beam_chan, regen_sun=regen_sun)

            # applying sun_filter to target_beam will return a ravelled array.
            # we need to be able to identify the original location of our peak sensitivity within target_beam
            ha_grid, beam_grid = np.meshgrid(
                np.arange(self.target_beam.shape[1]), np.arange(self.target_beam.shape[0])
            )

            try:
                flat_idx = np.nanargmax(self.target_beam[self.sun_filter & flag_filter & target_filter])
            except ValueError:
                print("Warning, no observation meets criteria for %s day %s" % (local_noon_str, c))
                continue
            ha_idx = ha_grid[self.sun_filter & flag_filter & target_filter][flat_idx]
            beam_idx = beam_grid[self.sun_filter & flag_filter & target_filter][flat_idx]
            obstime = self.all_obstimes[ha_idx]
            # FIXME: hard-coded indices
            target_filter[0, ha_idx - 75 : ha_idx + 75] = False
            self.flag_mask[ha_idx - 37 : ha_idx + 38] = 2

            out_dict["beam_%s" % c] = self.beams.df["beams"].dims[1][0][beam_idx]
            out_dict["ha_idx_%s" % c] = ha_idx
            # FIXME: hard-coded index
            out_dict["starttime_%s" % c] = obstime - 296
            out_dict["sun_attenuation_%s" % c] = np.log10(self.sun_beam[beam_idx, ha_idx])
            out_dict["target_sensitivity_%s" % c] = self.target_beam[beam_idx, ha_idx]

        for c in self.conf["priority"]:
            ha_idx = out_dict["ha_idx_%s" % c]
            if not np.isnan(ha_idx):
                # FIXME: hard-coded indices
                out_dict["unflagged_before_%s" % c] = num_unflagged(self.flag_mask, ha_idx - 38, False)
                out_dict["unflagged_after_%s" % c] = num_unflagged(
                    self.flag_mask.astype(bool), ha_idx + 38, True
                )
        return out_dict

class DayScheduler(Scheduler):
    def __init__(self, conf, flags=[]):
        self.conf = conf
        self.flags = flags
        self.beams = Beams(self.conf["files"]["beams"])


class SemesterScheduler(Scheduler):
    def __init__(self, conf_filename):
        self.conf_filename = conf_filename
        self.conf = safe_load(open(conf_filename))
        self.days = list(csv.DictReader(line for line in open(self.conf["files"]["targets"])))
        self.flags = json.load(open(self.conf["files"]["flags"]))
        self.flag_img = np.zeros(
            (len(self.days), self.conf["solarOffset"] * 2 // 8), dtype=np.uint8
        )
        self.beams = Beams(self.conf["files"]["beams"])
        self.observations = []

    def schedule(self):
        for d, day in enumerate(self.days):
            solar_noon_gps = int(round(float(day["local_noon_gps"])))
            out_dict = {}
            for key in ("local_noon_str", "local_noon_lst"):
                out_dict[key] = day[key]
            print("scheduling around local noon", day["local_noon_str"])
            out_dict = self.schedule_day(
                solar_noon_gps,
                day["local_noon_str"],
                has={
                    c: float(day["ha_%s" % c]) if day["ha_%s" % c] != "" else None
                    for c in self.conf["priority"]
                },
                decs={
                    c: float(day["dec_%s" % c]) if day["ha_%s" % c] != "" else None
                    for c in self.conf["priority"]
                },
                dec_sun=float(day["dec_sun"]),
                out_dict=out_dict,
            )
            self.observations.append(out_dict)

    def write_observations(self):
        with open(self.conf["files"]["observations"], "w") as csvfile:
            fieldnames = sorted(self.observations[0].keys(), key=lambda k: k[::-1])
            fieldnames.insert(0, fieldnames.pop(-1))
            fieldnames.insert(0, fieldnames.pop(-1))
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")

            writer.writeheader()
            for out_dict in self.observations:
                writer.writerow(out_dict)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    scheduler = SemesterScheduler(args.infile)
    scheduler.schedule()
    scheduler.write_observations()


if __name__ == "__main__":
    main()
