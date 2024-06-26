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
    def schedule_day(self, solar_noon_gps, day, out_dict=None):
        beam_chan = None
        if out_dict is None:
            out_dict = {}
        # Calculate and stop time based on solar noon
        mintime, maxtime = get_min_max(solar_noon_gps, self.conf["solarOffset"])
        all_obstimes = get_all_steps(solar_noon_gps, mintime, maxtime, self.conf["solarOffset"])

        obs_list = get_flags(self.flags, mintime, maxtime)
        flag_mask = flags_to_mask(obs_list, len(all_obstimes), mintime)
        flag_mask_convolved = np.convolve(flag_mask, np.ones(75, dtype=bool), "same")
        target_mask = np.zeros_like(flag_mask)
        target_filter = ~target_mask.reshape(1, -1)
        flag_filter = ~flag_mask_convolved.reshape(1, -1)
        for c in self.conf["priority"]:
            out_dict = add_out_dict_fields(c, out_dict)
            if day["ha_%s" % c] == "":
                continue
            if beam_chan is not None and self.conf["fields"][c]["beam_chan"] == beam_chan:
                # Only reconstruct sun_beam if we have switched frequency
                pass
            else:
                beam_chan = self.conf["fields"][c]["beam_chan"]
                sun_beam = self.beams.interpolate_beam_2d(
                    self.beams.beam_str_to_idx(beam_chan),
                    solar_noon_gps,
                    all_obstimes,
                    None,
                    float(day["dec_sun"]),
                )
            # Next, the target field grid
            target_beam = self.beams.interpolate_beam_2d(
                self.beams.beam_str_to_idx(beam_chan),
                solar_noon_gps,
                all_obstimes,
                float(day["ha_%s" % c]),
                float(day["dec_%s" % (c)]),
            )

            sun_filter = sun_beam < 10 ** self.conf["solarAttenuationCutoff"]
            # applying sun_filter to target_beam will return a ravelled array.
            # we need to be able to identify the original location of our peak sensitivity within target_beam
            ha_grid, beam_grid = np.meshgrid(
                np.arange(target_beam.shape[1]), np.arange(target_beam.shape[0])
            )

            try:
                flat_idx = np.nanargmax(target_beam[sun_filter & flag_filter & target_filter])
            except ValueError:
                print(
                    "Warning, no observation meets criteria for %s day %s"
                    % (day["local_noon_str"], c)
                )
                continue
            ha_idx = ha_grid[sun_filter & flag_filter & target_filter][flat_idx]
            beam_idx = beam_grid[sun_filter & flag_filter & target_filter][flat_idx]
            obstime = all_obstimes[ha_idx]
            # FIXME: hard-coded indices
            target_filter[0, ha_idx - 75 : ha_idx + 75] = False
            flag_mask[ha_idx - 37 : ha_idx + 38] = 2

            out_dict["beam_%s" % c] = self.beams.df["beams"].dims[1][0][beam_idx]
            out_dict["ha_idx_%s" % c] = ha_idx
            # FIXME: hard-coded index
            out_dict["starttime_%s" % c] = obstime - 296
            out_dict["sun_attenuation_%s" % c] = np.log10(sun_beam[beam_idx, ha_idx])
            out_dict["target_sensitivity_%s" % c] = target_beam[beam_idx, ha_idx]

        for c in self.conf["priority"]:
            ha_idx = out_dict["ha_idx_%s" % c]
            if not np.isnan(ha_idx):
                # FIXME: hard-coded indices
                out_dict["unflagged_before_%s" % c] = num_unflagged(flag_mask, ha_idx - 38, False)
                out_dict["unflagged_after_%s" % c] = num_unflagged(
                    flag_mask.astype(bool), ha_idx + 38, True
                )
        self.observations.append(out_dict)

#class DayScheduler(Scheduler):
#    def __init__()
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
            self.schedule_day(solar_noon_gps, day, out_dict)

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
