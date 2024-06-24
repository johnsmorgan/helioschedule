"""
Optimized scheduling of observations
"""

import argparse
import csv
import numpy as np
from h5py import File
from scipy.interpolate import interp1d
from yaml import safe_load
from helioschedule.beam_interp import Beams
from astropy.time import Time

round_to_8 = lambda x: int(x - x%8)

SOLAR_OFFSET = 3072*8 # just over 6.8 hours
SOLAR_OFFSET/3600.
print

def get_min_max(solar_noon_gps, offset):
    mintime = round_to_8(solar_noon_gps-offset)
    maxtime = mintime + round_to_8(offset*2)
    return mintime, maxtime

def get_all_steps(solar_noon_gps, mintime, maxtime, solar_offset):
    return np.linspace(mintime, maxtime, solar_offset*2//8, endpoint=False)

def get_flags(flags, mintime, maxtime):
    """
    Select relevant flags
    flags - start/stop in gps seconds
    solar noon (gps sections)
    offset (seconds)
    """
    return [(o[0], o[1]) for o in flags if o[1]>mintime if o[0]<maxtime]

def flags_to_mask(flags, n_steps, mintime):
    """
    convert flags (in gpstime) to boolean mask
    """
    flag_mask = np.zeros(n_steps, dtype=bool)
    for f in flags:
        flag_mask[(f[0]-mintime)//8:(f[1]-mintime)//8] = True
    return flag_mask

def num_unflagged(array, idx, forward=True):
    """
    find the length of run of False after/before idx
    inclusive of idx
    """
    if array[idx]==True:
        return 0
    if forward is True:
        arr = array[idx:]
    else:
        arr=array[:idx][::-1]
    if np.all(arr==False):
        return len(arr)
    return np.where(arr != False)[0][0]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    targets = list(csv.DictReader(line for line in open(conf["files"]["targets"])))
    df = File(conf["files"]["beams"], "r")
    #FIXME dummy flags
    flags = [[0, 0] for t in range(len(targets))]
    # FIXME replacement has to deal with no flags on a particular day

    # deal with flags
    flag_img = np.zeros((len(targets), conf['solarOffset']*2//8), dtype=np.uint8)
    flag_img2 = flag_img.astype(int).copy()

    beams = Beams(conf["files"]["beams"])
    obs_ha = []
    k=0

    # loop over each observing day
    for t, target in enumerate(targets):
        print("scheduling around local noon", target['local_noon_str'])
        beam_chan = None
        out_dict = {}
        for key in ("local_noon_str", "local_noon_lst"):
            out_dict[key] = target[key]

        # Calculate and stop time based on solar noon 
        solar_noon_gps = int(round(float(target["local_noon_gps"])))
        mintime, maxtime = get_min_max(solar_noon_gps, conf['solarOffset'])
        all_obstimes = get_all_steps(solar_noon_gps, mintime, maxtime, conf['solarOffset'])
        
        obs_list = get_flags(flags, mintime, maxtime)
        flag_mask = flags_to_mask(obs_list, len(all_obstimes), mintime)
        #print(f"mask before convolution {np.sum(flag_mask)}")
        flag_mask_convolved = np.convolve(flag_mask, np.ones(75, dtype=bool), 'same')
        target_mask = np.zeros_like(flag_mask)
        target_filter = ~target_mask.reshape(1, -1)
        flag_filter = ~flag_mask_convolved.reshape(1, -1)

        # loop over fields for that day
        for c in conf["priority"]:
            if target["ha_%s" % c] == "":
                out_dict["ha_idx_%s" % c] = np.nan
                out_dict["beam_%s" % c] = -1
                out_dict["sun_attenuation_%s" % c] = np.nan
                out_dict["starttime_%s" % c] = np.nan
                out_dict["target_sensitivity_%s" % c] = np.nan
                out_dict["unflagged_before_%s" % c] = np.nan
                out_dict["unflagged_after_%s" % c] = np.nan
                continue
            if beam_chan is not None and conf["fields"][c]["beam_chan"] == beam_chan:
                # Only reconstruct sun_beam if we have switched frequency
                pass
            else:
                beam_chan = beam_chan
                sun_beam = beams.interpolate_beam_2d(beams.beam_str_to_idx(conf["fields"][c]["beam_chan"]),
                                            solar_noon_gps,
                                            all_obstimes,
                                            None,
                                            float(target["dec_sun"]))
            # Next, the target field grid
            target_beam = beams.interpolate_beam_2d(beams.beam_str_to_idx(conf["fields"][c]["beam_chan"]),
                                            solar_noon_gps,
                                            all_obstimes,
                                            float(target["ha_%s" % c]),
                                            float(target["dec_%s" % (c)]))

            sun_filter = sun_beam < 10 ** conf["solarAttenuationCutoff"]
            # applying sun_filter to target_beam will return a ravelled array.
            # we need to be able to identify the original location of our peak sensitivity within target_beam
            ha_grid, beam_grid = np.meshgrid(np.arange(target_beam.shape[1]), np.arange(target_beam.shape[0]))

            try:
                flat_idx = np.nanargmax(target_beam[sun_filter & target_filter])
            except ValueError:
                print("Warning, no observation meets criteria for %s target %s" % (target["local_noon_str"], c))
                out_dict["ha_idx_%s" % c] = np.nan
                out_dict["beam_%s" % c] = -1
                out_dict["sun_attenuation_%s" % c] = np.nan
                out_dict["starttime_%s" % c] = np.nan
                out_dict["target_sensitivity_%s" % c] = np.nan
                out_dict["unflagged_before_%s" % c] = np.nan
                out_dict["unflagged_after_%s" % c] = np.nan
                break
            ha_idx = ha_grid[sun_filter & flag_filter & target_filter][flat_idx]
            beam_idx = beam_grid[sun_filter & flag_filter & target_filter][flat_idx]
            obstime = all_obstimes[ha_idx]
            #print(f"sum target filter: {np.sum(target_filter)}")
            target_filter[0, ha_idx-75:ha_idx+75] = False
            flag_img[t, ha_idx-37:ha_idx+38] = 2
            flag_img2[t, ha_idx-37:ha_idx+38] = 2*k
            #print(f"sum target filter: {np.sum(target_filter)}")
            #print()
            for i, o in enumerate(obs_list):
                if o[0]>obstime:
                    obs_list.insert(i, [int(obstime)-296, int(obstime)+304, c])
                    break
            
            out_dict["beam_%s" % c] = beams.df["beams"].dims[1][0][beam_idx]
            out_dict["ha_idx_%s" % c] = ha_idx
            out_dict["starttime_%s" % c] = obstime-296
            out_dict["sun_attenuation_%s" % c] = np.log10(sun_beam[beam_idx, ha_idx])
            out_dict["target_sensitivity_%s" % c] = target_beam[beam_idx, ha_idx]
        for j, c in enumerate(conf["priority"]):
            ha_idx = out_dict['ha_idx_%s' % c]
            if not np.isnan(ha_idx):
                out_dict["unflagged_before_%s" % c] = num_unflagged(flag_img[t].astype(bool), ha_idx-38, False)
                out_dict["unflagged_after_%s" % c] = num_unflagged(flag_img[t].astype(bool), ha_idx+38, True)
        obs_ha.append(out_dict)
        k+=1

    with open(conf["files"]["observations"], "w") as csvfile:
        fieldnames = sorted(obs_ha[0].keys(), key=lambda k: k[::-1])
        fieldnames.insert(0, fieldnames.pop(-1))
        fieldnames.insert(0, fieldnames.pop(-1))
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')

        writer.writeheader()
        for out_dict in obs_ha:
            writer.writerow(out_dict)

if __name__ == "__main__":
    main()