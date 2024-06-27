import csv
import json
import argparse
import numpy as np
from yaml import safe_load
from astropy.io import ascii
from astropy import units as u
from astropy.time import Time

S = slice(None, None, None)

SUN_OBS_STR = "schedule_observation.py --starttime={pre_time_comma} --stoptime=++16s --freq='{obs_chan}' --obsname={obs_name_prefix}Sun --source=Sun --mode=MWAX_CORRELATOR --inttime={inttime} --freqres={freqres} --creator={creator} --project={project}"
SUN_OBS_STR_POST = "schedule_observation.py --starttime={post_time_comma} --stoptime=++16s --freq='{obs_chan}' --obsname={obs_name_prefix}Sun --source=Sun --mode=MWAX_CORRELATOR --inttime={inttime} --freqres={freqres} --creator={creator} --project={project}"
OBSERVATION_STR = "schedule_observation.py --starttime={time_comma} --stoptime=++{duration}s --freq='{obs_chan}' --obsname={obs_name_prefix}{field} --shifttime={shifttime} --mode=MWAX_CORRELATOR --inttime={inttime} --freqres={freqres} --creator={creator} --project={project} --azimuth={az} --elevation={el}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    NO_WRITE = []

    azel = {int(k): v for k, v in json.load(open(conf["files"]["pointings"])).items()}

    obs_ha = ascii.read(conf["files"]["observations"])
    noons = Time(obs_ha["local_noon_str"][S])

    observations = []

    target_time = 0
    sun_target_time = 0
    schedule_time = 0

    for t in conf["priority"]:
        if "starttime_%s" % t in obs_ha.colnames:
            times_starttime = Time(
                np.where(np.isnan(obs_ha["starttime_%s" % t]), 0, obs_ha["starttime_%s" % t]),
                format="gps",
            )
        # times = (
        #    noons
        #    + Angle(np.nan_to_num(obs_ha["ha_%s" % t][S]) * u.deg).cycle * u.sday
        #    - u.second * conf["obs"]["duration"] / 2.0
        # )
        # round to nearest 8s (MWA observations must start and stop on these boundaries)
        # dt = (86400 * times.jd2 % 8) * u.second
        # times = np.where(dt < 4 * u.s, times - dt, times + (8 * u.s - dt))
        if "starttime_%s" % t in obs_ha.colnames:
            # print('using starttime column, min mean max difference',
            # f"{np.nanmin(np.abs(times_starttime.gps-times)):.2f}",
            # f"{np.nanmean(np.abs(times_starttime.gps-times)):.2f}",
            # f"{np.nanmax(np.abs(times_starttime.gps-times)):.2f}",
            # )
            times = times_starttime
        else:
            raise RuntimeError("missing start times, does not work with this version!")

        for j in range(len(noons)):
            if np.isnan(obs_ha["ha_idx_%s" % t][S][j]):
                continue
            out_dict = []
            out_dict = conf["obs"]
            for k in ('beam_chan', 'obs_chan', 'pre_time', 'duration', 'shifttime', 'inttime', 'freqres', 'creator', 'project', 'obs_name_prefix'):
                if k in conf["fields"][t].keys():
                    out_dict[k] = conf["fields"][t][k]
            out_dict["sweetspot"] = obs_ha["beam_%s" % t][S].data[j]
            out_dict["time"] = times[j].utc.isot[:19]
            out_dict["obsid"] = int(times[j].gps)
            assert out_dict["obsid"] == obs_ha["starttime_%s" % t][j], "Error in gps conversion!"
            out_dict["time_comma"] = times[j].utc.isot[:19].replace("T", ",")
            out_dict["pre_time_comma"] = (times[j] - 16 * u.second).utc.isot[:19].replace("T", ",")
            out_dict["post_time_comma"] = (
                (times[j] + conf["obs"]["duration"] * u.second).utc.isot[:19].replace("T", ",")
            )
            out_dict["az"], out_dict["el"] = azel[out_dict["sweetspot"]]
            out_dict["field"] = t
            out_dict["sun_attenuation"] = obs_ha["sun_attenuation_%s" % t][S].data[j]
            out_dict["target_sensitivity"] = obs_ha["target_sensitivity_%s" % t][S].data[j]
            out_dict["unflagged_before"] = obs_ha["unflagged_before_%s" % t][S].data[j]
            out_dict["unflagged_after"] = obs_ha["unflagged_after_%s" % t][S].data[j]
            observations.append(out_dict.copy())

    observations = sorted(observations, key=lambda o: o["time"])

    t1 = None
    with open(conf["files"]["schedule"], "w") as outfile:
        for o in observations:
            if o["field"] in NO_WRITE:
                t1 = Time(o["time"])
                continue
            t2 = Time(o["time"])
            if not t1 or (t1.datetime.day != t2.datetime.day):
                print("echo ", t2.isot[:10], file=outfile)
            # schedule Sun observations
            if o["unflagged_before"] >= 2:
                print(SUN_OBS_STR.format(**o), file=outfile)
                sun_target_time += 16
                schedule_time += 16
            elif o["unflagged_after"] >= 2:
                print(SUN_OBS_STR_POST.format(**o), file=outfile)
                sun_target_time += 16
                schedule_time += 16
            else:
                print(
                    "# no space around observation -- skipping Sun observation" % ((t2 - t1).sec),
                    file=outfile,
                )
            print(OBSERVATION_STR.format(**o), "#", o["obsid"], file=outfile)
            target_time += conf["obs"]["duration"]
            sun_target_time += conf["obs"]["duration"]
            schedule_time += conf["obs"]["duration"]
            t1 = t2

        print(
            "# target time: %d sun/target time: %d schedule time: %d"
            % (target_time, sun_target_time, schedule_time),
            file=outfile,
        )
    with open(conf["files"]["schedule_csv"], "w") as csvfile:
        fieldnames = (
            "time",
            "field",
            "sweetspot",
            "az",
            "el",
            "sun_attenuation",
            "target_sensitivity",
        )
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for out_dict in observations:
            writer.writerow(out_dict)


if __name__ == "__main__":
    main()
