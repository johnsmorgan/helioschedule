import argparse
import requests
import json
import csv
from astropy.time import Time
from time import sleep
from yaml import safe_load
from helioschedule.schedule import get_min_max

LEEWAY = 3600  # s otherwise mwa metadata searches will miss observations that start or stop *during* search period
SLEEP = 1  # only send max 1 request per second to avoid overloading server


def get_flags(noons, offset, warn_obscodes=["G0060"], skip_obscodes=[]):
    flags = []

    for t, time in enumerate(noons):
        solar_noon_gps = int(round(float(time["local_noon_gps"])))
        mintime, maxtime = get_min_max(solar_noon_gps, offset=offset)

        result = requests.get(
            "https://ws.mwatelescope.org/metadata/find",
            data={"mintime": f"{mintime-LEEWAY}", "maxtime": f"{maxtime+LEEWAY}", "extended": 1},
        )
        if result.text == "":
            print("no observations -- skipping")
            continue
        try:
            result_dict = json.loads(result.text)
        except json.JSONDecodeError:
            if "No observations found." in result.text:
                continue
            print("can't decode the following result--expecting valid json")
            print(result.text)
            raise

        start_times_unfiltered = [Time(o[0], format="gps").utc.isot for o in result_dict]
        print(f"{len(start_times_unfiltered)} unfiltered observations")
        print(set([o[4] for o in result_dict]))

        obs = [
            (o[0], o[1])
            for o in result_dict
            if o[1] > mintime
            if o[0] < maxtime
            if not o[4] in skip_obscodes
        ]
        obscodes = [
            o[4]
            for o in result_dict
            if o[1] > mintime
            if o[0] < maxtime
            if not o[4] in skip_obscodes
        ]
        for warn_obscode in warn_obscodes:
            if warn_obscode in obscodes:
                print(f"Warning! {warn_obscode} observations already scheduled!")
        duration = sum([o[1] - o[0] for o in obs]) / 60.0
        print(f"{len(obs)} unfiltered observations totalling {duration:.2f} minutes")
        print()
        flags = flags + obs
        sleep(SLEEP)
    return flags


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input yaml file")
    parser.add_argument(
        "--exclude_obscodes",
        help="Don't generate flags for these obscodes (comma-separated list)",
        default="",
    )
    args = parser.parse_args()
    conf = safe_load(open(args.infile))

    # allow flagging to be disabled
    if conf["files"]["flags"] == "":
        f = open(conf["files"]["flags"], "w")
        f.close()
        return

    noons = list(csv.DictReader(line for line in open(conf["files"]["noons"])))
    flags = get_flags(noons, conf["solarOffset"], skip_obscodes=args.exclude_obscodes.split(","))
    if not flags == []:
        with open(conf["files"]["flags"], "w") as f:
            json.dump(flags, f)


if __name__ == "__main__":
    main()
