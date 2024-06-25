# This test Makefile is pretty standard, except that when
# generating the flag file, it excludes G0060 observations.
# This is because the test is aimed at exactly reproducing
# those observations.
# 
# NB: Make is used to build targets from dependencies
# and track changes to depencies *only*
# yaml file controls filenames, and all other parameters
# it is up to the user to keep the yaml file and Makefile
# in sync!

CONF_FILE=test.yaml

test: schedule.sh observations.csv
	$(info Diffing solar_noons.csv)
	diff solar_noons.csv solar_noons_reference.csv
	$(info Diffing flags.json)
	diff flags.json flags_reference.json
	$(info Diffing targets.csv)
	diff targets.csv targets_reference.csv
	$(info Diffing observations.csv)
	diff observations.csv observations_reference.csv
	$(info Diffing schedule.sh)
	diff schedule.sh schedule_reference.sh

clean:
	rm solar_noons.csv flags.json targets.csv observations.csv schedule.sh

schedule.sh : observations.csv $(CONF_FILE)
	generate_observation_script $(CONF_FILE)
observations.csv : targets.csv flags.json $(CONF_FILE)
	schedule $(CONF_FILE)
targets.csv : solar_noons.csv $(CONF_FILE)
	convert_coordinates $(CONF_FILE)
flags.json: solar_noons.csv $(CONF_FILE)
	get_flags $(CONF_FILE) --exclude_obscodes="G0060"
solar_noons.csv : $(CONF_FILE)
	get_local_noons $(CONF_FILE)
