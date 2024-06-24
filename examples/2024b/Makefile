# NB: Make is used to build targets from dependencies
# and track changes to depencies *only*
# yaml file controls filenames, and all other parameters
# it is up to the user to keep the yaml file and Makefile
# in sync!


schedule.sh : observations.csv
	generate_observation_script example.yaml
observations.csv : targets.csv
	schedule example.yaml
targets.csv : solar_noons.csv
	convert_coordinates example.yaml
solar_noons.csv :
	get_local_noons example.yaml
