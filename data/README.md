# helioschedule data

This directory contains the following file
`azel.json`
which maps each beamformer setting to its nominal azimuth and elevation setting.

Additionally, a full lookup table of beams in the appropriate format for each frequency setup is required. This can be calculated *from* a beam that has been pre-calculated using `mwa_pb_lookup` using the script
`generate_beams.py`

Note that `mwa_pb_lookup` uses beams that have been generated for a set of standard frequencies marking the edges of GLEAM band. `generate_beams.py` weights these frequencies appropriately in order to approximate the beam averaged across the full field of view.

Note that the chosen channel strings (which denote the band edges) *must* align with the GLEAM bands. `121-132` will work, `122-133` will not.
