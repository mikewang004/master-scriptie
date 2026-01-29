from analyse7 import *

last_timestep_e5 = polymer_properties(atom_coords("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt"))
last_timestep_e5.gyration_radius()

print(last_timestep_e5.results.gyration_radius)