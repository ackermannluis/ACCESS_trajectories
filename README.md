# ACCESS_trajectories
Calculates particle trajectories (back and forward in time) using The Australian Community Climate and Earth-System Simulator (ACCESS) weather model (version APS1). This code was created to work directly on Gadi (supercomputer at Australia's NCI - National Computational Infrastructure).

The script works with python 3. As it is set it requires to import the U_Analysis_main, which is a package with multiple atmospheric sciences related functions. Only a small subset of the functions in U_Analysis_main are used.

The script is set to calculate multiple trajectories in parallel. It has been tested using the normal queue for calculating 1,000 trajectories. Using 37 cpus it takes about 1.5 hours. The script cannot work on multiple nodes at a time, so only set cpus up to the maximum allowed per node. It is recommended that if trying to calculate a large number of trajectories (multiple thousands), then multiple 1000 jobs are requested, instead of increasing the cpu count.

A sample qsub file is included. The script reads the starting trajectories locations and times from a csv file called ACCESS_trajectories_start_data, where each line represents one trajectory start variables. A sample ACCESS_trajectories_start_data file is included. The qsub file will look for the ACCESS_trajectories_start_data and copy it to the node that will be used to calculate the trajectories, so the ACCESS_trajectories_start_data needs to be prepared and uploaded to the corresponding folder before submitting the job to the queue. A file called build.env is included that is run by the qsub file and defines the environment; this is the environment I use for all my scripts so it has more things that is needed for this project.

Some lines from the script need to be modified by the user (lines 20-25), this is to define the folder in Gadi where the output will be placed, the folder where the status of the calculation can be observed, and other similar variables. Here the number of child processes is defined; it is recommended that it is set to one less than the number of cpus requested in the qsub file.

The trajectories are calculated by using the horizontal wind speeds at the location of the particle and calculating how much time (dt) would it take the particle to travel the distance of the mean grid size (calculated from the model's grid). Then the actual horizontal (dx and dy) and vertical displacement (dz) is calculated from the horizontal and vertical wind speeds at the particle's location. To get these speeds, a weighted average is created from the closest 8 grid points from the model to the particle's location. The weights are based on the three dimensional distance between the particle and the grid centres, in other words, a point inside a box, where each corner are the model's grid centers. The temporal distance is also considered, such that the closest 2 model hours are used. These results in 16 elements being used to calculate the each wind component, the 8 closest model's grids around the particle for the closest simulated hour, and the same 8 model points for the second closest simulated hour. Weights based on the temporal distance are calculated and implemented on top of the spatial weights. Once the dx, dy, dz, and dy are calculated, the new particle location is defined, then the loop is restarted from this new location. Note that if the particle happens to be located right at a grid centre (or hour) then the weight for that grid (hour) will be set to one and the other grid point's (or hour) weights will be set to zero (will not be used).


Description of ACCESS_trajectories_start_data:
Header:
  start_time_YYYYmmDDHHMM_str,hours_int,start_point_lat,start_point_lon,start_point_height_m_ASL,output_time_average_sec,array_instead_of_dict

Each line after the header will be used to start a new trajectory. The data from all these trajectories will be saved to the node, and once all are calculated, they will be zipped and the zip file will be stored as per line 21 (path_output) of the main script.

-start_time_YYYYmmDDHHMM_str: start time of the trajectory in UTC
-hours_int: integer number of hours to calculate, if the particle exits the domain before this time, the trajectory will end. If negative, the trajectory will be backwards in time. 
-start_point_lat: latitude in degrees where the trajectory will be started, negative is south of the equator.
-start_point_lon: longitude in degrees where the trajectory will be started, between 0 and 360.
-start_point_height_m_ASL: altitude in meters above sea level where the trajectory will be started
-output_time_average_sec: if not None, then the trajectory will be averaged timewise to discrete intervals separated by these number of seconds.
-array_instead_of_dict: if True, then the data will be saved as a numpy array instead of a dictionary, with the first column as time in seconds since the epoch UTC, second column as latitude, third column as longitude, and last column as altitude (m). If false, it will be saved as a numpy dictionary (can be uploaded using load_dictionary(filename_) from U_Analysis_main).

