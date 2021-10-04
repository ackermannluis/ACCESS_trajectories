#!/usr/bin/env python
# Copyright 2021
# author: Luis Ackermann <ackermann.luis@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from  U_Analysis_main import *

process_id = str(os.getpid())
################################################################ user defined variables #################
trajectories_start_data_file_str_lookup = 'ACCESS_trajectories_start_data'
path_output = '/scratch/k10/la6753/tmp/ACCESS_trajectories/'
path_log = '/g/data/k10/la6753/job_logs/'
log_msg_main_label = 'ACCESS_traj_calc_' + str(process_id)

processes_ = 36 # works good in a "normal" qsub, request 37 CPUs such that there is always one available for scheduling
################################################################ user defined variables #################
valid_time_start_stop_row = [2,8] # removes spin-up time and far predictions (where newer predictions are available)
path_data = '/g/data/lb4/ops_aps2/access-vt/1/'
wind_filename_U = 'wnd_ucmp.nc'
wind_filename_V = 'wnd_vcmp.nc'
wind_filename_W = 'vertical_wnd.nc'
height_filename = 'geop_ht.nc'
sfc_height_filename = path_data + '20180320/0000/fc/sfc/topog.nc'
lvl_filename = path_data + '20180320/0000/fc/pl/geop_ht.nc'
model_output_type = 'fc/pl/'
time_units_format = 'seconds since %Y-%m-%d %H'
path_program = os.path.dirname(os.path.realpath(sys.argv[0])) + '/'


# get lat, lon, sfc_height
topog_file = nc.Dataset(sfc_height_filename)
lat_    = topog_file.variables['lat'][:].filled(np.nan)
lon_    = topog_file.variables['lon'][:].filled(np.nan)
topog_  = topog_file.variables['topog'][0,:,:].filled(np.nan)
topog_file.close()

# get lvl
lvl_file = nc.Dataset(lvl_filename)
lvl_    = lvl_file.variables['lvl'][:].filled(np.nan)
lvl_file.close()

# calculate median grid distance
lat_delta_mean = np.mean(np.diff(lat_, axis=0))
lon_delta_mean = np.mean(np.diff(lon_, axis=0))
m_per_deg_lat, m_per_deg_lon = meter_per_degrees(np.mean(lat_))
lat_delta_mean_meters = np.abs(lat_delta_mean * m_per_deg_lat)
lon_delta_mean_meters = np.abs(lon_delta_mean * m_per_deg_lon)
mean_grid_size = np.mean([lat_delta_mean_meters, lon_delta_mean_meters])


# list ACCESS files and create file time series
ACCESS_U_file_list_full = sorted(list_files_recursive(path_data, wind_filename_U))
ACCESS_V_file_list_full = sorted(list_files_recursive(path_data, wind_filename_V))
ACCESS_W_file_list_full = sorted(list_files_recursive(path_data, wind_filename_W))
ACCESS_H_file_list_full = sorted(list_files_recursive(path_data, height_filename))

# remove filenames that are not from the model_output_type path
ACCESS_U_file_list = []
ACCESS_V_file_list = []
ACCESS_W_file_list = []
ACCESS_H_file_list = []
for filename_ in ACCESS_U_file_list_full:
    if model_output_type + wind_filename_U in filename_:
        ACCESS_U_file_list.append(filename_)
for filename_ in ACCESS_V_file_list_full:
    if model_output_type + wind_filename_V in filename_:
        ACCESS_V_file_list.append(filename_)
for filename_ in ACCESS_W_file_list_full:
    if model_output_type + wind_filename_W in filename_:
        ACCESS_W_file_list.append(filename_)
for filename_ in ACCESS_H_file_list_full:
    if model_output_type + height_filename in filename_:
        ACCESS_H_file_list.append(filename_)
ACCESS_U_file_list = sorted(ACCESS_U_file_list)
ACCESS_V_file_list = sorted(ACCESS_V_file_list)
ACCESS_W_file_list = sorted(ACCESS_W_file_list)
ACCESS_H_file_list = sorted(ACCESS_H_file_list)

file_list_times_sec = []
for filename_ in ACCESS_U_file_list:
    file_list_times_sec.append(time_str_to_seconds(filename_[len(path_data):len(path_data) + 13], '%Y%m%d/%H%M'))
file_list_times_sec = np.array(file_list_times_sec) + (valid_time_start_stop_row[0] * 3600)


def particle_trajectory_from_ACCESS(arg_list):
    global lat_, lon_, lvl_, topog_, mean_grid_size, \
        lat_delta_mean_meters, lon_delta_mean_meters, m_per_deg_lat, m_per_deg_lon, \
        ACCESS_U_file_list, ACCESS_V_file_list, ACCESS_W_file_list, ACCESS_H_file_list, file_list_times_sec

    start_time_YYYYmmDDHHMM_str,hours_int,start_point_lat,start_point_lon,start_point_height_m_ASL,\
    output_time_average_sec, array_instead_of_dict = arg_list

    # convert args to right data type
    hours_int = int(hours_int)
    start_point_lat = float(start_point_lat)
    start_point_lon = float(start_point_lon)
    start_point_height_m_ASL = float(start_point_height_m_ASL)
    output_time_average_sec = int(output_time_average_sec)
    if 't' in array_instead_of_dict or 'T' in array_instead_of_dict:
        array_instead_of_dict = True
    else:
        array_instead_of_dict = False

    # create back trajectory from wrf output
    start_time_sec = time_str_to_seconds(start_time_YYYYmmDDHHMM_str, '%Y%m%d%H%M')
    start_time_sec_original = start_time_sec * 1
    point_lat_lon = [start_point_lat, start_point_lon]
    point_height = start_point_height_m_ASL
    backtraj_hours = hours_int
    backtraj_secs = np.abs(backtraj_hours * 60 * 60)
    if hours_int < 0:
        traj_direction = -1
    elif hours_int > 0:
        traj_direction = 1
    stop_time_sec = start_time_sec + backtraj_secs * traj_direction


    # start output lists
    output_list_time = [start_time_sec]
    output_list_lat = [point_lat_lon[0]]
    output_list_lon = [point_lat_lon[1]]
    output_list_heigh = [point_height]


    if traj_direction == -1:
        while start_time_sec > stop_time_sec:

            # find spatial index of point
            point_r, point_c = find_index_from_lat_lon(lat_, lon_, point_lat_lon[0], point_lat_lon[1])

            # check if point is at edge of domain
            if point_r <= 3 or point_r >= lat_.shape[0] - 3 or point_c <= 3 or point_c >= lon_.shape[0] - 3:
                break


            # find starting files index
            file_index_1 = time_to_row_sec(file_list_times_sec, start_time_sec)
            # find starting time row index within file
            file_1_times = np.arange(file_list_times_sec[file_index_1],
                                     file_list_times_sec[file_index_1] +
                                     (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600, 3600)
            file_index_1_row = time_to_row_sec(file_1_times, start_time_sec)
            if file_1_times[file_index_1_row] > start_time_sec: # make the time row be to the left
                file_index_1_row -= 1
                if file_index_1_row < 0: # in case moving it to the left puts it on the previous file
                    file_index_1 -= 1
                    file_1_times = np.arange(file_list_times_sec[file_index_1],
                                             file_list_times_sec[file_index_1] +
                                             (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600, 3600)
                    file_index_1_row = file_1_times.shape[0]-1

                    # find the filename row index of the other time stamp
                    file_index_2 = file_index_1 + 1
                    file_index_2_row = 0
                    file_2_times = np.arange(file_list_times_sec[file_index_2],
                                             file_list_times_sec[file_index_2] +
                                             (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600,
                                             3600)

            elif file_index_1_row >= len(file_1_times)-1:  # needs the next file's first valid time row
                file_index_2 = file_index_1 + 1
                file_index_2_row = 0
                file_2_times = np.arange(file_list_times_sec[file_index_2],
                                         file_list_times_sec[file_index_2] +
                                         (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600, 3600)
            else:  # same file, just next row
                file_index_2 = file_index_1
                file_index_2_row = file_index_1_row + 1
                file_2_times = file_1_times



                # check that files are still available
                if file_index_1 == 0 or file_index_2 == 0 or \
                        file_index_1 >= len(file_list_times_sec) or file_index_2 >= len(file_list_times_sec):
                    break

            # calculate temporal weight for each file
            file_index_1_row_time_sec = file_1_times[file_index_1_row]
            file_index_2_row_time_sec = file_2_times[file_index_2_row]
            Weight_T_1 = np.abs(file_index_1_row_time_sec - start_time_sec)
            Weight_T_2 = np.abs(file_index_2_row_time_sec - start_time_sec)
            W_T_sum = Weight_T_2 + Weight_T_1
            Weight_T_1 = 1 - (Weight_T_1 / W_T_sum)
            Weight_T_2 = 1 - (Weight_T_2 / W_T_sum)

            # check that file is actually close to start time (in case model run missing)
            if Weight_T_1 > 2 * 3600 or Weight_T_2 > 2 * 3600:
                break

            # find the 4 grids closest to the point, place point_r and point_c on the top-left of the point
            if lat_[point_r] - point_lat_lon[0] < 0:
                point_r += 1
            if lon_[point_c] - point_lat_lon[1] > 0:
                point_c -= 1

            # create distance array to closest 4 points (2 by 2 square)
            lat_2D_arr = np.zeros((2,2))
            lon_2D_arr = np.zeros((2,2))
            lat_2D_arr[:,0] = lat_[point_r:point_r + 2]
            lat_2D_arr[:,1] = lat_[point_r:point_r + 2]
            lon_2D_arr[0,:] = lon_[point_c:point_c + 2]
            lon_2D_arr[1,:] = lon_[point_c:point_c + 2]
            D_H = distance_array_lat_lon_2D_arrays_degress_to_meters(lat_2D_arr,lon_2D_arr,
                                                                     point_lat_lon[0], point_lat_lon[1])

            # get point height array
            GPH_file = nc.Dataset(ACCESS_H_file_list[file_index_1])
            Z_ = GPH_file.variables['geop_ht'][file_index_1_row,:,
                 point_r:point_r+2, point_c:point_c+2].filled(np.nan) / gravity_
            GPH_file.close()


            # find height index, place Z_index on the bottom of the point
            Z_square_mean = np.mean(np.mean(Z_, axis=-1), axis=-1)
            Z_index = find_min_index_1d_array(np.abs(Z_square_mean - point_height))
            if Z_square_mean[Z_index] - point_height > 0:
                Z_index -= 1
            if Z_square_mean[Z_index] < topog_[point_r,point_c] or Z_index < 0: # in case it is bellow surface...
                Z_index += 1
            if Z_index >= lvl_.shape[0]-1: # in case it is at model's top
                Z_index -= 1


            # create vertical distance array to closest 2 layers
            D_Z = np.abs(Z_[Z_index:Z_index + 2] - point_height)

            # create absolute distance array to closest 8 points
            D_ = np.zeros(D_Z.shape)
            if D_[0, :, :].shape != (2,2) or D_Z[0, :, :].shape != (2,2) or D_H.shape != (2,2) or  D_.shape[0]!=2:
                print('D_.shape=',D_.shape)
                print('D_Z.shape=', D_Z.shape)
                print('D_H.shape=', D_H.shape)
                print('Z_.shape=', Z_.shape)
                print('Z_index=', Z_index)
                print('point_r=', point_r)
                print('point_c=', point_c)
                print('point_lat_lon=', point_lat_lon)
                print('point_height=', point_height)
                print('start_time_sec=', start_time_sec)
            D_[0, :, :] = (D_Z[0, :, :] ** 2 + D_H ** 2) ** 0.5
            D_[1, :, :] = (D_Z[1, :, :] ** 2 + D_H ** 2) ** 0.5

            # check if point is exactly at some grid center, if yes set weight to 1, else, distribute weights (avoids error x/0)
            if np.sum(D_ == 0) == 0:
                # calculate model points' weights
                D_reciprocal = 1 / D_
                Weights_D = D_reciprocal / np.sum(D_reciprocal)
            else:
                # set model points' weights to 1 for exact grid where point happens to be
                Weights_D = D_ * 0
                Weights_D[D_ == 0] = 1

            ################ get 3D wind data ############
            # file 1
            wrf_nc_1 = nc.Dataset(ACCESS_U_file_list[file_index_1])
            U_1_arr = wrf_nc_1.variables['wnd_ucmp'][file_index_1_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_1.close()
            wrf_nc_1 = nc.Dataset(ACCESS_V_file_list[file_index_1])
            V_1_arr = wrf_nc_1.variables['wnd_vcmp'][file_index_1_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_1.close()
            wrf_nc_1 = nc.Dataset(ACCESS_W_file_list[file_index_1])
            W_1_arr = wrf_nc_1.variables['vertical_wnd'][file_index_1_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_1.close()

            # file 0
            wrf_nc_2 = nc.Dataset(ACCESS_U_file_list[file_index_2])
            U_2_arr = wrf_nc_2.variables['wnd_ucmp'][file_index_2_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_2.close()
            wrf_nc_2 = nc.Dataset(ACCESS_V_file_list[file_index_2])
            V_2_arr = wrf_nc_2.variables['wnd_vcmp'][file_index_2_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_2.close()
            wrf_nc_2 = nc.Dataset(ACCESS_W_file_list[file_index_2])
            W_2_arr = wrf_nc_2.variables['vertical_wnd'][file_index_2_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_2.close()

            # calculate weighed mean wind components
            U_1_weighed = (Weight_T_1 * np.sum(U_1_arr * Weights_D)) + (Weight_T_2 * np.sum(U_2_arr * Weights_D))
            U_ = U_1_weighed * traj_direction

            V_1_weighed = (Weight_T_1 * np.sum(V_1_arr * Weights_D)) + (Weight_T_2 * np.sum(V_2_arr * Weights_D))
            V_ = V_1_weighed * traj_direction

            W_1_weighed = (Weight_T_1 * np.sum(W_1_arr * Weights_D)) + (Weight_T_2 * np.sum(W_2_arr * Weights_D))
            W_ = W_1_weighed * traj_direction

            # calculate horizontal wind speed
            WS_ = (U_ ** 2 + V_ ** 2) ** 0.5

            # if gotten to a dead point with absolutely no wind, exit loop as it will never leave this point
            if WS_ <= 0:
                break

            # calculate time to move to next grid
            time_grid_change_sec = int(mean_grid_size / WS_)
            # set next time stamp
            start_time_sec = start_time_sec + (time_grid_change_sec * traj_direction)

            # calculate vertical displacement
            vertical_delta = W_ * time_grid_change_sec

            # calculate horizontal displacement
            lat_delta_meter = V_ * time_grid_change_sec
            lon_delta_meter = U_ * time_grid_change_sec

            # get new point lat lon
            deg_per_m_lat, deg_per_m_lon = degrees_per_meter(point_lat_lon[0])
            point_lat_lon = [point_lat_lon[0] + (deg_per_m_lat * lat_delta_meter),
                             point_lat_lon[1] + (deg_per_m_lon * lon_delta_meter)]

            # get new point height
            point_height = point_height + vertical_delta
            if point_height < topog_[point_r,point_c]: # in case it ran into ground
                point_height = topog_[point_r, point_c]

            # store data
            output_list_time.append(start_time_sec)
            output_list_lat.append(point_lat_lon[0])
            output_list_lon.append(point_lat_lon[1])
            output_list_heigh.append(point_height)

    else:
        while start_time_sec < stop_time_sec:


            # find spatial index of point
            point_r, point_c = find_index_from_lat_lon(lat_, lon_, point_lat_lon[0], point_lat_lon[1])

            # check if point is at edge of domain
            if point_r == 0 or point_r > lat_.shape[0] - 2 or point_c == 0 or point_c > lon_.shape[0] - 2:
                break

            # find starting files index
            file_index_1 = time_to_row_sec(file_list_times_sec, start_time_sec)
            # find starting time row index within file
            file_1_times = np.arange(file_list_times_sec[file_index_1],
                                     file_list_times_sec[file_index_1] +
                                     (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600, 3600)
            file_index_1_row = time_to_row_sec(file_1_times, start_time_sec)
            if file_1_times[file_index_1_row] > start_time_sec: # make the time row be to the left
                file_index_1_row -= 1
                if file_index_1_row < 0: # in case moving it to the left puts it on the previous file
                    file_index_1 -= 1
                    file_1_times = np.arange(file_list_times_sec[file_index_1],
                                             file_list_times_sec[file_index_1] +
                                             (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600, 3600)
                    file_index_1_row = file_1_times.shape[0]-1

                    # find the filename row index of the other time stamp
                    file_index_2 = file_index_1 + 1
                    file_index_2_row = 0
                    file_2_times = np.arange(file_list_times_sec[file_index_2],
                                             file_list_times_sec[file_index_2] +
                                             (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600,
                                             3600)

            elif file_index_1_row >= len(file_1_times)-1:  # needs the next file's first valid time row
                file_index_2 = file_index_1 + 1
                file_index_2_row = 0
                file_2_times = np.arange(file_list_times_sec[file_index_2],
                                         file_list_times_sec[file_index_2] +
                                         (valid_time_start_stop_row[1] - valid_time_start_stop_row[0]) * 3600, 3600)
            else:  # same file, just next row
                file_index_2 = file_index_1
                file_index_2_row = file_index_1_row + 1
                file_2_times = file_1_times

            # calculate temporal weight for each file
            file_index_1_row_time_sec = file_1_times[file_index_1_row]
            file_index_2_row_time_sec = file_2_times[file_index_2_row]
            Weight_T_1 = np.abs(file_index_1_row_time_sec - start_time_sec)
            Weight_T_2 = np.abs(file_index_2_row_time_sec - start_time_sec)
            W_T_sum = Weight_T_2 + Weight_T_1
            Weight_T_1 = 1 - (Weight_T_1 / W_T_sum)
            Weight_T_2 = 1 - (Weight_T_2 / W_T_sum)

            # check that file is actually close to start time (in case model run missing)
            if Weight_T_1 > 2 * 3600 or Weight_T_2 > 2 * 3600:
                break

            # find the 4 grids closest to the point, place point_r and point_c on the top-left of the point
            if lat_[point_r] - point_lat_lon[0] < 0:
                point_r += 1
            if lon_[point_c] - point_lat_lon[1] > 0:
                point_c -= 1

            # create distance array to closest 4 points (2 by 2 square)
            lat_2D_arr = np.zeros((2,2))
            lon_2D_arr = np.zeros((2,2))
            lat_2D_arr[:,0] = lat_[point_r:point_r + 2]
            lat_2D_arr[:,1] = lat_[point_r:point_r + 2]
            lon_2D_arr[0,:] = lon_[point_c:point_c + 2]
            lon_2D_arr[1,:] = lon_[point_c:point_c + 2]
            D_H = distance_array_lat_lon_2D_arrays_degress_to_meters(lat_2D_arr,lon_2D_arr,
                                                                     point_lat_lon[0], point_lat_lon[1])

            # get point height array
            GPH_file = nc.Dataset(ACCESS_H_file_list[file_index_1])
            Z_ = GPH_file.variables['geop_ht'][file_index_1_row,:,
                 point_r:point_r+2, point_c:point_c+2].filled(np.nan) / gravity_
            GPH_file.close()


            # find height index, place Z_index on the bottom of the point
            Z_square_mean = np.mean(np.mean(Z_, axis=-1), axis=-1)
            Z_index = find_min_index_1d_array(np.abs(Z_square_mean - point_height))
            if Z_square_mean[Z_index] - point_height > 0:
                Z_index -= 1
            if Z_square_mean[Z_index] < topog_[point_r,point_c] or Z_index < 0: # in case it is bellow surface...
                Z_index += 1
            if Z_index >= lvl_.shape[0]-1: # in case it is at model's top
                Z_index -= 1


            # create vertical distance array to closest 2 layers
            D_Z = np.abs(Z_[Z_index:Z_index + 2] - point_height)

            # create absolute distance array to closest 8 points
            D_ = np.zeros(D_Z.shape)
            D_[0, :, :] = (D_Z[0, :, :] ** 2 + D_H ** 2) ** 0.5
            D_[1, :, :] = (D_Z[1, :, :] ** 2 + D_H ** 2) ** 0.5

            # check if point is exactly at some grid center, if yes set weight to 1, else, distribute weights (avoids error x/0)
            if np.sum(D_ == 0) == 0:
                # calculate model points' weights
                D_reciprocal = 1 / D_
                Weights_D = D_reciprocal / np.sum(D_reciprocal)
            else:
                # set model points' weights to 1 for exact grid where point happens to be
                Weights_D = D_ * 0
                Weights_D[D_ == 0] = 1

            ################ get 3D wind data ############
            # file 1
            wrf_nc_1 = nc.Dataset(ACCESS_U_file_list[file_index_1])
            U_1_arr = wrf_nc_1.variables['wnd_ucmp'][file_index_1_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_1.close()
            wrf_nc_1 = nc.Dataset(ACCESS_V_file_list[file_index_1])
            V_1_arr = wrf_nc_1.variables['wnd_vcmp'][file_index_1_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_1.close()
            wrf_nc_1 = nc.Dataset(ACCESS_W_file_list[file_index_1])
            W_1_arr = wrf_nc_1.variables['vertical_wnd'][file_index_1_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_1.close()

            # file 0
            wrf_nc_2 = nc.Dataset(ACCESS_U_file_list[file_index_2])
            U_2_arr = wrf_nc_2.variables['wnd_ucmp'][file_index_2_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_2.close()
            wrf_nc_2 = nc.Dataset(ACCESS_V_file_list[file_index_2])
            V_2_arr = wrf_nc_2.variables['wnd_vcmp'][file_index_2_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_2.close()
            wrf_nc_2 = nc.Dataset(ACCESS_W_file_list[file_index_2])
            W_2_arr = wrf_nc_2.variables['vertical_wnd'][file_index_2_row + valid_time_start_stop_row[0],
                      Z_index:Z_index + 2, point_r:point_r + 2, point_c:point_c + 2].data
            wrf_nc_2.close()

            # calculate weighed mean wind components
            U_1_weighed = (Weight_T_1 * np.sum(U_1_arr * Weights_D)) + (Weight_T_2 * np.sum(U_2_arr * Weights_D))
            U_ = U_1_weighed * traj_direction

            V_1_weighed = (Weight_T_1 * np.sum(V_1_arr * Weights_D)) + (Weight_T_2 * np.sum(V_2_arr * Weights_D))
            V_ = V_1_weighed * traj_direction

            W_1_weighed = (Weight_T_1 * np.sum(W_1_arr * Weights_D)) + (Weight_T_2 * np.sum(W_2_arr * Weights_D))
            W_ = W_1_weighed * traj_direction

            # calculate horizontal wind speed
            WS_ = (U_ ** 2 + V_ ** 2) ** 0.5

            # if gotten to a dead point with absolutely no wind, exit loop as it will never leave this point
            if WS_ <= 0:
                break

            # calculate time to move to next grid
            time_grid_change_sec = int(mean_grid_size / WS_)
            # set next time stamp
            start_time_sec = start_time_sec + (time_grid_change_sec * traj_direction)

            # calculate vertical displacement
            vertical_delta = W_ * time_grid_change_sec

            # calculate horizontal displacement
            lat_delta_meter = V_ * time_grid_change_sec
            lon_delta_meter = U_ * time_grid_change_sec

            # get new point lat lon
            deg_per_m_lat, deg_per_m_lon = degrees_per_meter(point_lat_lon[0])
            point_lat_lon = [point_lat_lon[0] + (deg_per_m_lat * lat_delta_meter),
                             point_lat_lon[1] + (deg_per_m_lon * lon_delta_meter)]

            # get new point height
            point_height = point_height + vertical_delta
            if point_height < topog_[point_r,point_c]: # in case it ran into ground
                point_height = topog_[point_r, point_c]

            # store data
            output_list_time.append(start_time_sec)
            output_list_lat.append(point_lat_lon[0])
            output_list_lon.append(point_lat_lon[1])
            output_list_heigh.append(point_height)


    # convert data to array
    output_ = np.column_stack((
        np.array(output_list_time),
        np.array(output_list_lat),
        np.array(output_list_lon),
        np.array(output_list_heigh)
    ))

    # average time wise?
    if output_time_average_sec is not None:
        if traj_direction == -1:
           time_mean, vals_mean = mean_discrete(output_[:,0], output_[:,1:], output_time_average_sec,
                                                stop_time_sec, last_index=start_time_sec_original,show_progress=False)
        else:
            time_mean, vals_mean = mean_discrete(output_[:, 0], output_[:, 1:], output_time_average_sec,
                                                 start_time_sec_original, last_index=stop_time_sec,show_progress=False)
        output_ = np.column_stack((time_mean, vals_mean))

    # convert lists to output dictionary
    if array_instead_of_dict:
        pass
    else:
        output_ = {
            'time'   : output_[:,0],
            'lat'    : output_[:,1],
            'lon'    : output_[:,2],
            'height' : output_[:,3]
        }

    # save data to local node
    filename_ = path_program + start_time_YYYYmmDDHHMM_str + '_' + str(hours_int) + '_' + str(start_point_lat) + \
                '_' + str(start_point_lon) + '_' + str(start_point_height_m_ASL) + '_'
    np.save(filename_, output_)



if __name__ == '__main__':
    log_msg('main process started', log_msg_main_label)

    # count processors and define pool
    from multiprocessing import Pool
    Pool_ = Pool(processes_)
    log_msg('started children processes', log_msg_main_label)


    # get trajectories arguments
    trajectories_start_data_file = list_files_recursive(path_program, trajectories_start_data_file_str_lookup)[0]
    trajectories_start_arr = open_csv_file(trajectories_start_data_file, delimiter=',', skip_head=1, dtype='<U128')
    # create list of arguments to pass to pool
    args_list = []
    for row_index in range(trajectories_start_arr.shape[0]):
        args_list.append(list(trajectories_start_arr[row_index,:]))
    log_msg('created argument list', log_msg_main_label)


    # start calculating trajectories
    if len(args_list) > 100:
        chunk_size = processes_ * 4
        chunks_ = math.ceil(len(args_list)/chunk_size)
        for chunk_index in range(chunks_):
            multip_out = Pool_.map(particle_trajectory_from_ACCESS,
                                   args_list[chunk_index*chunk_size:(chunk_index+1)*chunk_size])
            log_msg('finished calculations for chunk ' + str(chunk_index) + ' of ' + str(chunks_), log_msg_main_label)

    else:
        multip_out = Pool_.map(particle_trajectory_from_ACCESS, args_list)
        # report possible errors
        log_msg('finished calculations', log_msg_main_label)


    # get file list
    file_list_clean = list_files_recursive(path_program, filter_str='_.npy')
    # start the zip file
    current_time_struct = time.gmtime()
    current_time_str = str(current_time_struct[0]).zfill(4) + \
                       str(current_time_struct[1]).zfill(2) + \
                       str(current_time_struct[2]).zfill(2) + '_' + \
                       str(current_time_struct[3]).zfill(2) + \
                       str(current_time_struct[4]).zfill(2) + \
                       str(current_time_struct[5]).zfill(2)
    filename_zip = path_output + 'ACCESS_trajectories_' + trajectories_start_data_file[-6:-4] + '_' + \
                   current_time_str + '_' + str(process_id) + '.zip'
    zf = zipfile.ZipFile(filename_zip, mode='w')
    zf.close()
    # zip pngs and delete source
    for filename_ in file_list_clean:
        zf = zipfile.ZipFile(filename_zip, mode='a')
        zf.write(filename_)
        zf.close()
    log_msg('compressed clean files!', log_msg_main_label)












