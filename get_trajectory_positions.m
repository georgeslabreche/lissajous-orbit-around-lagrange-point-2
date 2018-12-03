function [X, Y, Z] = get_trajectory_positions(data_filename,...
    is_earth_origin,...
    include_earth_to_L2_orbit_trajectory,...
    include_decommissioning_trajectory)
    %GET_TRAJECTORY_POSITIONS Reads the given data file and returns X, Y, Z
    %coordinate arrays

    % Generate data table.
    date_format = 'dd MMMM yyyy HH:mm:ss.SSS';
    T = readtable(data_filename, 'Format',strcat('%{', date_format, '}D %f %f %f'));

    % Create index column for ease of referencing/debugging.
    Indices = [1:1:height(T)].';
    T = addvars(T,Indices,'Before',1);

    % Rename table columns.
    T.Properties.VariableNames = {'Index', 'Timestamp' 'X' 'Y' 'Z'};
    
    % Get key date times.
    orbit_insertion_datetime = '02 March 2019 02:40:15.978';
    orbit_decommissioning_datetime_str = '06 June 2020 20:14:06.891';

    date_strings = {orbit_insertion_datetime; orbit_decommissioning_datetime_str};
    dts = datetime(date_strings, 'InputFormat', date_format);

    % Scope the trajectory as desired, 
    orbit_start_index = 1;
    if include_earth_to_L2_orbit_trajectory == false
        orbit_start_index = T.Index(T.Timestamp == dts(1));
    end

    orbit_end_index = height(T);
    if include_decommissioning_trajectory == false
        orbit_end_index = T.Index(T.Timestamp == dts(2));
    end

    % Build X, Y, Z trajactor coordinate arrays.
    X = T(orbit_start_index:orbit_end_index, 3).X;
    Y = T(orbit_start_index:orbit_end_index, 4).Y;
    Z = T(orbit_start_index:orbit_end_index, 5).Z;

    % Shift coordinates if we want X, Y, Z origin to be Earth.
    % Shift by Earth's coordinate value which corresponds to the first
    % row of the data table.
    if is_earth_origin
        X = X - T(1, 3).X;
        Y = Y - T(1, 4).Y;
        Z = Z - T(1, 5).Z;
    end

end

