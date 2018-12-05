function [Earth_Position, L2_Position, insertion_distance_from_L2, Duration, X, Y, Z] =...
    get_trajectory_positions(data_filename,...
    is_earth_origin,...
    include_earth_to_L2_orbit_trajectory,...
    include_decommissioning_trajectory)
    %GET_TRAJECTORY_POSITIONS Reads the given data file and returns X, Y, Z
    %coordinate arrays

    % Generate data table.
    date_format = 'dd MMMM yyyy HH:mm:ss.SSS';
    T = readtable(data_filename, 'Format',strcat('%{', date_format, '}D %f %f %f %f'));
 
    % Create index column for ease of referencing/debugging.
    Indices = [1:1:height(T)].';
    T = addvars(T,Indices,'Before',1);
    
    % Rename table columns.
    T.Properties.VariableNames = {'Index', 'Timestamp', 'Duration', 'X' 'Y' 'Z'};
       
    % Set Earth and L2 positions.
    Earth_Position = [T(1, 4).X, T(1, 5).Y, T(1, 6).Z];
    L2_Position = [0 0 0];
    
    % Get key date times.
    orbit_insertion_datetime = '02 March 2019 02:40:15.978';
    orbit_decommissioning_datetime_str = '06 June 2020 20:14:06.891';

    date_strings = {orbit_insertion_datetime; orbit_decommissioning_datetime_str};
    dts = datetime(date_strings, 'InputFormat', date_format);
    
    % Get insertion distance from L2 at initial position of the orbit:
    L2_orbit_insertion_position = [T.X(T.Timestamp == dts(1)) T.Y(T.Timestamp == dts(1)) T.Z(T.Timestamp == dts(1))];
    insertion_distance_from_L2 = norm(L2_orbit_insertion_position);

    % Scope the trajectory as desired, 
    orbit_start_index = 1;
    if include_earth_to_L2_orbit_trajectory == false
        orbit_start_index = T.Index(T.Timestamp == dts(1));
    end

    orbit_end_index = height(T);
    if include_decommissioning_trajectory == false
        orbit_end_index = T.Index(T.Timestamp == dts(2));
    end
    
    % Build Duration array.
    Duration = T(orbit_start_index:orbit_end_index, 3).Duration;

    % Build X, Y, Z trajactory coordinate arrays.
    X = T(orbit_start_index:orbit_end_index, 4).X;
    Y = T(orbit_start_index:orbit_end_index, 5).Y;
    Z = T(orbit_start_index:orbit_end_index, 6).Z;
    
    if is_earth_origin
        % Shift coordinates if we want X, Y, Z origin to be Earth.
        % Shift by Earth's coordinate value which corresponds to the first
        % row of the data table.
        X = X - Earth_Position(1);
        Y = Y - Earth_Position(2);
        Z = Z - Earth_Position(3);
        
        L2_Position = [-Earth_Position(1), -Earth_Position(2), -Earth_Position(3)];
        Earth_Position = [0, 0, 0];
    end
end

