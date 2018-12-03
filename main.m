clearvars
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD TRAJECTORY DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_filename = 'Plex_Positions_With_Decommissioning_L2_Reference_Frame.csv';

% Determine whether the plot origin is that of earth or L2.
is_earth_origin = true;

% Indicate if the trajectory from earth to L2 orbit insertion should appear
% on the plot.
include_earth_to_L2_orbit_trajectory = true;

% Indicate if the decommissioning trajectory should appear on the plot.
include_decommissioning_trajectory = true;

[X, Y, Z] = get_trajectory_positions(data_filename,...
    is_earth_origin,...
    include_earth_to_L2_orbit_trajectory,...
    include_decommissioning_trajectory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT TRAJECTORY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
% Plot the Lissajous orbit viewed on yx-plane.
subplot(2,2,1);
plot(X, Y);
xlabel('X (km)') 
ylabel('Y (km)') 
title('yx-plane');

% Plot the Lissajous orbit viewed on yz-plane.
subplot(2,2,2);
plot(Z, Y);
xlabel('Z (km)') 
ylabel('Y (km)') 
title('yz-plane');

% Plot the Lissajous orbit viewed on zx-plane.
subplot(2,2,[3,4]);
plot(X, Z);
xlabel('X (km)') 
ylabel('Z (km)') 
title('zx-plane');

suptitle('Lissajous orbit viewed on different planes with respect to the L2 reference frame');

f2 = figure;
plot3(X,Y,Z, 'r', 'LineWidth', 1);
xlabel('X (km)') 
ylabel('Y (km)') 
zlabel('Z (km)')