%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all the previous data, plots and commands:
clear all
close all
clc

% Show or hide figures
fig_visible = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD TRAJECTORY DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_filename = 'data/Plex_Positions_With_Decommissioning_L2_Reference_Frame.csv';

% Determine whether the plot origin is that of earth or L2.
is_earth_origin = true;

% Indicate if the trajectory from earth to L2 orbit insertion should appear
% on the plot.
include_earth_to_L2_orbit_trajectory = true;

% Indicate if the decommissioning trajectory should appear on the plot.
include_decommissioning_trajectory = false;

[Earth_Position, L2_Position, insertion_distance_from_L2, Duration, X, Y, Z] = get_trajectory_positions(data_filename,...
    is_earth_origin,...
    include_earth_to_L2_orbit_trajectory,...
    include_decommissioning_trajectory);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Astronomic unit (distance Earh-Sun)(km):
AU=1.49598*10^8;

% Earth-L2 distance (km):
if is_earth_origin
    EL2=L2_Position(1);
else
    EL2=-Earth_Position(1);
end

% Sun-L2 distance (km):
SL2=AU+EL2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DISTANCES TO EARTH:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance between PLEX and Earth.
PE_distances = [];
for i = 1 : length(X)
    distance_plex_earth = norm(Earth_Position - [X(i), Y(i), Z(i)]);
    PE_distances(i) = distance_plex_earth;
end

fig_dist = figure('Name', 'Distance to Earth', 'visible', fig_visible);
set(fig_dist,'color','w');
set(fig_dist,'position',[10,10,1000,350])
plot(Duration / (3600 * 24), PE_distances);
legend('Earth')
xlabel('Time (days from launch)'); 
ylabel('Distance (km)'); 
title('Distance to Earth');
saveas(fig_dist,'img/stk/trajectory-distance-to-earth.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SUN-PLEX-EARTH (SPE) ANGLES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance between PLEX and L2.
PL2_distances = [];
for i = 1 : length(X)
    distance_plex_L2 = norm(L2_Position - [X(i), Y(i), Z(i)]);
    PL2_distances(i) = distance_plex_L2;
end

% Sun-L2-PLEX angles.
Sun_L2_Plex_Angles = [];
for i = 1 : length(X)
    angle_sun_l2_plex = acos((EL2^2 + PL2_distances(i)^2 - PE_distances(i)^2) / (2 * EL2 * PL2_distances(i)));
    Sun_L2_Plex_Angles(i) = rad2deg(angle_sun_l2_plex);
end

% Distance between PLEX and Sun.
PS_distances = [];
for i = 1 : length(X)
    distance_plex_sun = sqrt(SL2^2 + PL2_distances(i)^2 - 2 * SL2 * PL2_distances(i) * cos(deg2rad(Sun_L2_Plex_Angles(i))));
    PS_distances(i) = distance_plex_sun;
end

%{
% Sun-Earth-PLEX angles.
Sun_Earth_Plex_Angles = [];
for i = 1 : length(X)
    angles_sun_earth_plex = acos((AU^2 + PE_distances(i)^2 - PS_distances(i)^2) / (2 * AU * PE_distances(i)));
    Sun_Earth_Plex_Angles(i) = rad2deg(angles_sun_earth_plex);
end
%}

% Sun-PLEX-Earth angles.
Sun_Plex_Earth_Angles = [];
for i = 1 : length(X)
    angles_sun_plex_earth = acos((PS_distances(i)^2 + PE_distances(i)^2 - AU^2) / (2 * PS_distances(i) * PE_distances(i)));
    Sun_Plex_Earth_Angles(i) = rad2deg(angles_sun_plex_earth);
end

fig_sep_angles = figure('Name', 'Sun-Plex-Earth Angles', 'visible', fig_visible);
set(fig_sep_angles,'color','w');
set(fig_sep_angles,'position',[10,10,600,600])
plot(Duration / (3600 * 24), Sun_Plex_Earth_Angles);
legend('Earth')
xlabel('Time (days from launch)'); 
ylabel('Angle (deg)'); 
title('Sun-Plex-Earth Angles');
saveas(fig_sep_angles,'img/stk/trajectory-sun-plex-earth-angles.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT TRAJECTORY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_yx = figure('Name', 'yx-plane', 'visible', fig_visible);
set(fig_yx,'color','w');
set(fig_yx,'position',[10,10,1200,500])
% Plot the Lissajous orbit viewed on yx-plane.
%subplot(2,2,1);
plot(X, Y);
xlabel('X (km)'); 
ylabel('Y (km)'); 
title('yx-plane');
saveas(fig_yx,'img/stk/trajectory-yx-plane.png');

fig_yz = figure('Name', 'yz-plane', 'visible', fig_visible);
set(fig_yz,'color','w');
set(fig_yz,'position',[10,10,600,600])
% Plot the Lissajous orbit viewed on yz-plane.
%subplot(2,2,2);
plot(Z, Y);
xlabel('Z (km)'); 
ylabel('Y (km)');
title('yz-plane');
saveas(fig_yz,'img/stk/trajectory-yz-plane.png');

fig_zx = figure('Name', 'zx-plane', 'visible', fig_visible);
set(fig_zx,'color','w');
set(fig_zx,'position',[10,10,1200,500])
% Plot the Lissajous orbit viewed on zx-plane.
%subplot(2,2,[3,4]);
plot(X, Z);
xlabel('X (km)'); 
ylabel('Z (km)'); 
title('zx-plane');
saveas(fig_zx,'img/stk/trajectory-zx-plane.png');

%suptitle('Lissajous orbit viewed on different planes with respect to the L2 reference frame');

% Plot 3D Lissajous orbit.
fig_3d = figure('visible', fig_visible);
set(fig_3d,'color','w');
set(fig_3d,'position',[10,10,1200,500]);
plot3(X,Y,Z, 'LineWidth', 1);
view(160,30);
xlabel('X (km)'); 
ylabel('Y (km)'); 
zlabel('Z (km)');
grid on
saveas(fig_3d,'img/stk/trajectory-3d.png');

% Plot 3D Lissajous orbit with projections on 3 planes.
fig_3d_proj = figure('visible', fig_visible);
set(fig_3d_proj,'color','w');
set(fig_3d_proj,'position',[10,10,1200,500])
plot3(X,Y,Z, 'LineWidth', 1);
view(160,30);
xlabel('X (km)'); 
ylabel('Y (km)'); 
zlabel('Z (km)');
grid on
hold on
%plot(X, Y);
% Project on yx-plane.
plot3(X,Y,zeros(length(Z),1)-300000, 'r', 'LineWidth', 1, 'LineStyle', ':');

% Project on yz-plane.
plot3(zeros(length(X),1)-200000,Y,Z, 'r', 'LineWidth', 1, 'LineStyle', ':');

% Project on zx-plane.
plot3(X,zeros(length(Y),1)-1000000,Z, 'r', 'LineWidth', 1, 'LineStyle', ':');
hold off
saveas(fig_3d,'img/stk/trajectory-3d-with-projections-on-planes.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE ANIMATED TRAJECTORIES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Skip the section orbiting around Earth.
animation_start_index = 200;
%{
% Plot the Lissajous orbit viewed on yx-plane.
fig_anim_yx = figure('visible', fig_visible);
set(fig_anim_yx,'position',[10,10,1200,500])
set(fig_anim_yx,'color','w');
filename = 'img/stk/trajectory-animated-yx-plane.gif';
animate_2d_plot(fig_anim_yx, X, Y, 'yx-plane', 'X (km)', 'Y (km)', animation_start_index, [-200000 1800000], [-1000000 1000000], filename)

% Plot the Lissajous orbit viewed on yz-plane.
fig_anim_yz = figure('visible', fig_visible);
set(fig_anim_yz,'position',[10,10,600,600])
set(fig_anim_yz,'color','w');
filename = 'img/stk/trajectory-animated-yz-plane.gif';
animate_2d_plot(fig_anim_yz, Z, Y, 'yz-plane', 'Z (km)', 'Y (km)', animation_start_index, [-300000 200000], [-1000000 1000000], filename)

% Plot the Lissajous orbit viewed on zx-plane.
fig_anim_zx = figure('visible', fig_visible);
set(fig_anim_zx,'position',[10,10,1200,500])
set(fig_anim_zx,'color','w');
filename = 'img/stk/trajectory-animated-zx-plane.gif';
animate_2d_plot(fig_anim_zx, X, Z, 'zx-plane', 'X (km)', 'Z (km)', animation_start_index, [-200000 1800000], [-300000 200000], filename)
%}

%{
% Plot 3D Lissajous orbit with projections on 3 planes.
fig_anim_3d = figure('visible', fig_visible);
set(fig_anim_3d,'position',[10,10,1200,500]);
set(fig_anim_3d,'color','w');
filename = 'img/stk/trajectory-animated-3d.gif';
animate_3d_plot(fig_anim_3d, X, Y, Z, '', 'X (km)', 'Y (km)', 'Z (km)', animation_start_index, [-200000 1800000], [-1000000 1000000], [-300000 200000], 160, 30, false, filename);


% Plot 3D Lissajous orbit with projections on 3 planes.
fig_anim_3d_with_proj = figure('visible', fig_visible);
set(fig_anim_3d_with_proj,'position',[10,10,1200,500]);
set(fig_anim_3d_with_proj,'color','w');
filename = 'img/stk/trajectory-animated-3d-with-projections-on-planes.gif';
animate_3d_plot(fig_anim_3d_with_proj, X, Y, Z, '', 'X (km)', 'Y (km)', 'Z (km)', animation_start_index, [-200000 1800000], [-1000000 1000000], [-300000 200000], 160, 30, true, filename);
%}

