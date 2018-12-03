%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all the previous data, plots and commands:
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD TRAJECTORY DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_filename = 'Plex_Positions_With_Decommissioning_L2_Reference_Frame.csv';

% Determine whether the plot origin is that of earth or L2.
is_earth_origin = false;

% Indicate if the trajectory from earth to L2 orbit insertion should appear
% on the plot.
include_earth_to_L2_orbit_trajectory = false;

% Indicate if the decommissioning trajectory should appear on the plot.
include_decommissioning_trajectory = false;

[X, Y, Z] = get_trajectory_positions(data_filename,...
    is_earth_origin,...
    include_earth_to_L2_orbit_trajectory,...
    include_decommissioning_trajectory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT DATA:
%%%%%%%%%%%%%
% time length (s):
length = 6*364*24*3600;
% sample time (s):
s = 100;

% Inital coordinates y,z in the Lissajous orbit (km)
% (only necessary for calculating Lissajous phases):
Y0=-13000*sin(45*(pi/180));
Z0=13000*cos(45*(pi/180));

% DEFINING LISSAJOUS ORBIT PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lissajous orbit periods (days):
Txy=177.655*24*3600;
Tz=184.0*24*3600;

% Lissajous orbit amplitudes (km):
Az = max(Z) - min(Z);
Ay = max(Y) - min(Y);
Ax = max(X) - min(X);

%Az=4.5*10^5;
%Ay=1.5*10^6;
%Ax=5*10^5; 
%Ax=Ay/3.1872293; %Ax is related with Ay by the constant C2 = 3.1872293

%Maximum Lissajous orbit total time (with or without solar eclipses)(years)
years=8; %for example: 8 years

% Astronomic unit (distance Earh-Sun)(km):
AU=1.49598*10^8;

% Earth-L2 distance (km):
XL=1507683;

% Earth shadow radium (km):
r=13000;

% Creating Earth shadow in the yz-plane:
for tes=[0:1:10000];
    ry(tes+1)=r*-sin(tes);
    rz(tes+1)=r*cos(tes);
end

% Calculating Lissajous orbit rates:
Wxy=2*pi/Txy;
Wz=2*pi/Tz;

% Calculating Lissajous orbit phases:
Pxy=asin(Y0/-Ay);
Pz=acos(Z0/Az);

% DEFINING SCANNING LAW PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining motion rates:
w = ((2*pi)/8766)*(1./3600); % translation
wp = ((0.17*2*pi)/360)*(1./3600); % precession
ws = ((120*2*pi)/360)*(1./3600); % spinning

% Initial value for precession vector:
p1n = sqrt(0.5); p2n = -sqrt(0.5); p3n = 0;

% Initial value for spinning vector:
p1(1) = 0.0; p2(1) = 1.0; p3(1) = 0.0;

%initializing errors:
sigma_wp = ((0.1*2*pi)/(3*360))*(1./3600); % precession
sigma_ws = ((1.2*2*pi)/(3*360))*(1./3600); % spinning
sigma_point = ((5*2*pi)/(3*60*360)); % pointing

% INITIALIZING VECTORS:
%%%%%%%%%%%%%%%%%%%%%%%
% (This reserves enough memory space and optimizes the simulation process)
vl=(length/s)+1;
x=zeros(1,vl); y=zeros(1,vl); z=zeros(1,vl);
x1=zeros(1,vl); y1=zeros(1,vl); z1=zeros(1,vl);
x1L=zeros(1,vl); y1L=zeros(1,vl); z1L=zeros(1,vl);
x2=zeros(1,vl); y2=zeros(1,vl); z2=zeros(1,vl);
x2i=zeros(1,vl); y2i=zeros(1,vl); z2i=zeros(1,vl);
x3=zeros(1,vl); y3=zeros(1,vl); z3=zeros(1,vl);
x31=zeros(1,vl); y31=zeros(1,vl); z31=zeros(1,vl);
x3i=zeros(1,vl); y3i=zeros(1,vl); z3i=zeros(1,vl);
xp2=zeros(1,vl); yp2=zeros(1,vl); zp2=zeros(1,vl);
xp2i=zeros(1,vl); yp2i=zeros(1,vl); zp2i=zeros(1,vl);
nlx=zeros(1,vl); nly=zeros(1,vl); nlz=zeros(1,vl);
nx=zeros(1,vl); ny=zeros(1,vl); nz=zeros(1,vl);
nx1_nol=zeros(1,vl); ny1_nol=zeros(1,vl); nz1_nol=zeros(1,vl);
nx1=zeros(1,vl); ny1=zeros(1,vl); nz1=zeros(1,vl);
nx3=zeros(1,vl); ny3=zeros(1,vl); nz3=zeros(1,vl);
nx31=zeros(1,vl); ny31=zeros(1,vl); nz31=zeros(1,vl);
nx31L=zeros(1,vl); ny31L=zeros(1,vl); nz31L=zeros(1,vl);
nx3L=zeros(1,vl); ny3L=zeros(1,vl); nz3L=zeros(1,vl);
nx3i=zeros(1,vl); ny3i=zeros(1,vl); nz3i=zeros(1,vl);
nxL=zeros(1,vl); nyL=zeros(1,vl); nzL=zeros(1,vl);
nxi=zeros(1,vl); nyi=zeros(1,vl); nzi=zeros(1,vl);
lx=zeros(1,vl); ly=zeros(1,vl); lz=zeros(1,vl);
nxi_nol=zeros(1,vl); nyi_nol=zeros(1,vl); nzi_nol=zeros(1,vl);
nx3i_nol=zeros(1,vl); ny3i_nol=zeros(1,vl); nz3i_nol=zeros(1,vl);
anglenr1=zeros(1,vl); anlgenr2=zeros(1,vl); anglenr3=zeros(1,vl);
anglenr2noise=zeros(1,vl); anlgenr31noise=zeros(1,vl); anglenr3noise=zeros(1,vl);
new_theta=zeros(1,vl); new_fi=zeros(1,vl);
wpe=zeros(1,vl); wse=zeros(1,vl);
x=zeros(1,vl);
k=zeros(1,vl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LISSAJOUS ORBIT SIMULATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont=1;
ld=r+1;

for t = 0:s:length
    if (ld >= r) %(if the orbit didn't cross the Earth-shadow:)
        % Calculating position with the Lissajous orbit equations:
        lx(cont)=Ax*cos(Wxy*t+Pxy);
        ly(cont)=-Ay*sin(Wxy*t+Pxy);
        lz(cont)=Az*cos(Wz*t+Pz);
        
        % Normalizing vector:
        nlx(cont) = lx(cont)/sqrt(lx(cont)^2+ly(cont)^2+lz(cont)^2);
        nly(cont) = ly(cont)/sqrt(lx(cont)^2+ly(cont)^2+lz(cont)^2);
        nlz(cont) = lz(cont)/sqrt(lx(cont)^2+ly(cont)^2+lz(cont)^2);
        
        % calculating distance to Earth shadow
        % (only necessary for stop the Lissajous orbit when enters to Earth shadow)
        ld=sqrt(((ly(cont))^2)+((lz(cont))^2))+1;
        
        dlo=length; % duration of the Lissajous orbit
        cont=cont+1;
        
    else %(if the orbit crossed the Earh-shadow:)
        dlo=t; % duration of the Lissajous orbit
        break %stops the Lissajous orbit simulation
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCANNING LAW WITHOUT NOISE (WITH AND WITHOUT LISSAJOUS):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cont = 1;
for t = 0:s:length
    
    % TRANSLATION (WITHOUT LISSAJOUS):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1(cont) = (AU+XL)*cos(w*t);
    y1(cont) = (AU+XL)*sin(w*t);
    z1(cont) = 0*t;
    
    % Normalizing vector:
    nx1_nol(cont) = x1(cont)/sqrt(x1(cont)^2+y1(cont)^2+z1(cont)^2);
    ny1_nol(cont) = y1(cont)/sqrt(x1(cont)^2+y1(cont)^2+z1(cont)^2);
    nz1_nol(cont) = 0;
    
    % TRANSLATION (WITH LISSAJOUS):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1L(cont)=lx(cont)+x1(cont);
    y1L(cont)=ly(cont)+y1(cont);
    z1L(cont)=lz(cont)+z1(cont);
    
    % Normalizing vector:
    nx1(cont) = x1L(cont)/sqrt(x1L(cont)^2+y1L(cont)^2+z1L(cont)^2);
    ny1(cont) = y1L(cont)/sqrt(x1L(cont)^2+y1L(cont)^2+z1L(cont)^2);
    nz1(cont) = z1L(cont)/sqrt(x1L(cont)^2+y1L(cont)^2+z1L(cont)^2);
    % PRECESSION (WITHOUT LISSAJOUS):
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cont==1
        xp2i(cont) = 0;
        yp2i(cont) = 1;
        zp2i(cont) = 0;
    else
        xp2i(cont) = 0;
        yp2i(cont) = yp2i(cont-1).*cos(wp*s)-zp2i(cont-1).*sin(wp*s);
        zp2i(cont) = yp2i(cont-1).*sin(wp*s)+zp2i(cont-1).*cos(wp*s);
    end
    
    x2i(cont) = nx1_nol(cont)+xp2i(cont).*cos(w*t)-yp2i(cont).*sin(w*t);
    y2i(cont) = ny1_nol(cont)+xp2i(cont).*sin(w*t)+yp2i(cont).*cos(w*t);
    z2i(cont) = nz1_nol(cont)+zp2i(cont);

    % Normalizing vector:
    nxi_nol(cont) = x2i(cont)./sqrt(x2i(cont).^2+y2i(cont).^2+z2i(cont).^2);
    nyi_nol(cont) = y2i(cont)./sqrt(x2i(cont).^2+y2i(cont).^2+z2i(cont).^2);
    nzi_nol(cont) = z2i(cont)./sqrt(x2i(cont).^2+y2i(cont).^2+z2i(cont).^2);

    % PRECESION (WITH LISSAJOUS):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cont==1
        xp2i(cont) = 0;
        yp2i(cont) = 1;
        zp2i(cont) = 0;
    else
        xp2i(cont) = 0;
        yp2i(cont) = yp2i(cont-1).*cos(wp*s)-zp2i(cont-1).*sin(wp*s);
        zp2i(cont) = yp2i(cont-1).*sin(wp*s)+zp2i(cont-1).*cos(wp*s);
    end

    x2i(cont) = nx1(cont)+xp2i(cont).*cos(w*t)-yp2i(cont).*sin(w*t);
    y2i(cont) = ny1(cont)+xp2i(cont).*sin(w*t)+yp2i(cont).*cos(w*t);
    z2i(cont) = nz1(cont)+zp2i(cont);

    % Normalizing vector:
    nxi(cont) = x2i(cont)./sqrt(x2i(cont).^2+y2i(cont).^2+z2i(cont).^2);
    nyi(cont) = y2i(cont)./sqrt(x2i(cont).^2+y2i(cont).^2+z2i(cont).^2);
    nzi(cont) = z2i(cont)./sqrt(x2i(cont).^2+y2i(cont).^2+z2i(cont).^2);

    cont = cont+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LISSAJOUS ORBIT:
%%%%%%%%%%%%%%%%%%
figure('Name','Lissajous YZ-plane')
plot(ly,lz,ry,rz,'-')
%axis([-100000 100000 -100000 100000])
xlabel('Y (km)')
ylabel('Z (km)')
title('YZ-plane PLEX orbit')
grid on

figure('Name','Lissajous YX-plane')
plot(ly,lx,'-')
%axis([-100000 100000 -100000 100000])
xlabel('Y (km)')
ylabel('X (km)')
title('YX-plane PLEX orbit')
grid on

figure('Name','Lissajous XZ-plane')
plot(lx,lz,'-')
%axis([-100000 100000 -100000 100000])
xlabel('X (km)')
ylabel('Z (km)')
title('XZ-plane PLEX orbit')
grid on

figure('Name','Lissajous 3D')
plot3(lx,ly,lz,'-')
%axis([-100000 100000 -100000 100000 -100000 100000])
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('3D PLEX orbit')
grid on

% LISSAJOUS ORBIT TIME:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indicating the duration of the Lissajous orbit (years):
disp('* Duration of the Lissajous orbit (years):')
disp(dlo/(365.25*24*3600))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%