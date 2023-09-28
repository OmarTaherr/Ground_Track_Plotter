%% Plot the ground track on the planisphere and 3D Orbit.
%% BY:  Omar Essam AbuElTaher && Mahmmed Ibrahim AbdElHay
%% ------------------------------------------------------
clc
clear
close all
%% Constant
Step= 1;            % Step of time                              [sec]
RE = 6371;          % Earth's radius                            [km]
muE = 398600.4418;  % Earth gravitational parameter             [km^3/sec^2]
wE = (2*pi/86164);  % Earth rotation velocity aorund z-axis     [rad/sec]      86164 [sec] = 23:56:4 [hour:minute:second]
%% Choose type of input Six classical orbital element or R & V Vectors
disp(' ');
disp('           Ground track on the planisphere and 3D              ');
disp(' -------------------------------------------------------------- ');
disp(' ******************* Choose type of input ********************* ');
disp(' ------- Six classical orbital element or R & V vectors ------- ');
disp(' ');
disp(' 1) Six classical orbital element');
disp(' 2) R & V vectors');
grpp=0;
while (grpp~=1) && (grpp~=2)
    grpp = input('  Digit [1],[2] >> ');
end
%% Six classical orbital element Initial Value
if grpp==1
    disp(' Please Enter Six clasical orbital element Initial Value ');
    RAAN0         = input(' Right Ascension of Ascendent Node    [0,360]    RAAN0   (deg) = ');
    w0            = input(' Argument of perigee                  [0,360]    w0      (deg) = ');
    v0            = input(' True anomaly at the departure        [0,360]    v0      (deg) = ');
    i0            = input(' Inclination                          [0,180]    i0      (deg) = ');
    a            = input(' Semi-Major axis                      (>6371)     a      (km)  = ');
    ecc_max       = sprintf('%6.4f',1-RE/a);                  % maximum value of eccentricity allowed
    e           = input([' Eccentricity                       (<',ecc_max,')     e            = ']);
    % ------------------------------- CONV ---------------------------------%
    RAAN0  = RAAN0*pi/180;        % RAAN                          [rad]
    w0     = w0*pi/180;           % Argument of perigee           [rad]
    v0     = v0*pi/180;           % True anomaly at the departure [rad]
    i0     = i0*pi/180;           % inclination                   [rad]
    % ----------------------------------------------------------------------%
end
%% R & V vectors
if grpp==2
    R_vec=zeros(1,3);
    V_vec=zeros(1,3);
    disp(' ');
    disp(' Please Enter R vector cordinate ');
    R_vec(1,1)    = input(' X-cordinate of R vector   [km]  Rx = ');
    R_vec(1,2)    = input(' Y-cordinate of R vector   [km]  Ry = ');
    R_vec(1,3)    = input(' Z-cordinate of R vector   [km]  Rz = ');
    disp(' ');
    disp(' Please Enter V vector cordinate ');
    V_vec(1,1)    = input(' X-cordinate of V vector   [km/s]  Vx = ');
    V_vec(1,2)    = input(' Y-cordinate of V vector   [km/s]  Vy = ');
    V_vec(1,3)    = input(' Z-cordinate of V vector   [km/s]  Vz = ');
    R_mag=norm(R_vec);
    V_mag=norm(V_vec);
    % ------------------------------ VECTORS ---------------------------------%
    I_vec=[1,0,0];
    I_mag=norm(I_vec);
    k_vec=[0,0,1];
    k_mag=norm(k_vec);
    h_vec=cross(R_vec,V_vec);
    h_mag=norm(h_vec);
    n_vec=cross(k_vec,h_vec);
    n_mag=norm(n_vec);
    % ------------------------------------------------------------------------%
    Spacific_energy=((V_mag^2)/2)-(muE/R_mag);                                 % Spacific energy
    a=-muE/(2*Spacific_energy);                                                % Semi-Major axis                      [km]
    e_vec=(1/muE)*(((V_mag^2)-(muE/R_mag))*R_vec-(dot(R_vec,V_vec))*V_vec);
    e=norm(e_vec);                                                             % Eccentricity
    i0=acos(dot(k_vec,h_vec)/(k_mag*h_mag));                                   % Inclination                          [rad]
    RAAN0=acos(dot(I_vec,n_vec)/(I_mag*n_mag));                                % Right Ascension of Ascendent Node    [rad]
    w0=acos(dot(n_vec,e_vec)/(n_mag*e));                                       % Argument of perigee                  [rad]
    v0=acos(dot(e_vec,R_vec)/(e*R_mag));                                       % True anomaly at the departure        [rad]
    % PRINT Initial Values calculated from R & V Vectors
    fprintf('\n Right Ascension of Ascendent Node     RAAN0   [%5.0f deg]',RAAN0*180/pi);
    fprintf('\n Argument of perigee                   w0      [%5.0f deg]',w0*180/pi);
    fprintf('\n True anomaly at the departure         v0      [%5.0f deg]',v0*180/pi);
    fprintf('\n Inclination                           i0      [%5.0f deg]',i0*180/pi);
    fprintf('\n Semi-Major axis                       a       [%10.1f km]',a);
    fprintf('\n Eccentricity                          e       [%10.8f ]',e);
end
%% Animatied or not ??
disp(' ');
num_of_orbits =input(' Number of Orbits                           num_of_orbits      = ');
disp(' ------------------------------------------------------- ');
disp(' ************* Choose Planisphere Options ************** ');
disp(' -------- Animation of the Track to be printed --------- ');
disp(' 1) Not Animatied');
disp(' 2) Animatied');
grp=0;
while (grp~=1) && (grp~=2)
    grp = input('  Digit [1],[2] >> ');
end
if grp==2
    Step          =input(' Step                                 (sec)                    = ');
end
%% ORBIT COMPUTATION
rp = a*(1-e);                % radius of perigee             [km]
ra = a*(1+e);                % radius of apogee              [km]
Vp = sqrt(muE*(2/rp-1/a));   % velocity at the perigee       [km/s]
Va = sqrt(muE*(2/ra-1/a));   % velocity at the  apogee       [km/s]
n  = sqrt(muE/a^3);          % mean motion                   [rad/s]
P  = 2*pi/n;                 % period                        [s]
%% PRINT SOME DATAS
hours   = floor(P/3600);                   % hours   of the orbital period
minutes = floor((P-hours*3600)/60);        % minutes of the orbital period
seconds = floor(P-hours*3600-minutes*60);  % seconds of the orbital period
disp(' ');
fprintf('\n Radius of perigee [%10.3f km]       Altitude of perigee [%10.3f km]',rp,rp-RE);
fprintf('\n Radius of  apogee [%10.3f km]       Altitude of  apogee [%10.3f km]',ra,ra-RE);
fprintf('\n Velocity at the perigee [%6.4f km/s]   Velocity at the apogee [%6.4f km/s]',Vp,Va);
fprintf('\n Orbital Period    [%3d h: %3d m: %3d s] ',hours,minutes,seconds);
fprintf('   = [%10.2f s]\n',P);
%% Calculation of E0 & M0
E0=2*atan(tan(v0/2)*sqrt((1-e)/(1+e)));
if E0<0
    E0=E0+2*pi;
end
M0=E0-e*sin(E0);
% Time
t0=(P/(2*pi))*M0;
tf=num_of_orbits*P;
t=t0:Step:tf+t0;
M=(2*pi/P).*t;
% Loop to calculate E with Iteration
E=size(t);
for Loop=1:length(t)
    E(Loop)=anom_ecc(M(Loop),e);
end
% True anomaly
v= 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));          % true anomaly              [rad]
% preifocal plan
r=a*(1-e^2)./(1+e*cos(v));                       %radius to orbit with time
Co_0(1,:)=r.*cos(v);
Co_0(2,:)=r.*sin(v);
Co_0(3,:)=0;
% J2 Effect
RAAN_J2=-41708925*a^(-7/2)*cos(i0)*(1-e^2)^-2;                  % [rad/sec]
w_J2=20854462.5*a^(-7/2)*(4-5*sin(i0)^2)*(1-e^2)^-2;            % [rad/sec]
RAAN=zeros(1,length(t));
w=zeros(1,length(t));
RAAN=RAAN+RAAN0;
w=w+w0;
for Loop=1:length(t)
    RAAN(Loop)=RAAN(Loop)+RAAN_J2*t(Loop) ;
    w(Loop)=w(Loop)+w_J2*t(Loop) ;
end
% Transformation from preifocal
Co_1=size(Co_0);
Co_2=size(Co_0);
Co_3=size(Co_0);
for Loop=1:length(t)
    Co_1(1:3,Loop)=[cos(w(Loop)) -sin(w(Loop)) 0;sin(w(Loop)) cos(w(Loop)) 0;0 0 1]*Co_0(:,Loop);
    Co_2(1:3,Loop)=[1 0 0;0 cos(i0) -sin(i0);0 sin(i0) cos(i0)]*Co_1(:,Loop);
    Co_3(1:3,Loop)=[cos(RAAN(Loop)) -sin(RAAN(Loop)) 0;sin(RAAN(Loop)) cos(RAAN(Loop)) 0;0 0 1]*Co_2(:,Loop);
end
Co_4=Co_3;
theta=wE*(t-t0);
for Loop_2 =1:length(t)
    Co_4(:,Loop_2)=[cos(-theta(Loop_2)) -sin(-theta(Loop_2)) 0;sin(-theta(Loop_2)) cos(-theta(Loop_2)) 0;0 0 1]*Co_3(:,Loop_2);    % real x-y-z coordinates for orbit
end
%% Plot the orbit & earth lines of the surface in prefocal plane
figure('Name','Prefocal Plane')
plot(Co_0(1,:),Co_0(2,:),"r .",'LineWidth',1)         % Orbit
xlabel("X [km]",'interpreter','latex')
ylabel("Y [km]",'interpreter','latex')
axis equal
hold on
grid on
u=0:1:360;                                            % Earth's surface
plot(RE*cosd(u),RE*sind(u),'Color','#77AC30','LineWidth',2)
text(-4900,0,texlabel('Earth'),'Color',"#77AC30",'FontSize',25,'interpreter','latex');
%% Plot 3D orbit & Track
figure('Name','3D ORBIT & TRACK')
earth_sphere('km')
hold on
grid on
[x_equat,y_equat]=meshgrid(linspace(-1.8*a,1.8*a,25),linspace(-1.8*a,1.8*a,25));          % Equatorial Plane
surf(x_equat,y_equat,zeros(size(x_equat)),'FaceColor','none','EdgeColor','[0.4 0.4 0.4]')
text(-2*a,1.2*a,0,'Equatorial Plane','Color','[0.5 0.5 0.5]')
plot3([0,2*a],[0,0],[0,0],'--','Color',"#A2142F",'LineWidth',2);                          % X direction
plot3([0,0],[0,2*a],[0,0],'--','Color',"#77AC30",'LineWidth',2);                          % Y direction
plot3([0,0],[0,0],[0,2*a],'--','Color',"#0072BD",'LineWidth',2);                          % Z direction
text(2*a+1000,0,0,texlabel('X-axis'),'Color',"#A2142F",'FontSize',15,'interpreter','latex');
text(0,2*a+1000,0,texlabel('Y-axis'),'Color','#77AC30','FontSize',15,'interpreter','latex');
text(0,0,2*a+1000,texlabel('Z-axis'),'Color',"#0072BD",'FontSize',15,'interpreter','latex');
plot3(Co_4(1,:)./r*(RE+50),Co_4(2,:)./r*(RE+50),Co_4(3,:)./r*(RE+50),'Color',"r",'LineWidth',3)         % Track
plot3(Co_3(1,:),Co_3(2,:),Co_3(3,:),'Color',"[0.95 0.5 0.05]",'LineWidth',3)                            % Orbit
%% Plot track in 2D (Planisphere)
Latitude=size(v);
Longitude=size(v);
for Loop=1:length(v)
    R=(Co_4(1,Loop).^2+Co_4(2,Loop).^2+Co_4(3,Loop).^2).^0.5;
    Latitude(Loop)=asin(Co_4(3,Loop)./R);
    if Co_4(2,Loop)>=0
        Longitude(Loop)=acos(Co_4(1,Loop)./R./cos(Latitude(Loop)));
    else
        Longitude(Loop)=-acos(Co_4(1,Loop)./R./cos(Latitude(Loop)));
    end
end
% --------------------------------- CONV ---------------------------------%
Latitude  = Latitude*180/pi;        % RAAN                          [deg]
Longitude = Longitude*180/pi;       % Argument of perigee           [deg]
% ------------------------------------------------------------------------%
figure('Name','Planisphere Track');
hold on;
set(gca,'XTick',-180:90:180,'XTickMode','manual');
set(gca,'YTick',-90:45:90,'YTickMode','manual');
image_file = 'land_par_mer.jpg';
cdata      = imread(image_file);
imagesc([-180,180],[90,-90],cdata);
xlabel('Longtitude','interpreter','latex')
ylabel('Latitude','interpreter','latex')
hold on;
grid on
% Start Point
plot(Longitude(1),Latitude(1,1),'Color',"g",'Marker','.','MarkerSize',25)
text(Longitude(1)-15,Latitude(1,1)+5,texlabel('Start'),'Color',"g",'FontSize',20,'interpreter','latex');
% Track
if grp==1
    plot(Longitude,Latitude,'r .','MarkerSize',4);
elseif grp==2
    for Loop_3=1:length(Longitude)
        hold on
        plot(Longitude(Loop_3),Latitude(1,Loop_3),'r .','MarkerSize',4);
        pause(1e-1000000)
    end
end
% End Point
plot(Longitude(length(Longitude)),Latitude(1,length(Latitude)),'Color',"g",'Marker','.','MarkerSize',25)
text(Longitude(length(Longitude))-10,Latitude(1,length(Latitude))-5,texlabel('End'),'Color',"g",'FontSize',20,'interpreter','latex');
%% Function to calculate the mean Eccentricity
function [E] = anom_ecc(M,e)
E = M;
ratio = 1;
err = 1e-10;
while (ratio>err)
    o = M+e*sin(E);
    ratio = abs(abs(E)-abs(o));
    E = o;
end
end