function [A_h,E_h] = elevation_azimuthal(d_i,A_helio)
clc

%parameters
lttd = 26.28; %latitude
longitude=73.03;
slope = 0; %terrain slope
twrht = 10; %tower height
ht=twrht;
lr=0; %length of receiver
h_i=1; %height of heliostat

%incdnt %angle of incidennce of sun's rays
%d_i=20;
%A_helio=(90+90); %azimuth angle of heliostat
day=30; %day of year
td=5.5; %time zone in hours
time = 16; %local time in hours

lstm=15*td;
b=360/365*(day-81);
eot=9.87*sind(2*b)-7.53*cosd(b)-1.5*sind(b);
tc=4*(longitude-lstm)+eot;
lst = time + (tc / 60);
hr = 15*(lst - 12);

lambda = atand(d_i/(ht-(0.5*lr)-(d_i*tand(slope))-h_i)); %angle betweeen rflected ray and vertical of receiver

%declination angle
decl = 23.45*sind((360/365)*(284+day));
%decl= -23.45*cos((360/365)*(10+dy))*pi/180;
%decl = -asin(0.39779*cos(0.98565*(dy+10)+1.914*sin(0.98565*(dy-2))));

%elevation of sun
alfa = asind((cosd(lttd)*cosd(decl)*cosd(hr))+(sind(lttd)*sind(decl)));

%azimuth angle of sun
A_sun = acosd((sind(alfa)*sind(lttd)-sind(decl))/(cosd(alfa)*cosd(lttd)));

%azimuth of normal of heliostat's central pt
A_h = atand((sind(A_helio)*sind(lambda)-sind(A_sun)*cosd(alfa))/(cosd(A_helio)*sind(lambda)-cosd(A_sun)*cosd(alfa)));

%angle of incidence
incdnt_n = acosd(0.5*(2^0.5)*((sind(alfa)*cosd(lambda))-(cosd(A_helio-A_sun)*cosd(alfa)*sind(lambda))+1)^0.5);

%elevation of normal of heliostat's central pt
E_h = acosd((sind(alfa)+cosd(lambda))/(2*cosd(incdnt_n)));
end