clear;
clc;

%parameters
lttd = 26.3*pi/180; %latitude
slope = 0; %terrain slope
twrht = 100; %tower height
ht=twrht;
lnth=1;
wdth=1;
lr=0; %length of receiver
i_l=lnth; %length of mirror
i_w=wdth; %width of mirror
h_i=1; %height of heliostat
d_i=40; %distance of heliostat from tower
d_s=1; %safety distance between heliostats

%incdnt %angle of incidennce of sun's rays
A_helio=0; %azimuth angle of heliostat
d_i=ht;
dy=100; %day of year
time=1345; %time of day
hr=(time-1200)*pi/12; %hr angle
lambda = atan(d_i/(ht-(0.5*lr)-(d_i*tan(slope))-h_i)); %angle betweeen rflected ray and vertical of receiver

%declination angle
decl = 23.45*sin((360/365)*(284+dy));

%elevation of sun
alfa = asin(cos(lttd)*cos(decl)*cos(hr) + sin(lttd)*sin(decl));

%azimuth angle of sun
A_sun = acos((sin(alfa)*sin(lttd)-sin(decl))/cos(alfa)*cos(lttd));

%azimuth of normal of heliostat's central pt
A_h = atan((sin(A_helio)*sin(lambda)-sin(A_sun)*cos(alfa))/(cos(A_helio)*sin(lambda)-cos(A_sun)*cos(alfa)));

%{
%no blocking condition
%while(d_i<500)

    a = ((0.5*i_l)^2)-((ht-(h_i+d_i*tan(slope)))^2);
    b = 2*d_i*ht*(ht-(h_i+d_i*tan(slope)));
    c = (ht^2)*(((0.5*i_l)^2)-(d_i^2));
    d_i_0 = (-b-sqrt((b^2)-4*a*c))/(a*2);

    A = (ht+d_i_0*tan(slope))^2;
    B = -2*d_i_0*(ht-h_i)*(ht+d_i_0*tan(slope));
    C = (d_i_0^2)*(((ht-h_i)^2)-((0.5*i_l)^2)*((ht^2)+(d_i_0^2)));
    d_i_1 = (-B+sqrt((B^2)-4*A*C))/(A*2)
    d_i = d_i_1;

%end
%adjacent helios
DM = i_l*((sqrt(1+i_w/i_l))+d_s)
%}
%cos of angle of incidence
incdnt_n = (0.2*2^0.5)*(sin(alfa)*cos(lambda)-cos(A_helio-A_sun)*cos(alfa)*sin(lambda)+1)^0.5;

%elevation of normal of heliostat's central pt
E_h = acos((sin(alfa)*cos(lambda))/2*incdnt_n);
%{
%shading and blocking factor
%mirrors are divided into mxn grid identifying each as kxl
X_mkl = i_l*(l-0.5*n-0.5)/n;
Y_mkl = i_w*(0.5+0.5*m-k)/l;
%shading conddition
if((0.5-n*Ymi/i_w < k && k <n+0.5) && (0.5 < l && l < n+0.5+n*Xmi/i_l))
Ymi<0 && Xmi<0;
if((0.5 < k && k <n+0.5-n*Ymi/i_w) && (0.5 < l && l < n+0.5+n*Xmi/i_l))
Ymi>0 && Xmi<0;
if((0.5-n*Ymi/i_w < k && k <n+0.5) && (0.5+n*Xmi/i_l < l && l < n+0.5))
Ymi<0 && Xmi>0;
if((0.5< k && k <n+0.5-n*Ymi/i_w) && (0.5+n*Xmi/i_l < l && l < n+0.5))
Ymi>0 && Xmi>0;
end
%blocking condition
if((0.5-n*Ymr/i_w < k && k <n+0.5) && (0.5 < l && l < n+0.5+n*Xmr/i_l))
Ymr<0 && Xmr<0;
if((0.5 < k && k <n+0.5-n*Ymr/i_w) && (0.5 < l && l < n+0.5+n*Xmr/i_l))
Ymr>0 && Xmr<0;
if((0.5-n*Ymr/i_w < k && k <n+0.5) && (0.5+n*Xmr/i_l < l && l < n+0.5))
Ymr<0 && Xmr>0;
if((0.5< k && k <n+0.5-n*Ymr/i_w) && (0.5+n*Xmr/i_l < l && l < n+0.5))
Ymr>0 && Xmr>0;
end
%shading coordinates on central
Xmi= (Yo*cos(A)-X0*sin(A))*cos(A_h-A)+(Zo*cos(alfa)-sin(alfa)*(Xo*cos(A) - Yo*sin(A)))*sin(alfa)*sin(A_h-A)-(Zi - (sin(E_h)*(Xo*cos(A_h) - Yo*sin(A_h))+Zo*cos(E_h))/incdnt_n)*cos(alfa)*sin(A_h-A);
Ymi= -(Yo*cos(A)-X0*sin(A))*cos(E_h)*sin(A_h-A)+(Zo*cos(alfa)-sin(alfa)*(Xo*cos(A) - Yo*sin(A)))*(cos(E_h)*cos(A_h-A)*sin(theta_s)+sin(E_h)*cos(alfa))-(Zi - (sin(E_h)*(Xo*cos(A_h) - Yo*sin(A_h))+Zo*cos(E_h))/incdnt_n)*(cos(E_h)*cos(A_h-A)*cos(alfa)-sin(E_h)*sin(alfa));
%blocking coordinates on central
Xmr= (Yo*cos(A_helio)-X0*sin(A_helio))*cos(A_helio_h-A_h)+(cos(lambda)*(Xo*cos(A_helio) - Yo*sin(A_helio))+Zo*sin(lambda))*cos(lambda)*sin(A_helio-A_h)-(Zr- (sin(E_h)*(Xo*cos(A_h) - Yo*sin(A_h))+Zo*cos(E_h))/incdnt_n)*sin(lambda)*sin(A_helio-A_h);
Ymr= (Yo*cos(A_helio)-X0*sin(A_helio))*cos(E_h)*sin(A_helio-A_h)-(cos(lambda)*(Xo*cos(A_helio) - Yo*sin(A_helio))+Zo*sin(lambda))*(cos(E_h)*cos(lambda)*cos(A_helio-A_h)-sin(E_h)*sin(lambda))+(Zr- (sin(E_h)*(Xo*cos(A_h) - Yo*sin(A_h))+Zo*cos(E_h))/incdnt_n)*(cos(E_h)*cos(A_helio-A_h)*sin(lambda)+sin(E_h)*cos(lambda));

%interception factor
%incidence on mirror
del_tx= rho_sun*cos(theta_sun);
del_ty= rho_sun*sin(theta_sun);
del_mix= del_tx*cos(A_h-A)+del_ty*sin(alfa)*sin(A_h-A);
del_miy= -del_tx*cos(E_h)*sin(A_h-A)+del_ty*(cos(E_h)*sin(alfa)*cos(A_h-A)+sin(E_h)*cos(alfa));
def_rx=-(2*incdnt_n*(C*Xm+del_nwx)+del_mtx)*cos(A_helio-A_h)-(2*incdnt_n*(C*Ym+del_nwy)+del_mty)*cos(E_h)*sin(A_helio-A_h);
def_ry= (2*incdnt_n*(C*Xm+del_nwx)+del_mtx)*cos(lambda)*sin(A_helio-A_h)-(2*incdnt_n*(C*Ym+del_nwy)+del_mty)*(cos(E_h)*cos(lambda)*cos(A_helio-A_h)-sin(E_h)*sin(lambda));

Xr= -Xm*cos(A_helio-A_h)-Ym*cos(E_h)*sin(A_helio-A_h)-S*tan(def_rx)+S0*tan(def_trx);
Yr= Xm*sin(A_helio-A_h)*cos(lambda)+Ym*(-cos(E_h)*cos(lambda)*cos(A_helio-A_h)+sin(E_h)*sin(lambda))-S*tan(def_ry)+S0*tan(def_try);

Xt0= -Xr*cos(A_helio)-Yr*sin(A_helio)*cos(lambda);
Yt0= -Xr*sin(A_helio)*cos(del_r)-Yr*(cos(A_helio)*cos(del_r)*cos(lambda)-sin(del_r)*sin(lambda));
Zt0= Xr*sin(A_helio)*sin(del_r)+Yr*(cos(A_helio)*sin(del_r)*cos(lambda)+cos(del_r)*sin(lambda));

alfat=-sin(def_rx)*cos(A_helio)+sin(def_ry)*sin(A_helio)*cos(lambda)-cos(sqrt((def_rx^2)+(def_ry^2)))*sin(A_helio)*sin(lambda);
betat=-sin(def_rx)*sin(A_helio)*cos(del_r)-sin(def_ry)*(cos(A_helio)*cos(del_r)*cos(lambda)-sin(del_r)*sin(lambda))+cos(sqrt((def_rx^2)+(def_ry^2)))*(cos(A_helio)*cos(del_r)*sin(lambda)+sin(del_r)*cos(lambda));
gamat=sin(def_rx)*sin(A_helio)*sin(del_r)+sin(def_ry)*(cos(A_helio)*sin(del_r)*cos(lambda)+cos(del_r)*sin(lambda))-cos(sqrt((def_rx^2)+(def_ry^2)))*(cos(A_helio)*sin(del_r)*sin(lambda)-cos(del_r)*cos(lambda));

xt= -(Zt0*aflat/gamat)+Xt0;
yt= -(Zt0*betat/gamat)+Yt0;
%}
%efficiency of field

%{
if(d_t_i<=1000)
atmos_atten_n = 0.99321-(0.0001176*d_t_i)+(d_t_i^2)*(1.97*10^(-8));
else
atmos_atten_n = exp(-0.0001106*d_t_i);
field_n = mirror_n*incdnt_n*atmos_atten_n;%*intercept_n*blocking_shading_n
end
%}