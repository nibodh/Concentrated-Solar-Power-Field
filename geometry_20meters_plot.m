clear;
clc;

%inputs
slope = 0; %terrain slope
twrht = 10; %tower height
ht=twrht;
lnth=1; %length of mirror
wdth=1; %width of mirror
i_l=lnth;
i_w=wdth;
h_i=0.75; %height of heliostat
fa=0.98; %ratio of reflecting surface to total surface
%ratio of heliostat separation distance to heliostat length is dS
lr=0; %receiver height
Rmax=30; %maximum ring radius in the field
half_azi_angle_field=45; %maximum half azimuthal angle of the field
Vmax=half_azi_angle_field*pi/180; %maximum angular direction

%assigning intermediate variables
f=i_w/i_l; %ratio of width to length
dS= (2*f)-sqrt(1+(f*f)); %MINIMUM ratio of safety distance to length
Am=f*fa*i_l*i_l; %reflecting area on each mirror
Rmin=ht;   %in paper, first row at R = ht of tower
zl=ht-(0.5*lr);
rm=i_l/2;

DM = max((i_l*((sqrt(1+(f*f)))+dS)),(2*i_w)); %distance between adjoining mirrors in first row
del_rmin=DM*(3^0.5)*cos(slope)/2; %minimal radial increment from geometrical overview of circles
i=1; %ring counter as 1,2,3,...
j=1; %angular counter as 1,2,3...
d_i=Rmin;
vj=DM/(2*d_i); %HALF of angular separation between heliostats

t=i+1; 
s=j+1;
fill(s,1)=j;
fill(s,2)=vj*180/pi;

fill(t,3)=i;
fill(t,4)=d_i;
fill(t,5)=(2*(round(Vmax/(2*vj)))+1); %number of heliostats
fill(t,6)=d_i*2*vj;

while(d_i<=Rmax)
    i=i+1;
    d_i_1=(d_i*cos(vj))+sqrt(((DM*cos(slope))^2)+((d_i*sin(vj))^2)); %calculating second ring radius from overview
    while(i~=1)&&(d_i<=Rmax) %will return one ring more than the Rmax
        
        zm=h_i+(d_i*tan(slope));
        a = (zl^2)*((rm^2)-(d_i^2));
        b = 2*d_i*zl*(zl-zm);
        c = (rm^2)-((zm-zl)^2);
        d_i_0 = (-b+sqrt((b^2)-(4*a*c)))/(a*2);
        
        A = -(((2*zl*d_i_0)+tan(slope))*tan(slope)+(zl*zl*d_i_0*d_i_0));
        B = 2*(zl-h_i)*((zl*d_i_0)+tan(slope));
        C = ((rm^2)*(1+(zl*zl*d_i_0*d_i_0)))-((zl-h_i)^2);
        d_i_11 = max((-B-sqrt((B^2)-(4*A*C)))/(A*2),(d_i_1+del_rmin));
        
        if mod(i,2)==0
            Nm=2*(round(Vmax/(2*vj)))+1; %number of heliostats in essential rings
        else
            Nm=2*(round((Vmax-vj)/(2*vj)))+2; %number of heliostats in staggered rings
        end
        Af=Vmax*(((d_i_11+(0.5*DM))^2)-((d_i_1+(0.5*DM))^2)); %land area covered
        del=Am*Nm/Af; %mirror density= number of mirrors in each ring times  area of each mirror/area of land covered
        
        zm=h_i+(d_i_1*tan(slope));
        a = (zl^2)*((rm^2)-(d_i_1^2));
        b = 2*d_i_1*zl*(zl-zm);
        c = (rm^2)-((zm-zl)^2);
        d_i_0 = (-b+sqrt((b^2)-(4*a*c)))/(a*2);
        
        A = -(((2*zl*d_i_0)+tan(slope))*tan(slope)+(zl*zl*d_i_0*d_i_0));
        B = 2*(zl-h_i)*((zl*d_i_0)+tan(slope));
        C = ((rm^2)*(1+(zl*zl*d_i_0*d_i_0)))-((zl-h_i)^2);
        d_i_12 = max((-B-sqrt((B^2)-(4*A*C)))/(A*2),(d_i_11+del_rmin));
        
        vj_1=DM/(2*d_i_12);
        Nm1=2*(round(Vmax/(2*vj_1)))+1;
        Af1=Vmax*(((d_i_12+(0.5*DM))^2)-((d_i_1+(0.5*DM))^2));
        del1=Am*Nm1/Af1;
        
        if(del>=del1)
            t=t+1;
            fill(t,1)=j; %group number
            fill(t,3)=i; %row in group
            fill(t,2)=vj*180/pi; %HALF of angular separation in degrees
            fill(t,4)=d_i_1; %radial distance from foot of tower
            fill(t,6)=d_i_1*2*vj; %separation between adjacent heliostats
            fill(t,7)=fill(t,4)-fill(t-1,4); %distance from previous row
            if mod(i,2)==0
                fill(t,5)=(2*(round((Vmax-vj)/(2*vj)))+2); %number of heliostats in stagg rings
            else
                fill(t,5)=(2*(round(Vmax/(2*vj)))+1); %number of heliostats in ess rings
            end
            i=i+1;
            d_i=d_i_1;
            d_i_1=d_i_11;
        else
            NRGj=i-1; %number of groups
            j=j+1;
            vj=vj_1;
            i=1;
            t=t+1;
            s=t;
            fill(s,1)=j;
            fill(s,2)=vj*180/pi;
            fill(t,4)=d_i_12;
            fill(t,5)=Nm1;
            fill(t,6)=d_i_12*2*vj;
            fill(t,3)=i;
            fill(t,7)=fill(t,4)-fill(t-1,4); %distance from previous row
            d_i=d_i_12;
        end
    end
end

if(i~=1)
    NG=j; %number of ring in the group
    NRGj=i-2;
else
    NG=j-1;
end

fill(t,1)=0;
fill(t,2)=0;
fill(t,3)=0;
fill(t,4)=0;
fill(t,5)=0;
fill(t,6)=0;
fill(t,7)=0;

save radial.mat fill;

%for x,y,z coordinates
a=1;
p=0;
while (a<t-1)
    a=a+1;
    nmax=round(Vmax/(fill(a,2)*pi/180));
        if mod(fill(a,3),2)~=0
            for n=0:2:nmax
                p=p+1;
                Vm=n*(fill(a,2)*pi/180);
                coord(p,1)=fill(a,4)*sin(Vm);
                coord(p,2)=fill(a,4)*cos(Vm);
                [coord(p,4),coord(p,5)] = elevation_azimuthal(fill(a,4),(180-(n*(fill(a,2)))));
                %coord(p,3)=h_i+fill(a,4)*tan(slope);
                if (n~=0)
                    p=p+1;
                    coord(p,1)=-coord(p-1,1);
                    coord(p,2)=coord(p-1,2);
                    [coord(p,4),coord(p,5)] = elevation_azimuthal(fill(a,4),(180+(n*(fill(a,2)))));
                    %coord(p,3)=coord(p,3);
                end
            end
        else
            for n=1:2:nmax
                p=p+1;
                Vm=n*(fill(a,2)*pi/180);
                coord(p,1)=fill(a,4)*sin(Vm);
                coord(p,2)=fill(a,4)*cos(Vm);
                [coord(p,4),coord(p,5)] = elevation_azimuthal(fill(a,4),(180-(n*(fill(a,2)))));
                %coord(p,3)=h_i+fill(a,4)*tan(slope);
                p=p+1;
                    coord(p,1)=-coord(p-1,1);
                    coord(p,2)=coord(p-1,2);
                    [coord(p,4),coord(p,5)] = elevation_azimuthal(fill(a,4),(180+(n*(fill(a,2)))));
                    %coord(p,3)=coord(p,3);
            end
        end
end

save cartesian.mat coord;

a=1;
deld_min=fill(3,4)-fill(2,4);
deld_max=fill(3,4)-fill(2,4);
d_min=fill(2,6);
d_max=fill(2,6);
while (a<t-2)
    a=a+1;
    deld_min=min(deld_min,(fill(a+1,4)-fill(a,4)));
    deld_max=max(deld_max,(fill(a+1,4)-fill(a,4)));
    d_min=min(fill(a,6),d_min);
    d_max=max(fill(a,6),d_max);
end
limits(1,1)=deld_min; %max and min radial distances
limits(1,2)=deld_max;
limits(2,1)=d_min; %min and max separation between adjacent heliostats
limits(2,2)=d_max;

a=1;
E_min=coord(a,5);
E_max=coord(a,5);
A_min=coord(a,4);
A_max=coord(a,4);
while (a<p)
    a=a+1;
    E_min=min(coord(a,5),E_min);
    E_max=max(coord(a,5),E_max);
    A_min=min(coord(a,4),A_min);
    A_max=max(coord(a,4),A_max);
end
limits(3,1)=A_min; %max and min azimuthal angles
limits(3,2)=A_max;
limits(4,1)=E_min; %max and min elevation angles
limits(4,2)=E_max;

save test_limits.mat limits;

plot(0,0,'o',coord(:,1),coord(:,2),'o',0,-0.25,0,30.25,20.25,0,-20.25,0);
xlabel('meters');
ylabel('meters');
title('overview');