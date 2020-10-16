%% Ball and Plate Constants
%Ball specs                         
r_b = 19.05*(10^-3);    %radius of the rubber ball(m)> diameter = 38.1 mm
Vb = (4*pi*(r_b^3))/3;                  %volume of the ball
density_ball = 1522;                  %density of a rubber ball
m_b = Vb*density_ball;                %mass of the ball
Ib = (2/5)*m_b*r_b^2;                        % inertia of the ball (Kg.m^2)

K =  (m_b/(m_b+(Ib/r_b^2)));                 %constant (m/m+(Ib/r^2)) = 5/7
g = 9.81;  

%Plate specs
W  = 0.34;                                   %width plate(m)
L = 0.34;                                    %Length plate(m)
H = 0.005;                                  %height plate(m)
density = 1180;                         %Perspex_density (kg/m^3)
m_p = W*L*H*density ;                  %mass of the plate = volume*density
Ip = (1/12)*m_p*(L^2+W^2);                %inertia of the plate (Kg.m^2)

%Motor Constants
Kt = 0.055;                           %Torque Costant
Ra =0.08;                                   %Armature resistance
J = 10^-3;                                  %Rotor Inertia
f = 10^-8;                                  %Rotor Damping
Kb = 73*10^-6;                              %Back-emf Constant
La = 0;

%TF's:
G_px = tf(K*g,[1 0 0]); 
%% **** linear plant model*****
%
%  X(s)      Km/m
%  --- = -------------
%  I(s)    s(s+Cv/m)
%
% with x in [m] and i in [A]

%m=0.139;
%Km=11;
%Cv=2.5;
m=0.118+(1/3)*m_p;   %slider mass + 1/3 * Plate mass
Km=11;
Cv=16.5;

Ts=1/2000;

% **** non-linear aspects ********
% Commutation 

% distance in [m] between center of two adjacent magnets, pole pitch
TAU=0.01;

CurAngOffset=-0.00137*pi/TAU;



% system bandwidth
BW=Cv/m;

Hp=zpk([],[0 -Cv/m],Km/m);
[nump,denp]=tfdata(Hp,'v');
Hp_mi = zpk([],[-Cv/m],Km/m);
[nump_mi,denp_mi] = tfdata(Hp_mi,'v');
Kp=Km/m;
p=-Cv/m;
Ts_Outer= 0.1;

% controller BW 100 rad/s
% Hc=zpk([-100/3 -100/6],[-3*100 -6*100 0],330000);
% [numc,denc]=tfdata(Hc,'v');

% controller BW 200 rad/s (32 Hz)

B=300;
% choose the controller gain such that the open loop magnitude is 0dB at w=B
% the plant magnitude at w=B is:
P_magn=Kp/(B*sqrt(p^2+B^2));
% the controller magnitude at w=B is:
C_magn=1/(B*18);
%The required controller gain to make the open loop magnitude 0dB at w=B must then be:
Kc=1/(P_magn*C_magn);


Hc=zpk([-B/3 -B/6],[-3*B -6*B 0],Kc);
[numc,denc]=tfdata(Hc,'v');

% same controller minus the integrator for dSPACE implementation.
% integrator is added as separate block, where the saturation limits can be
% set to avoid integrator wind-up.

% Hc_mi=zpk([-100/3 -100/6],[-3*100 -6*100],330000);
% [numc_mi,denc_mi]=tfdata(Hc_mi,'v');

Hc_mi=zpk([-B/3 -B/6],[-3*B -6*B],Kc);
[numc_mi,denc_mi]=tfdata(Hc_mi,'v');

% m=0.367;
% Km=11;
% Cv=16.5;
% Hpreal=zpk([],[0 -Cv/m],Km/m);
% s = tf('s');
% Cont = 54*((s/(2*pi*10/3)+1)*(s+(2*pi*10/5)))/(s*(s/(2*pi*10*3)+1)*(s/(2*pi*1000)+1));
% Cont = 19*((s/(2*pi*5/3)+1)*(s+2*pi*5/5))/(s*(s/(2*pi*5*3)+1)*(s/(2*pi*1000)+1));
% Cont_mi = 50*((s/(2*pi*10/3)+1))/((s/(2*pi*10*3)+1));
% [numc_b,denc_b]=tfdata(Cont_mi,'v');

% Hc=zpk([-100/3],[-3*100 -6*100],330*7/6);
% [numc,denc]=tfdata(Hc,'v');

% same controller minus the integrator for dSPACE implementation.
% integrator is added as separate block, where the saturation limits can be
% set to avoid integrator wind-up.

% Hc_mi=zpk([-100/3],[-3*100 -6*100],330*7/6);
% [numc_mi,denc_mi]=tfdata(Hc_mi,'v');

%rltool(Hp,Hc)

%% %% Based on above calculations, below the tf of the plant and the controller 

Km = 1;
Tau_m = 0.005175;
G_m1=tf(Km,[Tau_m 1]);   % a similar TF of the motor eith same settling time
Plant = G_m1*G_px;   %Total Plant
Ts_Outer = 0.1;
z=tf('z',Ts_Outer);
% Case 2: design total system in discrete: 
Plant_dis = c2d(Plant,Ts_Outer,'ZOH'); %convert plant to discrete
[num_plant,den_plant] = tfdata(Plant_dis,'v');
zpk(Plant_dis);                         %gives the function in order for poles and zeroes 

       
C12_dis_designed = 1.9*(z+0.8546)*(z^2-1.974*z+0.974)/((z-1)*(z-0.5504)*(z+0.9429)); % chosen one
[numco12_dis_des,denco12_dis_des,Ts_Outer] = tfdata(C12_dis_designed,'v');
 
C13_dis_designed = 2.5037*(z-0.06981)*(z^2-1.92*z+0.9213)/((z-1)*(z-0.3774)*(z^2+0.3764*z+0.03567)); 
[numco13_dis_des,denco13_dis_des,Ts_Outer] = tfdata(C13_dis_designed,'v');

%% filter before actuators circle
F = tf(1,[(1/(2*pi*1)) 1]);
[numF,denF] = tfdata(F,'v'); 
zpk(F);
