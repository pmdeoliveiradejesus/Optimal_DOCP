% ODOCR - Dataset pm.deoliveira@uniandes.edu.co
% ---------------------------------------
global Dmin Co nr bdat ldat k1 k2 k3 Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf Dsol qmax lowerbound upperbound
%% Relay data %%%%%%
% relays R i=1,2,3,4,5,6
% lines L k=1,2,3
% nodes N j=1,2,3,4,5,6
%<-back-><-main-><-back-----><---main->
%relayi faultk i  k  nsend nrec nsend nrec
% 7	1	1	1	4	8	1
% 1	2	3	2	1	5	2
% 10	2	3	2	4	9	2
% 3	3	5	3	2	6	3
% 5	4	7	4	3	7	4
% 9	4	7	4	2	9	4
% 1	5	9	5	1	5	2
% 4	5	9	5	3	6	2
% 4	1	2	1	3	6	2
% 10	1	2	1	4	9	2
% 6	2	4	2	4	7	3
% 8	3	6	3	1	8	4
% 9	3	6	3	2	9	4
% 2	4	8	4	2	5	1
% 5	5	10	5	3	7	4
% 8	5	10	5	1	8	4


Co=.2;% allowable coordination interval (seconds)
%Relay curve settings, Standard Inverse (SI) IEC 60255
%k1=0.14;k2=0.02;k3=-1;
%Relay curve settings, Very Inverse (VI) IEC 60255
%k1=13.5;k2=1;k3=-1;
%Relay curve settings, Extremely Inverse (EI) IEC 60255
k1=80;k2=2;k3=-1;
%Pickup currents (kA)    
Ip=[0.2;0.2;0.25;0.25;0.2;0.3;0.3;0.2;0.2;0.25];
Dmin=0.1; % minimal dial setting
nr=length(Ip); 
%% Simulation/Optimization Data
lowerbound=.00000001; % from % of the lline
upperbound=.99999999; % to % of the line
%% System Data
Sbase=100;%MVA
Vbase=115;%kV
Zbase=Vbase^2/Sbase;%ohms
Ibase=Sbase/(sqrt(3)*Vbase);%kA
econv=0.00000001;
itermax=100;
% Generator rated powers and reactances   
SG1=150;%MVA
SG2=75;%MVA
SG3=50;%MVA
SG4=100;%MVA
XG1=.1;% pu
 XG2=.1;% pu
% XG3=.15;% pu
% XG4=.18;% pu
XG3=3;% pu
XG4=3;% pu
Xg1=(XG1*Vbase^2/SG1)/Zbase;
 Xg2=(XG2*Vbase^2/SG2)/Zbase;
 Xg3=(XG3*Vbase^2/SG3)/Zbase;
 Xg4=(XG4*Vbase^2/SG4)/Zbase;
Rg1=0;
Rg2=0;
Rg3=0;
Rg4=0;
Rt=[0 0 0 0];% pu
Xt=[.035 .045 .04 0.06];% pu
Xg=[Xg1 Xg2 Xg3 Xg4];%*Zbaselow/Zbase;% pu
Rg=[Rg1 Rg2 Rg3 Rg4];%pu
ngen=length(Xg);
% lines
Re=0.05;%Ohm/km
Xl=0.4;%Ohm/km
Length=[20 25 8 20 35];%km
R=Re*Length.*[1 1 1 1 1 ]/Zbase; %line resistances (ohms)
X=Xl*Length.*[1 1 1 1 1 ]/Zbase; %line reactances (ohms)
B=[0 0 0 0 0 ]*Zbase; %line total susceptance (siemens)
nlf=length(R);
ldat=[  
        1 5 10001 Rg1  Xg1   0     0   0;
        2 6 10002 Rg2  Xg2   0     0   0;
        3 7 10003 Rg3  Xg3   0     0   0;  
        4 8 10004 Rg4  Xg4   0     0   0;
        5 6 90001 R(1) X(1) B(1)  0   0;
        6 7 90002 R(2) X(2) B(2)  0   0;
        7 8 90003 R(3) X(3) B(3)  0   0;  
        8 5 90004 R(4) X(4) B(4)  0   0;  
        6 8 90005 R(5) X(5) B(5)  0   0; 
 ];
tdat=[ max(ldat(:,1))+1 ldat(1,2) 90031 Rt(1) Xt(1) 0  0   0;
       max(ldat(:,1))+2 ldat(2,2) 90031 Rt(2) Xt(2) 0  0   0;
       max(ldat(:,1))+3 ldat(3,2) 90036 Rt(3) Xt(3) 0  0   0;
       max(ldat(:,1))+4 ldat(4,2) 90036 Rt(4) Xt(4) 0  0   0];          
bdat=[  1 1001  1 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0    1.0  1.0   0  0      0 0  0 0;
        2 1001  2 0.75 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        3 1001  2 1.00 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        4 1001  2 0.50 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        5 9002  3 0    0 0.7500 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        6 9003  3 0    0 0.7500 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        7 9003  3 0    0 0.7500 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;       
        8 9003  3 0    0 0.7500 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
                 ];     
qmax=175; %relay polarization angle (deg) -to detect opposite currents

D=[0.3394;
    0.1748;
    0.2324;
    0.1000;
    0.1590;
    0.2822;
    0.1000;
    0.6673;
    0.1421;
    0.1000];%no TC-STABLE


D2=[0.2500;
0.3769;
0.2847;
0.1700;
0.3347;
0.8873;
0.4124;
0.2888;
0.1000;
0.1070;];%no TC-STABLE - No IBR


D3=[ 0.4695;
    0.1949;
    0.4988;
    0.1000;
    0.2181;
    0.3168;
    0.1000;
    0.7812;
    0.1941;
    0.1000];%TC-STABLE





    D=[0.4608;
    0.1797;
    0.5129;
    0.1000;
    0.2217;
    0.3029;
    0.1000;
    0.7917;
    0.1966;
    0.1000];%TC-STABLE (no IBR) must chand Xg3 and Xg4 = 30000




% Ibase*(inv(Xg1+Xt(1)))
% Ibase*(inv(Xg2+Xt(2)))
% sqrt(3)*115*Ibase*(inv(Xg3+Xt(3)))/1.4
% sqrt(3)*115*Ibase*(inv(Xg4+Xt(4)))/1.4