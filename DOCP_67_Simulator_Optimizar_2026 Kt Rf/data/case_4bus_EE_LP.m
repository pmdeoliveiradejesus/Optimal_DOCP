% ODOCR - Dataset pm.deoliveira@uniandes.edu.co
% ---------------------------------------
global Dmin Co nr bdat ldat k1 k2 k3 Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf Dsol qmax lowerbound upperbound
%% Relay data %%%%%%
Co=.2;% allowable coordination interval (seconds)
%Relay curve settings, Standard Inverse (SI) IEC 60255
k11=.14;k12=0.02;k13=-1; SI=[k11;k12;k13];
%Relay curve settings, Very Inverse (VI) IEC 60255
k21=13.5;k22=1;k23=-1;VI=[k21;k22;k23];
%Relay curve settings, Extremely Inverse (EI) IEC 60255
k31=80;k32=2;k33=-1;EI=[k31;k32;k33];
    % K=[SI,SI,SI,SI,SI,SI,SI,SI,SI,SI];
    % D=[0.359815594808536%OPTIMAL SETTINGS% Normal Inverse
    %     0.221829245472319
    %     0.372135798927781
    %     0.357967599030791
    %     0.360795485818768
    %     0.372195295569904
    %     0.224537423491399
    %     0.325001566911574
    %     0.447254631382051
    %     0.365716246107707
    %     0.371775407560915
    %     0.440524030190349
    %     0.493242931600303
    %     0.501201352706540];
% % 
% K=[VI,VI,VI,VI,VI,VI,VI,VI,VI,VI;
% % Very Inverse
% D=[0.402667571091758%Optimal settings
%     0.159320270876163
%     0.558888554664143
%     0.623131351593373
%     0.625240383294610
%     0.557814853728637
%     0.160873627723016
%     0.323330405718382
%     0.681734490479448
%     0.367438263300778
%     0.372009379474232
%     0.677918708176419
%     0.756175939064922
%     0.765236739134183];
K=[EI,EI,EI,EI,EI,EI,EI,EI,EI,EI];
%Pickup currents (kA)    
Ip=[0.2;0.2;0.2;0.2;0.2;0.2;0.2;0.2;0.2;0.2];
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
SG1=100;%MVA
SG2=75;%MVA
SG3=8.5;%MVA
SG4=17;%MVA
XG1=.066666;% pu
XG2=.1;% pu
In3=SG3/(sqrt(3)*Vbase); %kA
In4=SG4/(sqrt(3)*Vbase); %kA
IF3=1.4*In3; %kA
IF4=1.4*In4; %kA
Xg1=(XG1*Vbase^2/SG1)/Zbase;
Xg2=(XG2*Vbase^2/SG2)/Zbase;
Xg3=Vbase/sqrt(3)/IF3/Zbase;% pu
Xg4=Vbase/sqrt(3)/IF4/Zbase;% pu
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
ldat=[  1 5 10001 Rg1  Xg1   0     0   0;
        2 6 10002 Rg2  Xg2   0     0   0;
        3 7 10003 Rg3  Xg3   0     0   0;  
        4 8 10004 Rg4  Xg4   0     0   0;
        5 6 90001 R(1) X(1) B(1)  0   0;
        6 7 90002 R(2) X(2) B(2)  0   0;
        7 8 90003 R(3) X(3) B(3)  0   0;  
        8 5 90004 R(4) X(4) B(4)  0   0;  
        5 7 90005 R(5) X(5) B(5)  0   0; 
 ];
tdat=[ max(ldat(:,1))+1 ldat(1,2) 90031 Rt(1) Xt(1) 0  0   0;
       max(ldat(:,1))+2 ldat(2,2) 90031 Rt(2) Xt(2) 0  0   0;
       max(ldat(:,1))+3 ldat(3,2) 90036 Rt(3) Xt(3) 0  0   0
       max(ldat(:,1))+4 ldat(4,2) 90036 Rt(4) Xt(4) 0  0   0];       
bdat=[  1 1001  1 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0    1.0  1.0   0  0      0 0  0 0;
        2 1001  2 0.25 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        3 1001  2 0.50 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        4 1001  2 0.15 0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        5 9002  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        6 9003  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        7 9003  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;       
        8 9003  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
                 ];     
qmax=175; %relay polarization angle (deg) -to detect opposite currents
D=[0.3353;0.1733;0.3376;0.1000;0.1286;0.5856;0.1000;0.6454;0.1374;0.1000];%no TC-STABLE
Ip=[0.2;0.2;0.25;0.25;0.2;0.3;0.3;0.2;0.2;0.25];

Rf=0/Zbase;%fault resistance ohms
ktimes=0.4;% Dt > ktimes * Ti   
% D=[0.3394
%     0.1748
%     0.2324
%     0.1000
%     0.1590
%     0.2822
%     0.1000
%     0.6673
%     0.1421
%     0.1000];

% D=[0.2500;
% 0.2500;
% 0.3769;
% 0.2847;
% 0.1700;
% 0.3347;
% 0.8873;
% 0.4124;
% 0.2888;
% 0.1000;
% 0.1070];%noIBR
