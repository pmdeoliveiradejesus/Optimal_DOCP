% ---------------------------------------
% DOCR - Dataset pm.deoliveira@uniandes.edu.co
% ---------------------------------------
% Optimal coordination of directional overcurrent relays in interconnected power systems, 
% IEEE Transactions on Power Delivery, Volume 3, No. 3, July 1988
%
%  Modified (L2 is only 15% length the original one)
%
%             N1   L1    N2
%        G1---|-R1-*--R2-|---G2
%             R6   N4    R3
%               \      /
%          L3 N6 *    * N5 L2
%                 \  /
%                 R5 R4
%                  -- N3
%                  |__ G3
%
% relays R i=1,2,3,4,5,6
% lines L k=1,2,3
% nodes N j=1,2,3,4,5,6
%         <-back-><-main-><-back-----><---main->
%      relayi faultk i  k  nsend nrec nsend nrec
% index=[	5	 1	 1	1	3      6    1	4;
% 	        4	 1	 2	1	3      5    2   4;
% 	        1	 2   3	2	1      4    2	5;
% 	        6	 2	 4	2	1      6    3	5;
% 	        3	 3	 5	3	2      5    3	6;
% 	        2	 3	 6	3	2      4    1	6];
% global Dmin Co  nr bdat ldat K Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf  
%% Relay data %%%%%%
ktimes=0.5;% Dt > ktimes * Ti
Co=.2;% allowable coordination interval (seconds)
%Relay curve settings, Standard Inverse (SI) IEC 60255
k11=.14;k12=0.02;k13=-1; SI=[k11;k12;k13];
%Relay curve settings, Very Inverse (VI) IEC 60255
k21=13.5;k22=1;k23=-1;VI=[k21;k22;k23];
%Relay curve settings, Extremely Inverse (EI) IEC 60255
k31=80;k32=2;k33=-1;EI=[k31;k32;k33];
K=[SI,SI,SI,SI,SI,SI];
%K=[VI,VI,VI,VI,VI,VI];
%K=[EI,EI,EI,EI,EI,EI];
%qmax=175; %relay polarization angle (deg) -to detect opposite currents
Ip(1)=.3;
Ip(2)=.06;
Ip(3)=.2;
Ip(4)=.24;
Ip(5)=.08;
Ip(6)=.2;
Dmin=0.1; % minimal dial setting
nr=length(Ip); 
dictiorelays=[1 2 3 4 5 6 ; 
              1 2 3 4 5 6];  % original case numbering 
dictiolines=[1 2 3 ; 
             1 2 3 ];  % original case numbering 
%% Simulation/Optimization Data
 lowerbound=0.001; % from % of the line
 upperbound=0.999; % to % of the line
%neval=1000;%number of faults per line (simulation, evaluation)
%% System Data
Sbase=100;%MVA
Vbase=69;%kV
Zbase=Vbase^2/Sbase;%ohms
Ibase=Sbase/(sqrt(3)*Vbase);%kA
econv=0.00000001;
itermax=100;
% Generator rated powers and reactances   
SG1=100;%MVA
SG2=25;%MVA
SG3=50;%MVA
XG1=.2;% pu
XG2=.12;% pu
XG3=.18;% pu
Xg1=(XG1*Vbase^2/SG1)/Zbase;
Xg2=(XG2*Vbase^2/SG2)/Zbase;
Xg3=(XG3*Vbase^2/SG3)/Zbase;
% Xg1=0.7;%pu
% Xg2=2.8;%pu
% Xg3=1.4;%pu

Rg1=0;
Rg2=0;
Rg3=0;
Rt=[0 0 0];% pu
Xt=[.00001 .00001 .00001];% pu
Xg=[Xg1 Xg2 Xg3];%*Zbaselow/Zbase;% pu
Rg=[Rg1 Rg2 Rg3];%pu
ngen=length(Xg);
% lines
R=[ 5.5   0.15*4.4  7.6]/Zbase; %line resistances (ohms)
X=[22.85  0.15*18   27  ]/Zbase; %line reactance (ohms)
B=[ 0     0    0  ]*Zbase; %line total susceptance (siemens)
nlf=length(R);
ldat=[  
        1 4 10001 Rg1  Xg1  0     0   0;
        2 5 10002 Rg2  Xg2  0     0   0;
        3 6 10003 Rg3  Xg3  0     0   0;        
        4 5 90001 R(1) X(1) B(1)  0   0;
        5 6 90002 R(2) X(2) B(2)  0   0;
        6 4 90003 R(3) X(3) B(3)  0   0;];
tdat=[ max(ldat(:,1))+1 ldat(1,2) 90031 Rt(1) Xt(1) 0  0   0;
       max(ldat(:,1))+2 ldat(2,2) 90031 Rt(2) Xt(2) 0  0   0;
       max(ldat(:,1))+3 ldat(3,2) 90036 Rt(3) Xt(3) 0  0   0];                       
bdat=[  1 1001  1 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0    1.0  1.0   0  0      0 0  0 0;
        2 1001  2 0.25 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        3 1001  2 0.50 0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        4 9001  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        5 9002  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        6 9003  3 0    0 0.3000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
                 ]; 
%npr=6;
Rf=0/Zbase;%fault resistance
%qmax=175; %relay polarization angle (deg) -to detect opposite currents
% Initial guess
% D=[0.112926061725303;
% 0.212089923140955;
% 0.165537455017500;
% 0.183868905459942;
% 0.161014736142870;
% 0.153996946311168];%no TC-STABLE
% D=[0.112926061723708
% 0.212089923137583
% 0.165537455016614
% 0.183868905458061
% 0.161014736141476
% 0.153996946307331];%TC-STABLE with new constraints
D=ones(6,1)*.1;
%Curve display parameters
%line 1 2 pairs
main(1)=1;
back(1)=5;
main(2)=2;
back(2)=4;
linenumber(1)=1;
linenumber(2)=1;
%line 2 2 pairs
main(3)=3;
back(3)=1;
main(4)=4;
back(4)=6;
linenumber(3)=2;
linenumber(4)=2;
%line 3 2 pairs
main(5)=5;
back(5)=3;
main(6)=6;
back(6)=2;
linenumber(5)=3;
linenumber(6)=3;
reversest=[3 5 1];

