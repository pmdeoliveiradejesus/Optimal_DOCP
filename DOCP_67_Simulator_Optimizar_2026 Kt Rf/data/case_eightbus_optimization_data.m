% ODOCR - Dataset pm.deoliveira@uniandes.edu.co
% ---------------------------------------
% Braga, A. S., and J. Tome Saraiva. 
% "Coordination of overcurrent directional relays in meshed networks using the Simplex method." 
% Proceedings of 8th Mediterranean Electrotechnical Conference on Industrial Applications in Power Systems, 
% Computer Science and Telecommunications (MELECON 96). Vol. 3. IEEE, 1996.
% 8-bus test system 7-lines 2-generator/transformers
%
%            G1
%            |
%          ----- N14
%            |
%            T1
%            |           
%       N1 -----  N2 -----  N3 ---- 
%           |||  L1   | |  L2  | |
%          12||1__N7_2| |3_N8_4| 5
%           | \___N13__          |
%          N12         \14       N9
%        L6 |  __N11__ | _N10__  |L3
%          11 10  L5  9|8  L4  7 6
%           | |       |||      | |
%       N6 -----  N5 ----- N4 ----- 
%                     |
%                     T2
%                     |
%                   ----- N15
%                     |
%                     G2
% lines K k=1,2,3,4,5,6,7
% nodes B b=1,...,14
% relays R i=1,2,3,4,5,6,7,8,9,10,11,12,13,14
%        Backup	Primary	Backup	Primary		
%        Ri	Lk	Ri	Lk	Bs	Br	Bs  Br	
% Index=[11	1	1	1	6	12	1	7	;
%        14	1	1	1	5	13	1	7	;
%        1	2	3	2	1	7	2	8	;
%        3	3	5	3	2	8	3	9	;
%        5	4	7	4	3	9	4	10	;
%        7	5	9	5	4	10	5	11	;
%        13	5	9	5	1	13	5	11	;
%        9	6	11	6	5	11	6	12	;
%        2	7	13	7	2	7	1	13	;
%        11	7	13	7	6	12	1	13	;
%        4	1	2	1	3	8	2	7	;
%        6	2	4	2	4	9	3	8	;
%        8	3	6	3	5	10	4	9	;
%        10	4	8	4	6	11	5	10	;
%        13	4	8	4	1	13	5	10	;
%        12	5	10	5	1	12	6	11	;
%        2	6	12	6	2	7	1	12	;
%        14	6	12	6	5	13	1	12	;
%        7	7	14	7	4	10	5	13	;
%        10	7	14	7	6	11	5	13	];
%global Rf Dmin Co nr bdat ldat K Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf Dsol qmax lowerbound upperbound npr reversest linenumber back main dictiolines dictiorelays
%% Relay data 

Co=.3;% allowable coordination interval (seconds)
%Relay curve settings, Standard Inverse (SI) IEC 60255
k11=.14;k12=0.02;k13=-1; SI=[k11;k12;k13];
%Relay curve settings, Very Inverse (VI) IEC 60255
k21=13.5;k22=1;k23=-1;VI=[k21;k22;k23];
%Relay curve settings, Extremely Inverse (EI) IEC 60255
k31=80;k32=2;k33=-1;EI=[k31;k32;k33];
% K=[SI,SI,SI,SI,SI,SI,EI,SI,SI,SI,SI,SI,SI,SI];
% K=[VI,VI,VI,VI,VI,VI,VI,VI,VI,VI,VI,VI,VI,VI];
  K=[EI,EI,EI,EI,EI,EI,EI,EI,EI,EI,EI,EI,EI,EI];
%Pickup currents (kA)    
Ip=[0.5;0.12;0.24;0.12;0.12;0.24;0.12;0.6;0.18;0.12	;0.12;0.18;0.12;0.12];% 
Dmin=0.05; % minimal dial setting
nr=length(Ip); 
dictionodes=[1 2 3 4 5 6 7 8   9 10 11 12 13 14 15; 
             1 3 4 5 6 2 9 10 11 12 13 14 15  7 8];   % original case numbering     
dictiorelays=[1 2 3 4 5 6 7 8 9 10 11 12 13 14; 
              2 9 3 10 4 11 5 12 6 13 1 8 14 7];  % original case numbering  
dictiolines=[1 2 3 4 5 6 7 ; 
             2 3 4 5 6 1 7 ];  % original case numbering  
%% Simulation/Optimization Data
 lowerbound=0.00000001; % from % of the lline
 upperbound=0.99999999; % to % of the line
%% System Data
Sbase=150;%MVA
Vbase=150;%kV
Vbaselow=10;%kV
Zbase=Vbase^2/Sbase;%ohms
Zbaselow=Vbaselow^2/Sbase;%ohms
Ibase=Sbase/(sqrt(3)*Vbase);%kA
econv=0.00000001;
itermax=100;
% Generator rated powers and reactances    
Rt=[0 0];% pu
Xt=[.04 .04];% pu
Xg=[.15 .15];%*Zbaselow/Zbase;% pu
Rg=[0 0];%pu
ngen=length(Xg);
% lines
Le=[100 70 80 100 110 90 100]; 
R=[.004 .0057 .005 .005 .0045 .0044 .005].*Le/Zbase; %line resistances (ohms)
X=[ .05 .0714 .0563 .045 .0409 .05 .05].*Le/Zbase; %line reactance (ohms)
B=[0 0 0 0 0 0 0 0]*Zbase; %line total susceptance (siemens)
nlf=length(R);
ldat=[  1 3 10001 Rg(1) Xg(1) 0   0   0;
        2 7 10002 Rg(2) Xg(2) 0   0   0;
        3 4 90002 R(2) X(2) B(2)  0   0;
        4 5 90003 R(3) X(3) B(3)  0   0;
        5 6 90004 R(4) X(4) B(4)  0   0;        
        6 7 90005 R(5) X(5) B(5)  0   0;
        7 8 90005 R(6) X(6) B(6)  0   0; 
        8 3 90001 R(1) X(1) B(1)  0   0;
        3 7 90006 R(7) X(7) B(7)  0   0;  ];
tdat=[ max(ldat(:,1))+1 ldat(1,2) 90031 Rt(1) Xt(1) 0  0   0;
       max(ldat(:,1))+2 ldat(2,2) 90036 Rt(2) Xt(2) 0  0   0];
bdat=[  1 1001  1 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0    1.0  1.0   0  0      0 0  0 0;
        2 1001  2 150/150 0 0.0000 0.0000 1 0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        3 9001  3 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        4 9002  3 0    0 60/150 40/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        5 9003  3 0    0 70/150 40/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        6 9001  3 0    0 70/150 50/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        7 9002  3 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        8 9003  3 0    0 40/150 20/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;]; 
%npr=20;

%relay polarization angle (deg) -45 < angle(V.conj(I))<+135 deg -to detect opposite currents
tmax=45+360;%fault angle current tripping zone
tmin=225;%fault angle current tripping zone 
 
Rf=0/Zbase;%fault resistance ohms
ktimes=0.5;% Dt > ktimes * Ti   
D=0.10*ones(14,1);
% D=[0.666844297235227
% 0.140129428818729
% 1.50265594040710
% 2.33099070685522
% 2.26484506308189
% 1.52877072699728
% 0.145046517278131
% 0.465510931143944
% 1.87044662304501
% 0.564392990736045
% 0.545469828380725
% 1.88737790557358
% 1.75677146202653
% 1.91468006786445];
% Plotting data
%line 1 3 pairs
main(1)=1;
back(1)=14;
main(2)=1;
back(2)=11;
main(3)=2;
back(3)=4;
linenumber(1)=1;
linenumber(2)=1;
linenumber(3)=1;
%line 2 2 pairs
main(4)=3;
back(4)=1;
main(5)=4;
back(5)=6;
linenumber(4)=2;
linenumber(5)=2;
%line 3 2 pairs
main(6)=5;
back(6)=3;
main(7)=6;
back(7)=8;
linenumber(6)=3;
linenumber(7)=3;
%line 4  3 pairs
main(8)=7;
back(8)=5;
main(9)=8;
back(9)=10;
main(10)=8;
back(10)=13;
linenumber(8)=4;
linenumber(9)=4;
linenumber(10)=4;
%line 5 3pairs
main(11)=9;
back(11)=7;
main(12)=9;
back(12)=13;
main(13)=10;
back(13)=12;
linenumber(11)=5;
linenumber(12)=5;
linenumber(13)=5;
%line 6 3pairs
main(14)=11;
back(14)=9;
main(15)=12;
back(15)=2;
main(16)=12;
back(16)=14;
linenumber(14)=6;
linenumber(15)=6;
linenumber(16)=6;
%line 7 4 pairs
main(17)=13;
back(17)=11;
main(18)=13;
back(18)=2;
main(19)=14;
back(19)=7;
main(20)=14;
back(20)=10;
linenumber(17)=7;
linenumber(18)=7;
linenumber(19)=7;
linenumber(20)=7;
reversest=[3 5 7 9 10 13 15 16 19 20];


