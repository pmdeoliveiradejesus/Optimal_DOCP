% ODOCR (67 PHASE) - Optimal Directional 67 PHASE Overcurrent Coordination Relays Problem
% Iterative model with transient configurations
% Sorrentino, Elmer, and José Vicente Rodríguez.
% "A novel and simpler way to include transient configurations in 
% optimal coordination of directional overcurrent protections."
% Electric Power Systems Research 180 (2020): 106127.
% 
% First version: July 15, 2020 
% Second version: July 15, 2022
% Third version: Oct. 21, 2024 - Database unification
% Fourth version: March. 31, 2026 - Coordination matrix corrected, new
% selectivity constraint Dt > kt*Tp
% Fifth version, August 36, 2026 loads in Zbus and IBR effect (norton) SCR
% Calculation

clear all
close all
clc
warning('off', 'all');
% disp('DOCR (67 Phase) with Transient Configurations')
disp('Optimal Coordination of Directional Overcurrent Relays ')
disp('in Interconnected Power Systems')
disp('Including Transient Network Configurations')
disp('Version 5.0 (c) 2026')
disp('Paulo M. De Oliveira pdeoliv at gmail.com')
disp('Power and Energy Group - https://power.uniandes.edu.co/')
disp('Universidad de los Andes, Colombia')
% disp(' ')
% disp(' ')
time0=cputime;
%% Case study
addpath('./data/')
global loadsflag Dmin Co upperbound lowerbound nr bdat ldat K Rf Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf  reversest linenumber back main dictiolines dictiorelays
%% Select test case
%case_threebus_optimization_data; %Urdaneta/Perez/Nadira Modified Test Case
case_eightbus_optimization_data; %Braga Test Case
%case_eightbus_optimization_data_IBR; %Braga Test Case with full IBR
%case_eightbus_optimization_data_IBR2; %Braga Test Case with full IBR + 1 SC
%case_eightbus_optimization_data_IBR3; %Braga Test Case with full IBR + 2 SC
%% ------------------------------------------------
nrf=100  ; % Define Number of (Uniformly Distributed) Relevant Faults Location   
time000=cputime;
nf=10  ;%number of faults per line (only for simulation, evaluation)
loadsflag=1;% '1' means loads are included in Ybs as shunt admittances '0' means no loads
reply2 = 'y';% 'n' means no prefault conditions
i=0;
econverg=1;
criterion=0; 
while econverg > 0
    i=i+1
%% Build the relay coordination model (B matrix) from the simulator logic
Case0=zeros(15,1);% Initialize type pairs vector
%reply2='n'; %no PREFAULT LOADFLOW
ncase=0; 
halfSm2x=[];
Bc=[];
B_case=[];
iter=1;%initialize counter for separation times
%% Begins Iterative process
for k=1:nrf
h=k/(nrf+1);% Fault location inserted here, with neval=1 h=0.5 0<h<1  %
if k==1
h=lowerbound;
end
if k==nrf
h=upperbound;
end
%Uniform distributed faults, if neval=1000 x goes from 0.001 to 0.999 Distance with respect to relay i
[halfSm2,S,Sx,Mx,Case,nlf,Co,Tix,Tq,index,nr,Iq,Ipback,Ipq,Ipp,Ii,Ij,Iip,Ijp,Ipi,Ipj,CaseType,distance,Bcoord,Bcoord_case]=run_classification_Bcoord(h,ncase,reply2,nlf,D,Ip,K);%Runs the classification script
% h
% pause
Case0=Case0+Case;%All 15 pair types classified are aggregated here
halfSm2x=[halfSm2;halfSm2x];
Bc=[Bc;Bcoord];%Coordination matrix
B_case=[B_case;Bcoord_case];% Coordination matrix with case identification 
end 
halfSm=halfSm2x*ktimes; % New constraint
Bextended=[Bc,B_case];% coordination matrix,last column indicates the case number
%% -------------------
% Build the relay coordination model (B matrix)    
 for kk2=1:nrf
     %% Objective function: only primary times of for near-end faults 
%% Just like the OF ussed by Ezzadine and Sorrentino
    for jj=1:nlf %generate 2 fault locations per line 1 to nlf
    Hnf(1,jj)=lowerbound;
    Hnf(2,jj)=upperbound;
    end  
for kk=1:2
[index3]=run_shortcircuit_optimizer(Hnf(kk,:),reply2);  %invoke the shortcircuit program    
jj=1;
for ii=1:nr
   jj=find(index3(:,3)==ii);
   fx2(kk,ii)=index3(jj(1),13);
end
end  
    for kk=1:2:nr  
     f(kk)=fx2(1,kk);
    end
    for kk=2:2:nr  
     f(kk)=fx2(2,kk);
    end     
end
%% Optimization model   
x0=zeros(nr,1);
% %% LP solver
 LB=ones(nr,1)*Dmin;
 UB=[]; 
 Aeq=[];
 beq=[];
B=Bc;%Coordinated B with simulator
bneq=ones(length(B(:,1)),1)*Co;
% disp('Optimizing. Please wait...! ')
% % Optimization problem: min f st. B > bneq
options = optimoptions('linprog','Algorithm','interior-point','display','off');
%% Optimization model with the new constraint
[Dsol,FVAL,EXITFLAG]=linprog(f,-B,-bneq,Aeq,beq,LB,UB,options); 
%% Optimization model with the new constraint
% Bnew=[B;B];% new constraint 
% bneqnew=[bneq;halfSm];%new constraint
% [Dsol,FVAL,EXITFLAG]=linprog(f,-Bnew,-bneqnew,Aeq,beq,LB,UB,options);
%min(Bnew*Dsol-bneqnew)
%pause
%% 
if EXITFLAG <= 0  
%disp('****************************************************')
fprintf('DO NOT CONVERGE')  
%disp('****************************************************')
%pause
else 
%[Bast2] = run_checker_transients(Dsol,nr,qmax,H);%coordination matrix obtained from D (optimized)
Case0=zeros(15,1);% Initialize type pairs vector
%reply2='n';
ncase=0; 
Bast2=[];
%% Begins Iterative process
for k=1:nrf
h=k/(nrf+1);% Fault location inserted here, with neval=1 h=0.5 0<h<1  %
if k==1
h=lowerbound;
end
if k==nrf
h=upperbound;
end
%Uniform distributed faults, if neval=1000 x goes from 0.001 to 0.999 Distance with respect to relay i
[Bcoord2]=run_classification_Verifier(h,ncase,reply2,nlf,Dsol,Ip,K);%Runs the classification script
%Case0=Case0+Case;%All 15 pair types classified are aggregated here
Bast2=[Bast2;Bcoord2];
end
econverg=abs(sum(sum(Bast2))-sum(sum(B)))% Consistency verification
FVAL;
if econverg==criterion
disp('Converged! Tripping sequence is correct!')
Bxx=B;
for k=1:nr
Dstr(i,k)=Dsol(k) ;
end
%pause
clear B
break 
break 
end
D=Dsol;
for k=1:nr
Dstr(i,k)=Dsol(k); 
end
end
end
elapsedtime=cputime-time000; 
%end
%%% investigating the quality of the solution
%Dsim=Dstr(i-1,:); % solution for the new constraint%
Dsim=Dstr(i,:); % solution for the base case
D=Dsim
%Rf=00/Zbase;
main_simulator_2026_internal% Run the simulator  

%Only for 8-bus case 
if nlf==7
dictiorelays=[1 2 3 4 5 6 7 8 9 10 11 12 13 14; 
               2 9 3 10 4 11 5 12 6 13 1 8 14 7];  % original case numbering 
   for kk=1:length(Dsim)  
   TAPs(dictiorelays(2,kk),1)= Dsim(kk); 
   end
   %% only IBR
   MWibr=[150 150]';%MW
   % for i=1:ngen
   % If(i)=1.4*MWibr(i)/(sqrt(3)*Vbase)*Ibase;%kA
   % end
 [Icc,mIcc,XR,Isources] =  run_SCRcalc(nlf,reply2);
 MVAcc=mIcc'*sqrt(3)*Vbase;
 for i=1:2
 SCRf(i)=MVAcc(i)/MWibr(i); 
 end
 SCRf % short circuit ratio
 mIcc = mIcc %short circuit in generating busbars in kA
 MVAcc=MVAcc %MVAcc in  generating busbars in MVA
 kIn=Isources/Ibase % currents from sources in per unit
 % Xginstall=[inf Xg(2) Xg(3) Xg(4) inf  Xg(6)];
 % for k=1:6
 %     Qc(k)=0.2*Vbase^2/(Xginstall(k)*Zbase);
 % end 
 % sum(Qc)
end