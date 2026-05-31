
warning('off', 'all');
time000=cputime;
%% Initial screen
% disp('*******************************************************')
 disp('DOCR (67 Phase)                                        ')
 disp('67 Relay System Graphical Simulator                    ')
 disp('Version 1.0 (c) 2022                                   ')
 disp('Power and Energy Group - https://power.uniandes.edu.co/')
 disp('Universidad de los Andes, Colombia')
 disp('*******************************************************')
addpath('./data/')
global Dmin Rf ktimes Co nr bdat ldat K Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf qmax lowerbound upperbound npr reversest linenumber back main dictiolines dictiorelays
clear Case0
Case0=zeros(15,1);% Initialize type pairs vector
B=[];%Initialize the coordination matrix
B_case=[];%Initialize the coordination vector case   
T=[0];%Initialize primary times
i=1;%Initialize flag for time progress
iter=1;%initialize counter for separation times
iterk=1;%initialize counter short circuit currents
%% Begins Iterative process
B22=[];Tppal=[];
for k=1:nf
i=i+1;
if i==nf/4
fprintf('Wait! simulating... Progress: 25%%')
elseif i==nf/2
fprintf(', 50%%')
elseif i==3*nf/4
fprintf(', 75%%')
elseif i==nf
fprintf(', 100%%\n')
end
h=k/(nf+1);% Fault location inserted here, with neval=1 h=0.5 0<h<1  %
if k==1
h=lowerbound;
end
if k==nf
h=upperbound;
end
%Uniform distributed faults, if neval=1000 x goes from 0.001 to 0.999 Distance with respect to relay i
[Tppalx, S,Sx,Mx,Case,nlf,Co,Tix,Tq,index,nr,Iq,Ipback,Ipq,Ipp,Ii,Ij,Iip,Ijp,Ipi,Ipj,CaseType,distance,Bcoord3,Bcoord_case3]=run_classification4(h,ncase,reply2,nlf,Dsim,Ip,K);
Tppal=[Tppal;Tppalx];
Case0=Case0+Case;%All 15 pair types classified are aggregated here
B22=[B22;Bcoord3];%Coordination matrix
%B_case=[B_case;Bcoord_case];%Coordination matrix with case identification 
T=unique(vertcat(unique([Tix';Tq']),T));%All nr primary times are aggregated here
T(T==0)=[];
sizeS(k)=length(S);%Set length of each separation vector
for k=1:length(S)
SepTime(iter,1)=S(k);
SepTime(iter,2)=Sx(k,1);
SepTime(iter,3)=Sx(k,2);
SepTime(iter,4)=Sx(k,3);
SepTime(iter,5)=Sx(k,4);
SepTime(iter,6)=Sx(k,5);
SepTime(iter,7)=Sx(k,6);%Tip

iter=iter+1;
end%All calculated separation times are aggregated here
for k=1:length(Iq)
    iterk=iterk+1;
    ShortC(iterk,1)=Iq(k);
ShortC(iterk,2)=Mx(k,1);
ShortC(iterk,3)=Mx(k,2);
ShortC(iterk,4)=Mx(k,3);
ShortC(iterk,5)=Mx(k,4);
ShortC(iterk,6)=Mx(k,5);
ShortC(iterk,7)=Ipback(k);
ShortC(iterk,8)=Ii(k);
ShortC(iterk,9)=Ij(k);
ShortC(iterk,10)=Iip(k);
ShortC(iterk,11)=Ijp(k);
ShortC(iterk,12)=Ipp(k);
ShortC(iterk,13)=Ipi(k);
ShortC(iterk,14)=Ipj(k);
% ShortC(iterk,15)=Ipq(k);
end%
end
%% Iterative process ends
AAA=B22*Dsim'-Tppal*ktimes;
for k=length(AAA(:,1))/2+1:length(AAA(:,1))/2
if AAA(k,1) < 0 
kki=kki+1;
end
end
%determine number of separation times below specified Co (CTI)
Bextended=[B,B_case];% coordination matrix,last column indicates the case number
elapsedtime000=cputime-time000;% Set simulation time
ki=0;
kki=0;
for k=1:length(SepTime(:,1))
if SepTime(k,1) < Co
ki=ki+1;
end
end%determine number of separation times below specified Co (CTI)

% Types 1 to 6 calculation
result(1,1)= Case0(1); %Number of relay pairs Type 1
result(2,1)= Case0(2); %Number of relay pairs Type 2
result(3,1)= Case0(3)+Case0(4)+Case0(5); %Number of relay pairs Type 3
result(4,1)= Case0(6); %Number of relay pairs Type 4
result(5,1)= Case0(7)+Case0(8)+Case0(9); %Number of relay pairs Type 5
result(6,1)= Case0(10)+Case0(11)+Case0(12)+Case0(13)+Case0(14)+Case0(15); %Number of relay pairs Type 6
Nf=result(1,1)+result(2,1)+result(3,1)+result(4,1)+result(5,1);% Number of calculable sep time backup-main relay pair
Nnf=result(6,1);% Number of Non-Feasible relay pairs
N=Nf+Nnf;% Total pairs
Nnosen=Case0(3)+Case0(7)+Case0(12)+Case0(15);
Nnosel=ki;% Pairs with loss of selectivity
Nnosel2=kki;% Pairs with loss of selectivity
%% Performance indexes
% selectivity
sel=(1-Nnosel/Nf)*100;%selectivity level index
sel2=(1-Nnosel2/Nf)*100;%selectivity level index
minSepTime=(min(SepTime(:,1)));%Minimum separation time (seconds)
maxSepTime=(max(SepTime(:,1)));%Minimum separation time (seconds)
meanSepTime=mean(SepTime(:,1));%Mean Separation Time (seconds)
devSepTime=std(SepTime(:,1));%Std. Dev. Separation Time (seconds)
% sensitivity
sen=100*(1-Nnosen/(N));%sensitivity level index
% speed
T(T==0) = [];
size(T)
meanPrimTime=mean(T);% Average primary operation time (seconds)
%AvgPrimSpeed=1/meanPrimTime; %Average primary speed (1/seconds)
AvgPrimSpeed=mean(T.^-1); %Average primary speed (1/seconds)
devPrimTime=std(T);% Average primary operation time (seconds)
meanBackTime=meanPrimTime+meanSepTime;% Average backup operation time (seconds)
% Screen output
 
disp('*******************************************************')
%fprintf('Case study:%s\n',casestudy)
fprintf('Simulation results:\n')
fprintf('Prefault power flow included?: %s\n',reply2)
fprintf('Number of fauls per line %d\n',nf)
fprintf('Simulated primary-backup pairs %d\n',(N))
fprintf('Simulated primary relays %d\n',(N))
fprintf('Simulated backup pairs %d\n',(N))
fprintf('Relay polarization angle -45 deg  < angle(VI^*)< +135 deg\n')
fprintf('___________________________________________________________________________________\n');
fprintf('Relay response classification:\n');
fprintf('Type 1 Normal operation pairs p-q  :  %4d  %4.1f %%\n',result(1,1), 100*result(1,1)/(N));
fprintf('Type 2 Normal operation pairs j-i  :  %4d  %4.1f %%\n',result(2,1), 100*result(2,1)/(N));
fprintf('Type 3 Partial operation relay  j  :  %4d  %4.1f %%\n',result(3,1), 100*result(3,1)/(N));
fprintf('Type 4 Partial operation relay  i  :  %4d  %4.1f %%\n',result(4,1), 100*result(4,1)/(N));
fprintf('Type 5 Partial operation relays j-i:  %4d  %4.1f %%\n',result(5,1), 100*result(5,1)/(N));
fprintf('Type 6 No operation                :  %4d  %4.1f %%\n',result(6,1), 100*result(6,1)/(N));
fprintf('___________________________________________________________________________________\n');
fprintf('Back-main pairs with calculated separation time                  :  %4d,  %4.1f %%\n',Nf, 100*Nf/(N));
fprintf('Back-main pairs where no separation time can be calculated       :  %4d,  %4.1f %%\n',Nnf, 100*Nnf/(N));
fprintf('Backup relays with no loss of sensitivity                        :  %4d,  %4.1f %%\n',N-Nnosen, sen);
fprintf('Mean   Operation Times with specified TDSs and Ips               : %7.4f (ms)\n',meanPrimTime*1000)
fprintf('StdDev Operation Times with specified TDSs and IPs               : %7.4f (ms)\n',devPrimTime*1000)
fprintf('Mean   Separation Time with specified TDSs and Ips               : %7.4f (ms)\n',meanSepTime*1000)
fprintf('StdDev Separation with specified TDSs and IPs               : %7.4f (ms)\n',devSepTime*1000)
fprintf('Result: %5.4f %% of the faults accomplishes the allowable coordination interval C=%7.4f\n',sel,Co)
fprintf('___________________________________________________________________________________\n');
fprintf('Backup relays with loss of sensitivity                            :  %4d      \n',Nnosen);
fprintf('Relay pairs with loss of selectivity (Original method, Dt > Sm)   :  %4d      \n',Nnosel);
fprintf('Relay pairs with loss of selectivity (Dt > Sm & Dt > k * Tprimary):  %4d      \n',Nnosel2);
fprintf('Relay pairs with feasible interval time calculation               :  %4d      \n',Nf);
fprintf('Sensitivity Index                                                 :  %4.2f %% \n',sen);
fprintf('Selectivity Index 1 (Original method, Dt > Sm)                    :  %4.2f %% \n',sel);
fprintf('Selectivity Index 2 (Dt > Sm & Dt > kt Tprimary)                  :  %4.2f %% \n',sel2);
fprintf('Average Speed Index                                               :  %4.3f 1/s\n',AvgPrimSpeed);
fprintf('___________________________________________________________________________________\n');
fprintf('Elapsed simulation time: %6.2f s \n',elapsedtime000)
disp('****************************************************************************************')
close all
% Figure 1 - Histogram of Separation Times
figure('name','Separation Times Histogram (s)','position',[0, 300, 400, 200])
xbins = -0.50:0.01:2;
h=histogram(SepTime(:,1),xbins,'FaceColor','yellow');
set(gcf,'color','w')
counts = h.Values';
% Figure 2 - Histogram of System Primary Relay Operation Times
figure('name','SPROT Histogram (s)','position',[0, 0, 400, 200])
% SROT Histogram
h=histogram((T.^-1),200,'FaceColor','yellow');
set(gcf,'color','w')
counts = h.Values;
% %Save results for latex table
% res(1,1)=N; % Total pairs
% res(2,1)=Nf; % Total feasible pairs
% res(3,1)=Nnf; % Total non-feasible pairs
% res(4,1)=sen; % percentage
% res(5,1)=sel; % percentage
% res(6,1)=AvgPrimSpeed;%1/s
% res(7,1)=meanPrimTime*1000; %ms
% res(8,1)=devPrimTime*1000;%ms
% res(9,1)=meanSepTime*1000;%ms
% res(10,1)=minSepTime*1000;%ms
% res(11,1)=maxSepTime*1;%s
% res(12,1)=elapsedtime000;% CPU time
% Case0=Case0';
% ll=length(T);
% if ll < nr*nf
% for kiter=ll+1:nr*nf
%     T(kiter)=meanPrimTime;
% end
% end
% display curves according to De Oliveira-De Jesus, Paulo M., and Elmer Sorrentino. "A graphical tool to examine the coordination details of directional overcurrent protections as a function of fault location." IEEE Transactions on Power Delivery 40.1 (2024): 659-662.
DisplayShortCircuitCurves
%DisplayShortCircuitCurves2
DisplaySelectivityCurves
%leg = legend({'Relay-pair'}, 'FontSize', 12, 'Location','best');
%set(leg,'Interpreter','latex');
% size(Bcoordc)
% sum(sum(Bcoordc))
% sum(SepTime(:,1))/length(SepTime(:,1))
% Aresp=[AvgPrimSpeed meanPrimTime*1000 sel sel2 Nf];
% dictiorelays=[1 2 3 4 5 6 7 8 9 10 11 12 13 14; 
%               2 9 3 10 4 11 5 12 6 13 1 8 14 7];  % original case numbering 
%   for kk=1:length(Dsol)  
%   TAPs(dictiorelays(2,kk),1)= Dsol(kk); 
%   end
% Aresp2=[Rf; iterxx; FVAL; Nf; 100*Nf/(N); AvgPrimSpeed; meanPrimTime*1000; devPrimTime*1000; meanSepTime*1000; devSepTime*1000;minSepTime*1000; maxSepTime*1000; sel; sel2; TAPs;length(Bnew);Case0];