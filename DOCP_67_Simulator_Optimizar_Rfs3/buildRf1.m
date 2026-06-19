%% Build the relay coordination model (B matrix) from the simulator logic
global Rf
Case0=zeros(15,1);% Initialize type pairs vector
%reply2='n'; %no PREFAULT LOADFLOW
ncase=0; 
halfSm2x=[];Bc=[];B_case=[];iter=1;%initialize counter for separation times
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

Case0=Case0+Case;%All 15 pair types classified are aggregated here
halfSm2x=[halfSm2;halfSm2x];
Bc=[Bc;Bcoord];%Coordination matrix
B_case=[B_case;Bcoord_case];%Coordination matrix with case identification 
% T=unique(vertcat(unique([Tix';Tq']),T));%All nr primary times are aggregated here
% T(T==0)=[];
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
end 


halfSm=halfSm2x*ktimes;
sum(halfSm2x);
%sum(Tppal)
 % size(halfSmx0)
 % pause
%Bextended=[Bc,B_case];% coordination matrix,last column indicates the case number
