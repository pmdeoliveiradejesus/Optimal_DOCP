function [halfSm2,S,Sx,Mx,fl,nlf,Co,Tix,Tq,index,nr,Iq,Ipback,Ipq,Ipp,Ii,Ij,Iip,Ijp,Ipi,Ipj,CaseType,distance,Bcoord,Bcoord_case] =  run_classification4(x,ncase,reply2,nlf,D,Ip,K)
global Dmin Co nr reply2 bdat ldat K Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf qmax lowerbound upperbound npr reversest linenumber back main dictiolines dictiorelays
%relay polarization angle (deg) -45 < angle(V.conj(I))<+135 deg -to detect opposite currents
tmax=45+360;%fault angle current tripping zone
tmin=225;%fault angle current tripping zone
fl=zeros(15,1);flag=0;flag2=0;
for jj=1:nlf
L(1,jj)=x;% fault distance xmin < x < xmax
end
[index]=run_shortcircuit(L(1,:),reply2);%Run external short-circuit program
npr=length(index(:,1));
%%all angles > 0
for k=1:length(index(:,1))
    index(k,17)= index(k,17)+360;
    index(k,18)= index(k,18)+360;
    index(k,19)= index(k,19)+360;
    index(k,20)= index(k,20)+360;
end 
for k=1:length(index(:,1)) %number of relay pairs
beta(index(k,3),index(k,4))= index(k,13);  %beta i main first interval due to fault in line k
beta(index(k,1),index(k,2))= index(k,14);  %beta j backup first interval due to fault in line k
betap(index(k,3),index(k,4))= index(k,15); %beta i' main transient due to fault in line k when q is open
betap(index(k,1),index(k,2))= index(k,16); %beta j' backup transient due to fault in line k when q is open
 theta(index(k,3),index(k,4))= index(k,17); %theta i main first interval due to fault in line k
 theta(index(k,1),index(k,2))= index(k,18); %theta j backup first interval due to fault in line k
 thetap(index(k,3),index(k,4))= index(k,19);%theta i' main transient due to fault in line k when q is open
 thetap(index(k,1),index(k,2))= index(k,20);%theta j' backup transient due to fault in line k when q is open
gammap(index(k,22),index(k,1),index(k,3),index(k,2))= index(k,21);% gamma' q j i k
gammapp(index(k,22),index(k,1),index(k,3),index(k,2))= index(k,24);% gamma'' q j i k  relay j
gammappp(index(k,22),index(k,1),index(k,3),index(k,2))= index(k,25);% gamma''' q j i k relaj i
end
i=zeros(1,length(index(:,1)));
q=zeros(1,length(index(:,1)));
p=zeros(2,length(index(:,1)));
j=zeros(2,length(index(:,1)));
Tix=zeros(1,length(index(:,1)));
Tq=zeros(1,length(index(:,1)));
 
for k=1:length(index(:,1))
 if   beta(index(k,22),index(k,4))*D(index(k,22)) > 0
 if  beta(index(k,3),index(k,4))*D(index(k,3)) > beta(index(k,22),index(k,4))*D(index(k,22)) %Who trips first
 q(k)=index(k,22);
 i(k)=index(k,3);
 J=find(index(:,22)==q(k));
 for ii=1:length(J)
     j(ii,k)=index(J(ii),1);
 end
  P=find(index(:,3)==q(k));
 for ii=1:length(P)
     p(ii,k)=index(P(ii),1);
 end
 end
 else
 q(k)=index(k,3);
 i(k)=index(k,22);
 J=find(index(:,3)==q(k));
 for ii=1:length(J)
     j(ii,k)=index(J(ii),1);
 end
  P=find(index(:,22)==q(k));
 for ii=1:length(P)
     p(ii,k)=index(P(ii),1);
 end
end
end % Identifies for every fault who is q, i, j and p
Ad=vertcat(i,j);
Bd=unique(Ad','rows')';
Bd(:,1) = [];
Ai=vertcat(q,p);
Bi=unique(Ai','rows')';
Bi(:,1) = [];
iter=1;
for h=1:nlf
for m=2:3
if  Bd(m,h) > 0
C1(iter,1)=Bd(m,h);% This is j
C1(iter,2)=Bd(1,h);% This is i
C1(iter,3)=h;%This is faulted line k
C1(iter,4)=Bi(1,h); %This is q
C1(iter,5)=Bi(m,h); %This is p
if Bi(m,h) == 0
C1(iter,5)=Bi(m-1,h); %This is p
end
iter=iter+1;
end
end
end%C1 indicates tripping sequence ij and identifies j i q p
iter=1;
for h=1:nlf
for m=2:3
if  Bi(m,h) > 0
C2(iter,1)=Bi(m,h);%This is p
C2(iter,2)=Bi(1,h);%This is q
C2(iter,3)=h;%This is faulted line k
C2(iter,4)=Bd(m,h);% This is j
C2(iter,5)=Bd(1,h);% This is i
iter=iter+1;
end
end
end%C2 indicates reverse tripping sequence
Ipp=zeros(npr,1);
Ipj=zeros(npr,1);
Ipi=zeros(npr,1);
Ipq=zeros(npr,1);
Ipback=zeros(npr,1);
Iq=zeros(npr,1);
Ii=zeros(npr,1);
Ij=zeros(npr,1);
Iip=zeros(npr,1);
Ijp=zeros(npr,1);
% Separation times calculation
%%---------------------------------------------------------
% Type 1  - Case 1 - Normal Operation qp
for g=1:length(C2(:,1))
if theta(C2(g,1),C2(g,3)) < tmax &...%period 1: no reverse current relay p
   theta(C2(g,1),C2(g,3)) > tmin &... %period 1: no reverse current relay p
        beta(C2(g,1),C2(g,3)) > 0
%period 1: no loss of sensitivity relay p
flag=flag+1;
flag2=flag2+1;
fl(1)=fl(1)+1;
S(flag)=beta(C2(g,1),C2(g,3))*D(C2(g,1))-beta(C2(g,2),C2(g,3))*D(C2(g,2));
C2(g,1);
C2(g,2);
C2(g,3);
flag; 
Bcoord(flag,C2(g,1))=beta(C2(g,1),C2(g,3));
Bcoord(flag,C2(g,2))=-beta(C2(g,2),C2(g,3));
Bcoord_case(flag,1)=1;
Sx(flag,1)=x;
Sx(flag,2)=C2(g,3);%line
Sx(flag,3)=C2(g,2);%relay q
Sx(flag,4)=C2(g,1);%relay p
Sx(flag,5)=1;
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
Tp=beta(C2(g,1),C2(g,3))*D(C2(g,1));
Sx(flag,6)=Tq(flag);
halfSm2(flag,1)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
%disp('back p')
C2(g,1);
%disp('main q')
C2(g,2);
distance(flag2)=x;
Iq(flag2)=Ip(C2(g,2))*(K(1,C2(g,2))/beta(C2(g,2),C2(g,3))+1)^(1/K(2,C2(g,2)));
Ipback(flag2)=Ip(C2(g,1))*(K(1,C2(g,1))/beta(C2(g,1),C2(g,3))+1)^(1/K(2,C2(g,1)));
CaseType(flag2)=1;
Mx(flag2,1)=x;%line
Mx(flag2,2)=C2(g,3);%line
Mx(flag2,3)=C2(g,2);%relay q
Mx(flag2,4)=C2(g,1);%relay p
Mx(flag2,5)=1;%case
Ii(flag2)=0;
Ij(flag2)=0;
Iip(flag2)=0;
Ijp(flag2)=0;
end
end% End Type 1  - Case 1 - Normal Operation qp:
%Relays q and p are sensitive: $\beta_{qkh}>0$,$\beta_{pkh}>0$. There is no reverse current at relay p: $|\theta'_{pkh}-\theta'_{qkh}|<\phi$

%%---------------------------------------------------------
% Type 2  - Case 2   - Normal Operation ij
for g=1:length(C1(:,1))
if theta(C1(g,1),C1(g,3)) < tmax &...%period 1: no reverse current relay j
   theta(C1(g,1),C1(g,3)) > tmin &... %period 1: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0  %period 2: no loss of sensitivity relay i
flag=flag+1;
flag2=flag2+1;
fl(2)=fl(2)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammap(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord(flag,C1(g,4))=-gammap(C1(g,4),C1(g,1),C1(g,2),C1(g,3));
Bcoord_case(flag,1)=2;

% disp('Type 2')
% betap(C1(g,1),C1(g,3))
% D(C1(g,1))
% -betap(C1(g,2),C1(g,3))
% D(C1(g,2))
% -gammap(C1(g,4),C1(g,1),C1(g,2),C1(g,3))
% D(C1(g,4))
% pause

Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);%line 
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=2;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tqq+Tpj*(1-Tqq/Tj);
Tjx(flag)-Tix(flag);
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tqq+Tpi*(1-Tqq/Ti);
%disp('back j')
C1(g,1);
%disp('slower main i')
C1(g,2);
%disp('faster main q')
C1(g,4);
distance(flag2)=x;
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
% Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=2;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);%line 
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=2;
end
end% End of Type 2  - Case 2
%Normal Operation ij:
%Relays i and j are sensitive in both periods: $\beta_{ikh}>0$, $\beta_{jkh}>0$, $\beta'_{ikh}>0$, $\beta'_{jkh}>0$.
%There is  no reverse current at relay j in both periods: $|\theta_{jkh}-\theta_{iqkh}|<\phi$ and $|\theta'_{jkh}-\theta'_{ikh}|<\phi$. A separation time can be calculated.


%%---------------------------------------------------------
% Type 3a - Case 3
for g=1:length(C1(:,1))
if theta(C1(g,1),C1(g,3)) < tmax &...%period 1: no reverse current relay j
   theta(C1(g,1),C1(g,3)) > tmin &... %period 1: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0%period 2: no loss of sensitivity relay i
flag=flag+1;
flag2=flag2+1;
fl(3)=fl(3)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord(flag,C1(g,4))=-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3));
Bcoord_case(flag,1)=3;



Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=3;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tix(flag)+S(flag);
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tqq+Tpi*(1-Tqq/Ti);
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
Ipj(flag2)=Ip(C1(g,1));
CaseType(flag)=3;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=3;
end
end% End  Type 3a - Case 3
% Partial Operation ij:
% Relay j is not sensitive at period 1: $\beta_{jkh}<0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i:is sensitive at both periods $\beta_{ikh}>0$, $\beta'_{ikh}>0$.
% There is  no reverse current at relay j in both periods: $|\theta_{jkh}-\theta_{iqkh}|<\phi$ and $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.

%%---------------------------------------------------------
% Type 3b - Case 4
for g=1:length(C1(:,1))
if or(theta(C1(g,1),C1(g,3)) > tmax,theta(C1(g,1),C1(g,3)) < tmin) &... %period 1: yes reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0 %period 2: no loss of sensitivity relay i
   flag=flag+1;
      flag2=flag2+1;
fl(4)=fl(4)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
for k=1:nr
Bcoord(flag,k)=0;
end
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord(flag,C1(g,4))=-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3));
Bcoord_case(flag,1)=4;

% disp('CASE 4')
% betap(C1(g,1),C1(g,3))
% D(C1(g,1))
% -betap(C1(g,2),C1(g,3))
% D(C1(g,2))
% -gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))
% D(C1(g,4))
% pause
Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=4;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tix(flag)+S(flag);
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tqq+Tpi*(1-Tqq/Ti);
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
%Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=4;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=4;
end
end% End Type 3b - Case 4
% Partial Operation ij:
% Relay j sees a reverse current at period 1: $|\theta_{jkh}-\theta_{ikh}|>phi$
% Relay j is sensitive at period 1: $\beta_{jkh}>0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i: is sensitive at both periods $\beta_{ikh}>0$, $\beta'_{ikh}>0$.
% There is  no reverse current at relay j in period 2: $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.

%%---------------------------------------------------------
% Type 3c - Case 5
for g=1:length(C1(:,1))
if or(theta(C1(g,1),C1(g,3)) > tmax,theta(C1(g,1),C1(g,3)) < tmin) &... %period 1: yes reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0 %period 2: no loss of sensitivity relay iflag=flag+1;
fl(5)=fl(5)+1;
   flag=flag+1;
      flag2=flag2+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord(flag,C1(g,4))=-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3));
Bcoord_case(flag,1)=5;

% disp('CASE 5')
% betap(C1(g,1),C1(g,3))
% D(C1(g,1))
% -betap(C1(g,2),C1(g,3))
% D(C1(g,2))
% -gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))
% D(C1(g,4))
% pause
Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=5;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tix(flag)+S(flag);
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tqq+Tpi*(1-Tqq/Ti);
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;
% disp('Iq')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
%Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=5;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=5;
end
end
% End 3c - Case 5
% Partial Operation ij:
% Relay j sees a reverse current at period 1: $|\theta_{jkh}-\theta_{iqkh}|>\phi$
% Relay j is no sensitive at period 1: $\beta_{jkh}<0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i: is sensitive at both periods $\beta_{ikh}>0$, $\beta'_{ikh}>0$.
% There is  no reverse current at relay j in period 2: $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.
    %%---------------------------------------------------------
    % Type 4  - Case 6
    for g=1:length(C1(:,1))
if theta(C1(g,1),C1(g,3)) < tmax &...%period 1: no reverse current relay j
   theta(C1(g,1),C1(g,3)) > tmin &... %period 1: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) >0  %period 2: no loss of sensitivity relay iif abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) < qmax &  abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax & beta(C1(g,1),C1(g,3)) > 0 & betap(C1(g,1),C1(g,3)) > 0 & beta(C1(g,2),C1(g,3)) < 0 & betap(C1(g,2),C1(g,3)) > 0  & beta(C1(g,4),C1(g,3))*D(C1(g,4))
flag=flag+1;
flag2=flag2+1;
fl(6)=fl(6)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord(flag,C1(g,4))=-gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3));
Bcoord_case(flag,1)=6;

Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=6;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tjx(flag)=Tqq+Tpj*(1-Tqq/Tj);
Tix(flag)=Tjx(flag)-S(flag);
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tqq+Tpj*(1-Tqq/Tj)-(betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4)));
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=6;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=6;
end
end
% End  Type 4  - Case 6
% Partial Operation ij:
% Relay j is sensitive at period 1: $\beta_{jkh}>0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i is not sensitive at period 1: $\beta_{ikh}<0$
% Relay i is sensitive at period 2: $\beta'_{ikh}>0$
% There is  no reverse current at relay j in both periods: $|\theta_{jkh}-\theta_{iqkh}|<\phi$ and $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.

%%---------------------------------------------------------
% Type 5a - Case 7
for g=1:length(C1(:,1))
if theta(C1(g,1),C1(g,3)) < tmax &...%period 1: no reverse current relay j
   theta(C1(g,1),C1(g,3)) > tmin &... %period 1: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0  %period 2: no loss of sensitivity relay i
fl(7)=fl(7)+1;
flag=flag+1;
flag2=flag2+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-0*gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord_case(flag,1)=7;

Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=7;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tpi;
Tjx(flag)=Tpj;
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tpi;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=7;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=7;
end
end
% End  Type 5a - Case 7
% Partial Operation ij:
% Relay j is not sensitive at period 1: $\beta_{jkh}<0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i is not sensitive at period 1: $\beta_{ikh}<0$
% Relay i is sensitive at period 2: $\beta'_{ikh}>0$
% There is  no reverse current at relay j period 1: $|\theta_{jkh}-\theta_{iqkh}|<\phi$.
% There is  no reverse current at relay j period 2: $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.

%%---------------------------------------------------------
% Type 5b - Case 8
for g=1:length(C1(:,1))
if or(theta(C1(g,1),C1(g,3)) > tmax ,theta(C1(g,1),C1(g,3))< tmin) &... %period 1: yes reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) >0  %period 2: no loss of sensitivity relay iif abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &  abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax & beta(C1(g,1),C1(g,3)) > 0 & betap(C1(g,1),C1(g,3)) > 0 & beta(C1(g,2),C1(g,3)) < 0 & betap(C1(g,2),C1(g,3)) > 0  & beta(C1(g,4),C1(g,3))*D(C1(g,4))
flag=flag+1;
flag2=flag2+1;
fl(8)=fl(8)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-0*gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord_case(flag,1)=8;
Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=8;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tpi;
Tjx(flag)=Tpj;
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tpi;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
% disp('Iq')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
%Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=8;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=8;
end
end
% End  Type 5b - Case 8
% Partial Operation ij:
% Relay j is sensitive at period 1: $\beta_{jkh}>0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i is not sensitive at period 1: $\beta_{ikh}<0$
% Relay i is sensitive at period 2: $\beta'_{ikh}>0$
% There is reverse current at relay j period 1: $|\theta_{jkh}-\theta_{iqkh}|>\phi$.
% There is  no reverse current at relay j period 2: $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.

%%---------------------------------------------------------
% Type 5c - Case 9
for g=1:length(C1(:,1))
if or(theta(C1(g,1),C1(g,3)) > tmax ,theta(C1(g,1),C1(g,3))< tmin) &... %period 1: yes reverse current relay j
   thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin &...  %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) >0  %period 2: no loss of sensitivity relay iif abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &  abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax & beta(C1(g,1),C1(g,3)) > 0 & betap(C1(g,1),C1(g,3)) > 0 & beta(C1(g,2),C1(g,3)) < 0 & betap(C1(g,2),C1(g,3)) > 0  & beta(C1(g,4),C1(g,3))*D(C1(g,4))
flag=flag+1;
flag2=flag2+1;
fl(9)=fl(9)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-0*gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
%Bcoord=zeros(flag,nr);
Bcoord(flag,C1(g,1))=betap(C1(g,1),C1(g,3));
Bcoord(flag,C1(g,2))=-betap(C1(g,2),C1(g,3));
Bcoord_case(flag,1)=9;
Sx(flag,1)=x;
Sx(flag,2)=C1(g,3);
Sx(flag,3)=C1(g,2);%relay i
Sx(flag,4)=C1(g,1);%relay j
Sx(flag,5)=9;
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tpi;
Tjx(flag)=Tpj;
Sx(flag,6)=Tix(flag);
halfSm2(flag,1)=Tpi;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
%Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=9;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=9;
end
end
% End  Type 5c - Case 9
% Partial Operation ij:
% Relay j is not sensitive at period 1: $\beta_{jkh}<0$
% Relay j is sensitive at period 2: $\beta'_{jkh}>0$
% Relay i is not sensitive at period 1: $\beta_{ikh}<0$
% Relay i is sensitive at period 2: $\beta'_{ikh}>0$
% There is reverse current at relay j period 1: $|\theta_{jkh}-\theta_{iqkh}|>\phi$.
% There is no reverse current at relay j period 2: $|\theta'_{jkh}-\theta'_{ikh}|<\phi$.
% A separation time can be calculated.

%---------------------------------------------------------
%% Type 6a - Case 10
for g=1:length(C2(:,1))
if or(theta(C2(g,1),C2(g,3)) > tmax ,theta(C2(g,1),C2(g,3))< tmin) &... %period 1: yes reverse current relay j
        beta(C2(g,1),C2(g,3)) > 0  %period 1: no  loss of sensitivity relay p
fl(10)=fl(10)+1;
flag=flag+1;
flag2=flag2+1;
distance(flag2)=x;% disp('Iq')
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
Iq(flag2)=Ip(C2(g,2))*(K(1,C2(g,2))/beta(C2(g,2),C2(g,3))+1)^(1/K(2,C2(g,2)));
Ipq(flag2)=Ip(C2(g,2));
CaseType(flag2)=10;
Mx(flag2,1)=x;%line
Mx(flag2,2)=C2(g,3);%line
Mx(flag2,3)=C2(g,2);%relay q
Mx(flag2,4)=C2(g,1);%relay p
Mx(flag2,5)=10;%case
end
end%
% End  Type 6a - Case 10
% No operation, p does not operate
% There is reverse current at relay p: $|\theta_{pkh}-\theta_{qkh}|>\phi$.
% A separation time can not be calculated.

%%---------------------------------------------------------
% Type 6b - Case 11
for g=1:length(C1(:,1))
    if or(thetap(C1(g,1),C1(g,3)) > tmax ,thetap(C1(g,1),C1(g,3))< tmin) &... %period 2: yes reverse current relay j
betap(C1(g,1),C1(g,3)) > 0 %period 2: no loss of sensitivity relay j
    fl(11)=fl(11)+1;
flag2=flag2+1;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Ipq')
Ipq(flag2)=Ip(C1(g,4));
Ipi(flag2)=Ip(C1(g,2));
Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=11;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=11;
end
end
% End Type 6b - Case 11
% No operation, j does not operate
% There is reverse current at relay j at period 2: $|\theta'_{jkh}-\theta'_{ikh}|>\phi$.
% A separation time can not be calculated.

% Type 6c - Case 12
for g=1:length(C2(:,1))
if theta(C2(g,1),C2(g,3)) < tmax &...%period 1: no reverse current relay p
   theta(C2(g,1),C2(g,3)) > tmin &... %period 1: no reverse current relay p
        beta(C2(g,1),C2(g,3)) < 0  %period 1: yes  loss of sensitivity relay p

%         tmin
%     theta(C2(g,1),C2(g,3))
%     tmax
%     beta(C2(g,1),C2(g,3))
% pause
% 
    
    fl(12)=fl(12)+1;
flag=flag+1;
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
flag2=flag2+1;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
Iq(flag2)=Ip(C2(g,2))*(K(1,C2(g,2))/beta(C2(g,2),C2(g,3))+1)^(1/K(2,C2(g,2)));
Ipback(flag2)=Ip(C2(g,1))*(K(1,C2(g,1))/beta(C2(g,1),C2(g,3))+1)^(1/K(2,C2(g,1)));
Ipq(flag2)=Ip(C2(g,2));
Ipp(flag2)=Ip(C2(g,1));
CaseType(flag2)=12;
Mx(flag2,1)=x;%line
Mx(flag2,2)=C2(g,3);%line
Mx(flag2,3)=C2(g,2);%relay q
Mx(flag2,4)=C2(g,1);%relay p
Mx(flag2,5)=12;%case
end
end
% End Type 6c - Case 12
% No operation, p does not operate
% Relay p is no sensitive: $\beta_{pkh}<0$
% There is no reverse current at relay p: $|\theta_{pkh}-\theta_{qkh}|<\phi$.
% A separation time can not be calculated.



%%---------------------------------------------------------
% Type 6d - Case 13
for g=1:length(C2(:,1))
if or(theta(C2(g,1),C2(g,3)) > tmax ,theta(C2(g,1),C2(g,3))< tmin) &... %period 1: yes reverse current relay j
        beta(C2(g,1),C2(g,3)) < 0  %period 1: yes  loss of sensitivity relay p
fl(13)=fl(13)+1;
flag=flag+1;
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
flag2=flag2+1;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
Iq(flag2)=Ip(C2(g,2))*(K(1,C2(g,2))/beta(C2(g,2),C2(g,3))+1)^(1/K(2,C2(g,2)));
Ipback(flag2)=Ip(C2(g,1))*(K(1,C2(g,1))/beta(C2(g,1),C2(g,3))+1)^(1/K(2,C2(g,1)));
Ipq(flag2)=Ip(C2(g,2));
CaseType(flag2)=13;
Mx(flag2,1)=x;%line
Mx(flag2,2)=C2(g,3);%line
Mx(flag2,3)=C2(g,2);%relay q
Mx(flag2,4)=C2(g,1);%relay p
Mx(flag2,5)=13;%case
end
end
% End of 6d - Case 13
% No operation, p does not operate
% Relay p is no sensitive: $\beta_{pkh}<0$
% There is reverse current at relay p: $|\theta_{pkh}-\theta_{qkh}|>\phi$.
% A separation time can not be calculated.
%%---------------------------------------------------------
% Type 6e - Case 14
for g=1:length(C1(:,1))
if or(thetap(C1(g,1),C1(g,3)) > tmax ,thetap(C1(g,1),C1(g,3))< tmin) &... %period 2: yes reverse current relay j
   betap(C1(g,1),C1(g,3)) < 0 %period 2: yes loss of sensitivity relay j
fl(14)=fl(14)+1;
flag2=flag2+1;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,2))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=14;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=14;
end
end
% End of Type 6e - Case 14
% No operation, j does not operate
% Relay ji is not sensitive at period 2: $\beta'_{ikh}<0$
% A separation time can not be calculated.

%%---------------------------------------------------------
% Type 6f - Case 15
for g=1:length(C1(:,1))
 if  betap(C1(g,1),C1(g,3)) < 0 &...%period 2: yes loss of sensitivity relay j
    thetap(C1(g,1),C1(g,3)) < tmax &... %period 2: no reverse current relay j
   thetap(C1(g,1),C1(g,3)) > tmin   %period 2: no reverse current relay j
 fl(15)=fl(15)+1;
 flag2=flag2+1;
% disp('back')
% C1(g,1)
% disp('main')
% C1(g,2)
distance(flag2)=x;% disp('Iq')
%disp('Ii')
Ii(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/beta(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ij')
Ij(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/beta(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,1)));
%disp('Iip')
Iip(flag2)=Ip(C1(g,2))*(K(1,C1(g,2))/betap(C1(g,2),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ijp')
Ijp(flag2)=Ip(C1(g,1))*(K(1,C1(g,1))/betap(C1(g,1),C1(g,3))+1)^(1/K(2,C1(g,2)));
%disp('Ipq')
% Ipq(flag2)=Ip(C1(g,4));
% Ipi(flag2)=Ip(C1(g,2));
Ipj(flag2)=Ip(C1(g,1));
CaseType(flag2)=15;
Mx(flag2,1)=x;
Mx(flag2,2)=C1(g,3);
Mx(flag2,3)=C1(g,2);%relay i
Mx(flag2,4)=C1(g,1);%relay j
Mx(flag2,5)=15;
 end
end
% End of Type 6f - Case 15

 



 
 