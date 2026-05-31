global Dmin  Co nr bdat ldat K Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf D qmax lowerbound upperbound npr reversest linenumber back main dictiolines dictiorelays
 
figure ('Color','w','units','normalized','outerposition',[0 0 1 1],'name','Coordination Results','numbertitle','off')
%for jj=18:18 
    for jj=1:npr 
    clear xx yy zz ww
iterxx=0;
for kj=1:length(SepTime(:,1))
if SepTime(kj,4)==main(jj) & SepTime(kj,5)==back(jj) & SepTime(kj,3)==linenumber(jj) 
 iterxx=iterxx+1;
xx(iterxx,1)=SepTime(kj,2);%localization
% if SepTime(kj,1) < .8
yy(iterxx,1)=SepTime(kj,1);%septime
if SepTime(kj,7) > Co/(ktimes)
yy2(iterxx,1)=(ktimes)*SepTime(kj,7);%septime
else
yy2(iterxx,1)=Co;%septime
end

ww(iterxx,1)=SepTime(kj,3);%

if SepTime(kj,6)==1
zz(iterxx,1)=1;%pair type
elseif SepTime(kj,6)==2
zz(iterxx,1)=2;    
elseif SepTime(kj,6)==3
zz(iterxx,1)=3;    
elseif SepTime(kj,6)==4
zz(iterxx,1)=3;   
elseif SepTime(kj,6)==5
zz(iterxx,1)=3; 
else
zz(iterxx,1)=6;  
end
uu=Co*ones(length(ww),1);
line=SepTime(kj,3);
end
end
if nf-length(xx)==0 
B=[xx,yy,zz,ww,uu,yy2];% loc, septime, type, line
else
    yy1=zeros(nf,1);
    zz1=ones(nf,1)*6;
    ww1=zeros(nf,1);
    uu1=ones(nf,1)*0;
    uu2=ones(nf,1)*0;
for k=1:nf
    xx1(k,1)=k/(nf+1);
end 
for k=1:nf
   for j=1:length(xx)
   if abs(xx(j,1)-xx1(k,1))< 0.00001
       yy1(k,1)=yy(j,1);
       zz1(k,1)=zz(j,1);
       ww1(k,1)=ww(j,1);
       uu1(k,1)=uu(j,1);
       uu2(k,1)=yy2(j,1);
   end  
   end
end
B=[xx1,yy1,zz1,ww1,uu1,uu2];
end 
Bx=sortrows(B,1,'descend');
subplot(5,4,jj)
fsize=12;
yyaxis left
if any(reversest==jj)
plot(B(:,1),Bx(:,2),'color', [0 .5 0],'LineWidth',2) %
hold on
plot(B(:,1),B(:,5),'color', [0.5,0.5,0.5],'LineWidth',2)  
hold on
plot(B(:,1),Bx(:,6),'color',  [0.5,0.0,0.5],'LineWidth',2) %
else
plot(B(:,1),B(:,2),'color', [0 .5 0],'LineWidth',2) %
 hold on
 plot(B(:,1),B(:,5),'color',  [0.5,0.5,0.5],'LineWidth',2) %
hold on
plot(B(:,1),B(:,6),'color',  [0.5,0.0,0.5],'LineWidth',2) %
end
set(gca,'FontSize',fsize)
 
tope=max([max(B(:,2));max(Bx(:,2));max(B(:,6));max(Bx(:,6))]); 
if tope  > 10  
set(gca,'YScale', 'log')
yticks([1 10 100])
yaxis = ylabel({'{\Delta}t (s, log scale)'}, 'FontSize', fsize,'color',  [0 0 0]);
else
ylim([0.05 tope])
yaxis = ylabel({'{\Delta}t (s)'}, 'FontSize', fsize,'color',  [0 0 0]);
end

xaxis = xlabel({['Fault location distance (d) in pu from main DOCP '  num2str(dictiorelays(2,main(jj)))]}, 'FontSize', fsize);


set(yaxis,'color',  [0 0 0]);
ax = gca;
ax.YColor = '[0 0 0]';
hold on
yyaxis right
if any(reversest==jj)
plot(B(:,1),Bx(:,3),'k-','LineWidth',2) %
else
plot(B(:,1),B(:,3),'k-','LineWidth',2) %
end
yyaxis = ylabel({'DOCP pair Type'}, 'FontSize', fsize,'Color','k');
set(yyaxis,'Color','k');
ax = gca;
ax.YColor = 'k';
ylim([1 10])
yticks([1 2 3 4 5 6])
%yticklabels({'1','2','3','4','5','6'})

ticklabels2 = {'1','2','3','4','5','6',' ',' ',' ',' '};
% prepend a color for each tick label
ticklabels_new = cell(size(ticklabels2));
for i = 1:length(ticklabels2)
    ticklabels_new{i} = ['\color{black} ' ticklabels2{i}];
end
% set the tick labels
set(gca, 'YTickLabel', ticklabels_new);



 

txt = ['Line: ' num2str(dictiolines(2,line)) ' DOCP Pair: main: '  num2str(dictiorelays(2,main(jj))) ' - backup: ' num2str(dictiorelays(2,back(jj)))];
text(0.05,10+0.5,txt)





a=[B(:,1),B(:,2),B(:,4),B(:,3)];
a(a(:,2)==0,:) = [];
a(a(:,1)>=0.98,:) = [];
a(a(:,1)<=0.02,:) = [];
% a(a(:,1)>=upperbound*0.99,:) = [];
% a(a(:,1)<=lowerbound*1.01,:) = [];
idx(jj)=find(a(:,2)==min(a(:,2)));
Localiz(jj,1)=a(idx(jj),1);% localiz
Localiz(jj,2)=a(idx(jj),2);% selectivity
Localiz(jj,3)=a(idx(jj),3);% line
Localiz(jj,4)=a(idx(jj),4);% type
Localiz(jj,5)=jj;
end