% Print short circuit results
% according to De Oliveira-De Jesus, Paulo M., and Elmer Sorrentino. "A graphical tool to examine the coordination details of directional overcurrent protections as a function of fault location." IEEE Transactions on Power Delivery 40.1 (2024): 659-662.



global Dmin Co nr bdat ldat K Ip Sbase Vbase Zbase Ibase econv itermax tdat nlf D qmax lowerbound upperbound npr reversest linenumber back main dictiolines dictiorelays
figure ('Color','w','units','normalized','outerposition',[0 0 1 1],'name','Coordination Results','numbertitle','off')
for jj=1:npr
    clear xx yy zz ww
iterxx=0;
for kj=1:length(ShortC(:,1))
if ShortC(kj,4)==main(jj) & ShortC(kj,5)==back(jj) & ShortC(kj,3)==linenumber(jj) 
 iterxx=iterxx+1;
xx(iterxx,1)=ShortC(kj,2);%distance
yy(iterxx,1)=ShortC(kj,7);%Ip,back
ww(iterxx,1)=ShortC(kj,9);%Ij,back
zz(iterxx,1)=ShortC(kj,6);%Realy pair Case
uu(iterxx,1)=ShortC(kj,12);%Ip,pick
ww2(iterxx,1)=ShortC(kj,11);%Ij',back
uu2(iterxx,1)=ShortC(kj,14);%Ij,pick
line=ShortC(kj,3);
end
end
B=[xx,yy,zz,ww,uu,ww2,uu2];
Bx=sortrows(B,1,'descend');
subplot(5,4,jj)
fsize=10;
yyaxis left
if any(reversest==jj)
h1=plot(B(:,1),Bx(:,2),'-','color','r','LineWidth',1);%
% hold on
% plot(B(:,1),Bx(:,5),'g.','LineWidth',2) 
hold on
h2=plot(B(:,1),Bx(:,4),'-','color','r','LineWidth',1);%
% hold on
% plot(B(:,1),Bx(:,7),'g.','LineWidth',2) 
hold on
h3=plot(B(:,1),Bx(:,6),'b-','LineWidth',2);% 
hold on
h3x=plot(B(:,1),zeros(1000,1),'k-','LineWidth',1);% 
else
h1=plot(B(:,1),B(:,2),'color','r','LineWidth',1);%
% hold on
% plot(B(:,1),B(:,5),'g.','LineWidth',2) %
hold on
h2=plot(B(:,1),B(:,4),'-','color','r','LineWidth',1);%
% hold on
% plot(B(:,1),B(:,7),'g.','LineWidth',2) 
hold on
h3=plot(B(:,1),B(:,6),'b-','LineWidth',2);% 
hold on
h3x=plot(B(:,1),zeros(1000,1),'k-','LineWidth',1);% 
end
set(gca,'FontSize',fsize)
xticks([.1 .2 .3 .4 .5 .6 .7 .8 .9 1.0 ])
ax = gca;
ax.YColor = 'k';

if  max(B(:,4)) >0
%ylim([0 max(B(:,6))]) 
else
    ylim([0 6])
end
%ylim([0 5])
%yticks([.1 .2 .3 .4 .5 .6 .7 .8  ])
xaxis = xlabel({['Fault location distance (d) in pu from main DOCP '  num2str(dictiorelays(2,main(jj)))]}, 'FontSize', fsize,'Color','k');
%set(xaxis,'Interpreter','latex');
yaxis = ylabel({'Fault current (kA)'}, 'FontSize', fsize);
%set(yaxis,'Interpreter','latex');
hold on
yyaxis right
set(gca,'FontSize',fsize)
if any(reversest==jj)
h4=plot(B(:,1),Bx(:,3),'k-','LineWidth',2); %
else
h4=plot(B(:,1),B(:,3),'k-','LineWidth',2); %
end
yyaxis = ylabel({'DOCP-pair case (1-15)'}, 'FontSize', fsize);
%set(yyaxis,'Interpreter','latex','Color','r');
ymax=15;
 ylim([1 ymax])
 yticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])

 
 txt = ['Line: ' num2str(dictiolines(2,line)) ' DOCP Pair: main: '  num2str(dictiorelays(2,main(jj))) ' - backup: ' num2str(dictiorelays(2,back(jj)))];
text(0.05,ymax*1.05,txt)

ax = gca;
ax.YColor = 'k';
ax.XColor = 'k';
%legend([h1 h2 h3 h4], 'First line', 'Second line', 'Third line', '4th
%line');`
%legend([h2 h3 h4], ['Backup Isc',num2str(dictiorelays(2,back(jj)))],['Backup Isc`',num2str(dictiorelays(2,back(jj)))],'Case number', 'FontSize', 8, 'Location','best');
end



% % Plot some data
% figure;
% h1 = plot(rand(1,7), rand(1,7));
% hold on; 
% h2 = plot(rand(1,7), rand(1,7));
% h3 = plot(rand(1,7), rand(1,7));
% % Label each line in the legend
% legend([h1 h2 h3], 'First line', 'Second line', 'Third line');


