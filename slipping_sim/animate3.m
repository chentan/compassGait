function [sout,wout,bout,dsout,swout] = animate3(name)
scrsz = get(groot,'ScreenSize');
mov = figure('Name','Animation','Position',[1 1 scrsz(3) scrsz(4)]);
pause;

Leg1_color = [0, 0.4470, 0.7410];
Leg2_color = [0.9290, 0.6940, 0.1250];
COM_color = 'k';

%------------------------------------------------------------------%

subplot(3,2,4)
scoord_all = plot(name.tout,name.xout(:,4),'-','color',Leg2_color);
hold on
wcoord_all = plot(name.tout,name.xout(:,5),'-','color',Leg1_color);
plot([name.tout(1) name.tout(end)],[0 0],'k:')
ylabel('Joint vel. (rad/s)')
set(gca,'FontSize',24);
set(scoord_all,'LineWidth',2)
set(wcoord_all,'LineWidth',2)

subplot(3,2,6)
dxs_all = plot(name.tout,name.xout(:,6),'r-');
hold on
plot([name.tout(1) name.tout(end)],[0 0],'k:')
xlabel('Time (s)')
ylabel('Slip vel. (m/s)')
set(gca,'FontSize',24);
set(dxs_all,'LineWidth',2)
dxs_upbound = max(name.xout(:,6));
dxs_lowbound = min(name.xout(:,6));
set(gca,'Ylim',[dxs_lowbound*1.1-1e-9 dxs_upbound*1.1+1e-9]);

%-----------------------------------------------------------%

[pFoot1,pFoot2,pHip,~,pCOM] = limb_pos(name.xout(1,:),name.p_horz(1));

figure(gcf)
subplot(3,2,[1 5])
hold on
cla
Leg1 = line([pFoot1(1);pHip(1)],[pFoot1(2);pHip(2)]);
Leg2 = line([pFoot2(1);pHip(1)],[pFoot2(2);pHip(2)]);
% COM = plot(pCOM(1),pCOM(2),'+');
Leg1COM = plot(0.8*pHip(1)+0.2*pFoot1(1),0.8*pHip(2)+0.2*pFoot1(2),'o');
Leg2COM = plot(0.2*pFoot2(1)+0.8*pHip(1),0.2*pFoot2(2)+0.8*pHip(2),'o');
tdata=text(pCOM(1)+0.4,1.15,[num2str(name.tout(1)),'s']);
ground = line([-10;10],[0;0]);
set(Leg1,'LineWidth',6,'Color',Leg1_color);
set(Leg2,'LineWidth',6,'Color',Leg2_color);
% set(COM,'LineWidth',3,'MarkerSize',8','Color',COM_color);
set(Leg1COM,'MarkerSize',15','Color',Leg1_color,'MarkerFaceColor',Leg1_color);
set(Leg2COM,'MarkerSize',15','Color',Leg2_color,'MarkerFaceColor',Leg2_color);
set(tdata,'FontSize',24);
set(ground,'LineWidth',1.5,'Color',[0.7 0.7 0.7]);
axis square
axis equal
set(gca,'XLim',[pCOM(1)-1 pCOM(1)+1],'Ylim',[-.1 1.2],'FontSize',24);

subplot(3,2,2)
cla
Leg1_close = line([pFoot1(1);pHip(1)],[pFoot1(2);pHip(2)]);
Leg2_close = line([pFoot2(1);pHip(1)],[pFoot2(2);pHip(2)]);
hold on;
COM_close = plot(pCOM(1),pCOM(2),'+');
ground2 = line([-10;10],[0;0]);
set(Leg1_close,'LineWidth',6,'Color',Leg1_color);
set(Leg2_close,'LineWidth',6,'Color',Leg2_color);
set(COM_close,'LineWidth',3,'MarkerSize',8','Color',COM_color);
set(ground2,'LineWidth',1.5,'Color',[0.7 0.7 0.7]);
set(gca,'XLim',[pFoot1(1)-0.01 pFoot1(1)+0.01],'Ylim',[-.01 0.01],'FontSize',24);
xlabel('Closer look of the stance foot')

subplot(3,2,4)
scoord = plot(name.tout(1),name.xout(1,4),'x','color',Leg2_color);
wcoord = plot(name.tout(1),name.xout(1,5),'x','color',Leg1_color);
set(scoord,'MarkerSize',10,'LineWidth',3);
set(wcoord,'MarkerSize',10,'LineWidth',3);
set(gca,'FontSize',24);

subplot(3,2,6)
dxs = plot(name.tout(1),name.xout(1,6),'rx');
set(dxs,'MarkerSize',10,'LineWidth',3);
set(gca,'FontSize',24);

v = VideoWriter('foot_slip_anima.avi');
open(v);

for k=2:20:length(name.xout)
    
    [pFoot1,pFoot2,pHip,~,pCOM]=limb_pos(name.xout(k,:),name.p_horz(k));
    
    subplot(3,2,[1 5]);
    set(Leg1,'XData',[pFoot1(1);pHip(1)],'YData',[pFoot1(2);pHip(2)]);
    set(Leg2,'XData',[pFoot2(1);pHip(1)],'YData',[pFoot2(2);pHip(2)]);
%     set(COM,'XData',pCOM(1),'YData',pCOM(2));
    set(Leg1COM,'XData',0.8*pHip(1)+0.2*pFoot1(1),'YData',0.8*pHip(2)+0.2*pFoot1(2));
    set(Leg2COM,'XData',0.2*pFoot2(1)+0.8*pHip(1),'YData',0.2*pFoot2(2)+0.8*pHip(2));
    set(tdata,'String',[num2str(name.tout(k)),'s'],'Position',[pCOM(1)+0.4 1.15 0]);
    set(tdata,'FontSize',24);
    set(Leg1,'LineWidth',6,'Color',Leg1_color);
    set(Leg2,'LineWidth',6,'Color',Leg2_color);
%     set(COM,'LineWidth',3,'MarkerSize',8','Color',COM_color);   
    set(Leg1COM,'MarkerSize',15','Color',Leg1_color,'MarkerFaceColor',Leg1_color);
    set(Leg2COM,'MarkerSize',15','Color',Leg2_color,'MarkerFaceColor',Leg2_color);
    axis square
    axis equal
    set(gca,'XLim',[pCOM(1)-1 pCOM(1)+1],'Ylim',[-.1 1.2],'FontSize',24);
    drawnow;
    
    subplot(3,2,2);
    set(Leg1_close,'XData',[pFoot1(1);pHip(1)],'YData',[pFoot1(2);pHip(2)]);
    set(Leg2_close,'XData',[pFoot2(1);pHip(1)],'YData',[pFoot2(2);pHip(2)]);
    set(COM_close,'XData',pCOM(1),'YData',pCOM(2));
    set(Leg1_close,'LineWidth',6,'Color',Leg1_color);
    set(Leg2_close,'LineWidth',6,'Color',Leg2_color);
    set(COM_close,'LineWidth',3,'MarkerSize',8','Color',COM_color);
    set(gca,'XLim',[pFoot1(1)-0.01 pFoot1(1)+0.01],'Ylim',[-.01 0.01],'FontSize',24);
    drawnow;
    
    subplot(3,2,4)
    set(scoord,'XData',name.tout(k),'YData',name.xout(k,4));
    set(wcoord,'XData',name.tout(k),'YData',name.xout(k,5));
    drawnow;
    
    subplot(3,2,6)
    set(dxs,'XData',name.tout(k),'YData',name.xout(k,6));
    drawnow;
    
    frame = getframe(mov);
    writeVideo(v,frame);
    
end
close(v);