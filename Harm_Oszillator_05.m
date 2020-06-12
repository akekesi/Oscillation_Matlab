%% Harmonischer Oszillator
%  2 x Pendel
%  1 x Feder
clear all
clc

global g;
global m1;
global m2;
global L1;
global L2;
global l1;
global l2;
global k;
global d
global x0
global x

vid = 0;        % als Video zu speichern: vid = 1

steps = 300;
h = 0.1;        % Zeitschritt
g = 9.81;       % Gravkonst.
m1 = 1;         % Gewicht 1
m2 = 1;         % Gewicht 2
L1 = 3;         % Laenge 1
L2 = 3;         % Laenge 2
l1 = 3*L1/5;    % Laenge_Feder 1
l2 = 4*L2/5;    % Laenge_Feder 2
k = 1;          % Federkoeff. 1
d = 4;          % Abstand zw. Pendel

w1 = 0.0/180*pi();  % Winkel_0 1
w2 = 45/180*pi();   % Winkel_0 2
dw1 = 0.0/180*pi(); % Winkelgeschw_0 1
dw2 = 0.0/180*pi(); % Winkelgeschw_0 2

x = @(l1,l2,w1,w2) sqrt((d+l2*sin(w2)-l1*sin(w1)).^2+(l2*cos(w2)-l1*cos(w1)).^2);
x0 = x(l1,l2,0,0);

y0 = [w1;w2;dw1;dw2];
Y = zeros(length(y0),steps);

for n = 1:1:steps
    if n == 1
        y = y0;
    else
        y = ruku(@abl,y,h);
    end
    Y(:,n) = y;
end

[A1,A2,C,omega_n] = analfunc(y0);
x1_a = @(t) A1(1,1)*(C(1,1)*cos(omega_n(1,1)*t)+C(2,1)*sin(omega_n(1,1)*t))+A2(1,1)*(C(3,1)*cos(omega_n(2,1)*t)+C(4,1)*sin(omega_n(2,1)*t));
x2_a = @(t) A1(2,1)*(C(1,1)*cos(omega_n(1,1)*t)+C(2,1)*sin(omega_n(1,1)*t))+A2(2,1)*(C(3,1)*cos(omega_n(2,1)*t)+C(4,1)*sin(omega_n(2,1)*t));

%% Plot
% Hilfswerte fuer Plot
xfeder = x(l1,l2,Y(1,:),Y(2,:));
xmax = max(xfeder-x0);
xmin = min(xfeder-x0);
xMax = max(xmax,abs(xmin));

wmax1 = max(Y(1,:));
wmin1 = min(Y(1,:));
wMax1 = max(wmax1,abs(wmin1));

wmax2 = max(Y(2,:));
wmin2 = min(Y(2,:));
wMax2 = max(wmax2,abs(wmin2));
wMax = max(wMax1,wMax2);

dwmax1 = max(Y(3,:));
dwmin1 = min(Y(3,:));
dwMax1 = max(dwmax1,abs(dwmin1));

dwmax2 = max(Y(4,:));
dwmin2 = min(Y(4,:));
dwMax2 = max(dwmax2,abs(dwmin2));
dwMax = max(dwMax1,dwMax2);

yMax = max(max(wMax,dwMax),xMax);
if yMax < 0.5
    yMax = fix(yMax*10)/10+0.2;
else
    yMax = fix(yMax*20)/10;
end

Epf = 1/2*k*(xfeder(1,:)-x0).^2;
Ep1 = m1*g*(L1-L1*cos(Y(1,:)));
Ep2 = m2*g*(L2-L2*cos(Y(2,:)));
Ek1 = 1/2*m1*(L1*Y(3,:)).^2;
Ek2 = 1/2*m2*(L2*Y(4,:)).^2;
E = Epf+Ep1+Ep2+Ek1+Ek2;
EMax = Epf(1,1)+Ep1(1,1)+Ep2(1,1)+Ek1(1,1)+Ek2(1,1);
if EMax < 1
    EMaxy = fix(EMax*10)/10+0.2;
else
    EMaxy = fix((Epf(1,1)+Ep1(1,1)+Ep2(1,1)+Ek1(1,1)+Ek2(1,1))*12)/10;
end

LW0 = 3;   % Linewidth 1
LWmin = 1.5;
LWd = (LW0-LWmin)/(xMax);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_05.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps

    subplot(2,2,1)
    hw = 0.1;
    geschw_d1 = 60/180*pi();
    geschw_dt1p = 0:hw:geschw_d1;
    geschw_dt1n = 0:-hw:-geschw_d1;
    geschw_1 = geschw_d1*Y(3,n)/dwMax1;
    geschw_y1 = 0.5;

    geschw_d2 = geschw_d1;
    geschw_dt2p = 0:hw:geschw_d2;
    geschw_dt2n = 0:-hw:-geschw_d2;
    geschw_2 = geschw_d2*Y(4,n)/dwMax2;
    geschw_y2 = 0.5;

    if geschw_1 > 0
        geschw_t1 = 0:hw:geschw_1;
        mark_geschw1 = '>';
    elseif geschw_1 < 0
        geschw_t1 = 0:-hw:geschw_1;
        mark_geschw1 = '<';
    else
        geschw_t1 = 0;
        mark_geschw1 = '.';
    end

    if geschw_2 > 0
        geschw_t2 = 0:hw:geschw_2;
        mark_geschw2 = '>';
    elseif geschw_2 < 0
        geschw_t2 = 0:-hw:geschw_2;
        mark_geschw2 = '<';
    else
        geschw_t2 = 0;
        mark_geschw2 = '.';
    end

    plot([-5 5],[0.4 0.4],'LineWidth',15,'Color','k')                                                                                                   % Wand
    hold on
    plot([-d/2 d/2],[0 0],'.','MarkerSize',12,'Color','k')                                                                                              % Aufhaeng. 1
    plot([-d/2 d/2],[0.17 0.17],'v','MarkerSize',12,'LineWidth',2,'Color','k')                                                                          % Aufhaeng. 2
    plot([-d/2 -d/2+L1*sin(Y(1,n))],[0 -L1*cos(Y(1,n))],'k')                                                                                            % Faden 1
    plot([d/2 d/2+L2*sin(Y(2,n))],[0 -L2*cos(Y(2,n))],'k')                                                                                              % Faden 2
    p13 = plot([-d/2+l1*sin(Y(1,n)) d/2+l2*sin(Y(2,n))],[-l1*cos(Y(1,n)) -l2*cos(Y(2,n))],'LineWidth',LW0-LWd*(xfeder(1,n)-x0),'Color','#77AC30');      % Feder
    plot([-d/2+l1*sin(Y(1,n)) d/2+l2*sin(Y(2,n))],[-l1*cos(Y(1,n)) -l2*cos(Y(2,n))],'.','MarkerSize',12,'Color','#77AC30')                              % Feder Aufhaeng.
    
    plot([-d/2+L1*sin(Y(1,n)) -d/2+(L1+geschw_y1)*sin(Y(1,n))],[-L1*cos(Y(1,n)) -(L1+geschw_y1)*cos(Y(1,n))],'--','Color',[0.7 0.7 0.7])
    plot(-d/2+L1*sin(Y(1,n))+geschw_y1*sin(Y(1,n)+geschw_dt1p),-L1*cos(Y(1,n))-geschw_y1*cos(Y(1,n)+geschw_dt1p),'LineWidth',3,'Color',[0.7 0.7 0.7])   % Center: Masse
    plot(-d/2+L1*sin(Y(1,n))+geschw_y1*sin(Y(1,n)+geschw_dt1n),-L1*cos(Y(1,n))-geschw_y1*cos(Y(1,n)+geschw_dt1n),'LineWidth',3,'Color',[0.7 0.7 0.7])   % Center: Masse
    p12 = plot(-d/2+L1*sin(Y(1,n))+geschw_y1*sin(Y(1,n)+geschw_t1),-L1*cos(Y(1,n))-geschw_y1*cos(Y(1,n)+geschw_t1),'LineWidth',3,'Color','#D95319');    % Center: Masse
    
    plot([d/2+L2*sin(Y(2,n)) d/2+(L2+geschw_y2)*sin(Y(2,n))],[-L2*cos(Y(2,n)) -(L2+geschw_y2)*cos(Y(2,n))],'--','Color',[0.7 0.7 0.7])
    plot(d/2+L2*sin(Y(2,n))+geschw_y2*sin(Y(2,n)+geschw_dt2p),-L2*cos(Y(2,n))-geschw_y2*cos(Y(2,n)+geschw_dt2p),'LineWidth',3,'Color',[0.7 0.7 0.7])    % Center: Masse
    plot(d/2+L2*sin(Y(2,n))+geschw_y2*sin(Y(2,n)+geschw_dt2n),-L2*cos(Y(2,n))-geschw_y2*cos(Y(2,n)+geschw_dt2n),'LineWidth',3,'Color',[0.7 0.7 0.7])    % Center: Masse
    plot(d/2+L2*sin(Y(2,n))+geschw_y2*sin(Y(2,n)+geschw_t2),-L2*cos(Y(2,n))-geschw_y2*cos(Y(2,n)+geschw_t2),'LineWidth',3,'Color','#D95319')            % Center: Masse

    p11 = plot([-d/2+L1*sin(Y(1,n)) d/2+L2*sin(Y(2,n))],[-L1*cos(Y(1,n)) -L2*cos(Y(2,n))],'.','MarkerSize',35,'Color','k');                             % Masse 1-2

    xlim([-5 5])
    ylim([-5 0.55])
    yticks(-5:1:1)
    legend([p11 p12 p13],{'Masse','Geschw.','Feder'},'location','SouthEast')
    text(-4.4,-3.5,'\alpha_0 =')
    text(-3.9,-3.5,{[num2str(w1/pi*180) '\circ'] [num2str(w2/pi*180) '\circ']})
    text(-4.4,-4.5,'\beta_0 =')
    text(-3.9,-4.5,{[num2str(dw1/pi*180) '\circ/s'] [num2str(dw2/pi*180) '\circ/s']})
    title('Animation','FontSize',15,'FontWeight','normal')
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,2,2)
    p25 = plot([0 0],[0 Ep2(1,n)],'LineWidth',150,'Color',[0 0.2470 0.5410]);
    hold on
    p24 = plot([0 0],[Ep2(1,n) Ep2(1,n)+Ek2(1,n)],'LineWidth',150,'Color',[0.8500 0.3250 0.0980]);
    p23 = plot([0 0],[Ep2(1,n)+Ek2(1,n) Ep2(1,n)+Ek2(1,n)+Epf(1,n)],'LineWidth',150,'Color',[0 0.4470 0.7410]);
    p22 = plot([0 0],[Ep2(1,n)+Ek2(1,n)+Epf(1,n) Ep2(1,n)+Ek2(1,n)+Epf(1,n)+Ek1(1,n)],'LineWidth',150,'Color',[1.0 0.5250 0.2980]);
    p21 = plot([0 0],[Ep2(1,n)+Ek2(1,n)+Epf(1,n)+Ek1(1,n) Ep2(1,n)+Ek2(1,n)+Epf(1,n)+Ek1(1,n)+Ep1(1,n)],'LineWidth',150,'Color',[0.2 0.6470 0.9410]);
    plot([-EMax EMax],[EMax EMax],'--k')
    xlim([-EMax EMax])
    xticks([])
    ylim([0 EMaxy])
    legend([p21 p22 p23 p24 p25],{'E_{Pot1}','E_{Kin1}','E_{PotF}','E_{Kin2}','E_{Pot2}'},'Location','southeast')
    text(-0.7*EMax,EMax*1.05,'E_{Summe}')
    title('Energie','FontSize',16,'FontWeight','normal')
    ylabel('E [J]')
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,2,3)
    t=0:h:h*(n-1);
    plot(t,Y(1,1:n),'Color','#0072BD')
    hold on
    plot(t,Y(2,1:n),'--','Color','#0072BD')
    plot(t,Y(3,1:n),'Color','#D95319')
    plot(t,Y(4,1:n),'--','Color','#D95319')
    plot(t,xfeder(1,1:n)-x0,'Color','#77AC30')
    xlim([0 h*steps])
    ylim([-yMax yMax])
    xticks(0:5:h*steps)
    legend({'$\alpha_1$','$\alpha_2$','$\beta_1$','$\beta_2$','$\Delta x$'},'Interpreter', 'latex','location','NorthEast')
    title('$\alpha-t \,\,/\,\, \beta-t \,\,/\,\, F_{Feder}-t$','Interpreter', 'latex','FontSize',16,'FontWeight','normal')
    xlabel('$t \, [s]$','Interpreter', 'latex')
    ylabel('$\alpha \, [rad], \,\, \beta \, [rad/s], \,\, F_{Feder} \, [N]$','Interpreter', 'latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
%    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,2,4)
    t=0:h:h*(n-1);
    plot(t,Y(1,1:n),'Color','#0072BD')
    hold on
    plot(t,Y(2,1:n),'Color','#0072BD')
    plot(t,x1_a(t),'--','Color','k');
    plot(t,x2_a(t),'--','Color','k');
    xlim([0 h*steps])
    ylim([-yMax yMax])
    xticks(0:5:h*steps)
    legend({'$\alpha_1 \, (RK4)$','$\alpha_2 \, (RK4)$','$\alpha_1 \, (ana)$','$\alpha_2 \, (ana)$'},'Interpreter', 'latex','location','NorthEast')
    text(0.05*steps*h,0.75*yMax,'A_1 =')
    text(0.1*steps*h,0.75*yMax,num2str(A1))
    text(0.15*steps*h,0.75*yMax,'A_2 =')
    text(0.2*steps*h,0.75*yMax,num2str(A2))
    title('$\alpha-t \,\, (num \,/\, ana)$','Interpreter', 'latex','FontSize',16,'FontWeight','normal')
    xlabel('$t \, [s]$','Interpreter', 'latex')
    ylabel('$\alpha \, [rad]$','Interpreter', 'latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
%    daspect([1 1 1])
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

%% Abl.
function dy = abl(y)
global g;
global m1;
global m2;
global L1;
global L2;
global l1;
global l2;
global k;
global d;
global x0
global x;
w1 = y(1,1);
w2 = y(2,1);
dw1 = y(3,1);
dw2 = y(4,1);
xfeder = x(l1,l2,w1,w2);

ddw1 = -m1*g*L1*sin(w1)/(m1*L1^2)-k*(xfeder-x0)/(m1*L1^2)*((d+l2*sin(w2)-l1*sin(w1))*(-l1*cos(w1))+(l2*cos(w2)-l1*cos(w1))*(l1*sin(w1)))/(xfeder);
ddw2 = -m2*g*L2*sin(w2)/(m2*L2^2)-k*(xfeder-x0)/(m2*L2^2)*((d+l2*sin(w2)-l1*sin(w1))*(l2*cos(w2))+(l2*cos(w2)-l1*cos(w1))*(-l2*sin(w2)))/(xfeder);
dy = [dw1;dw2;ddw1;ddw2];
end

%% Runge-Kutta
function y_neu = ruku(funk,y,h)
k1 = funk(y);
k2 = funk(y+h/2*k1);
k3 = funk(y+h/2*k2);
k4 = funk(y+h*k3);
y_neu = y+h*(k1/6+k2/3+k3/3+k4/6);
end

%% analytische Funktion mit Kleinwinkelnaeherung (1.Grad)
function [A1,A2,C,omega_n] = analfunc(y0)
global g;
global m1;
global m2;
global L1;
global L2;
global l1;
global l2;
global k;
omega_n = sqrt(roots([m1*m2*L1^2*L2^2 -((m1*g*L1+k*l1*l2)*m2*L2^2+(m2*g*L2+k*l1*l2)*m1*L1^2) (m1*g*L1+k*l1*l2)*(m2*g*L2+k*l1*l2)-k^2*l1^2*l2^2]));
a1_1 = 1;
a1_2 = (-m1*L1^2*omega_n(1,1)^2+m1*g*L1+k*l1*l2)/(k*l1*l2);
a2_1 = 1;
a2_2 = (-m1*L1^2*omega_n(2,1)^2+m1*g*L1+k*l1*l2)/(k*l1*l2);
A1 = [a1_1;a1_2];
A2 = [a2_1;a2_2];
C = [a1_1 0 a2_1 0;a1_2 0 a2_2 0;0 a1_1*omega_n(1,1) 0 a2_1*omega_n(2,1);0 a1_2*omega_n(1,1) 0 a2_2*omega_n(2,1)]\y0;
end