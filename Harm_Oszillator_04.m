%% Harmonischer Oszillator
%  horizontal
%  2 x Punktmasse
%  3 x Feder
clear all
clc

global m1;
global m2;
global k1;
global k2;
global k3;

vid = 0;    % als Video zu speichern: vid = 1

steps = 150;
h = 0.1;    % Zeitschritt
m1 = 2;     % Gewicht 1
m2 = 1;     % Gewicht 2
k1 = 2;     % Federkoeff. 1
k2 = 1;     % Federkoeff. 2
k3 = 1;     % Federkoeff. 3

x1 = 0;     % Anfangspos. 1
x2 = 0;     % Anfangspos. 2
v1 = 0;     % Anfangsgeschw. 1
v2 = 1;     % Anfangsgeschw. 2

y0 = [x1;x2;v1;v2];
Y = zeros(length(y0),steps);

for n = 1:1:steps
    if n == 1
        y = y0;
    else
        y = ruku(@abl,y,h);
    end
    Y(:,n) = y;
end

[A1,A2,C,omega_n] = anafunc(y0);
x1_a = @(t) A1(1,1)*(C(1,1)*cos(omega_n(1,1)*t)+C(2,1)*sin(omega_n(1,1)*t))+A2(1,1)*(C(3,1)*cos(omega_n(2,1)*t)+C(4,1)*sin(omega_n(2,1)*t));
x2_a = @(t) A1(2,1)*(C(1,1)*cos(omega_n(1,1)*t)+C(2,1)*sin(omega_n(1,1)*t))+A2(2,1)*(C(3,1)*cos(omega_n(2,1)*t)+C(4,1)*sin(omega_n(2,1)*t));

%% Plot
% Hilfswerte fuer Plot
maxvx = max(Y,[],(2));
minvx = min(Y,[],(2));

xmax1 = maxvx(1,1);
xmin1 = minvx(1,1);
xMax1 = max(xmax1,abs(xmin1));
L01 = fix(30*xMax1)/10;
if L01 < 1
    L01 = 1;
end

xmax2 = max(Y(2,:)-Y(1,:));
xmin2 = min((Y(2,:)-Y(1,:)));
xMax2 = max(xmax2,abs(xmin2));
L02 = fix(30*xMax2)/10;
if L02 < 1
    L02 = 1;
end
d = L02+1;  % Abstand zw. Punktmassen

L03 = fix(30*xMax2)/10;
if L03 < 1
    L03 = 1;
end

xMax = max(xMax1,xMax2);

vmax1 = maxvx(3,1);
vmin1 = minvx(3,1);
vMax1 = max(vmax1,abs(vmin1));

vmax2 = maxvx(4,1);
vmin2 = minvx(4,1);
vMax2 = max(vmax2,abs(vmin2));

vMax = max(vMax1,vMax2);

kraft1 = -1*Y(1,:)*k1+(Y(2,:)-Y(1,:))*k2;
kraftMax1 = max(abs(kraft1));
kraft2 = -(Y(2,:)-Y(1,:))*k2-Y(2,:)*k3;
kraftMax2 = max(abs(kraft2));

yMax = max(max(xMax,vMax),max(kraftMax1,kraftMax2));

Ep1 = 1/2*k1*Y(1,:).^2;
Ep2 = 1/2*k2*(Y(2,:)-Y(1,:)).^2;
Ep3 = 1/2*k3*Y(2,:).^2;
Ek1 = 1/2*m1*Y(3,:).^2;
Ek2 = 1/2*m2*Y(4,:).^2;
EMax = Ep1(1,1)+Ep2(1,1)+Ep3(1,1)+Ek1(1,1)+Ek2(1,1);

LW01 = 3;   % Linewidth 1
LWmin1 = 1.5;
LWd1 = (LW01-LWmin1)/(xMax1);

LW02 = 3;   % Linewidth 2
LWmin2 = 1.5;
LWd2 = (LW02-LWmin2)/(xMax2);

LW03 = 3;   % Linewidth 3
LWmin3 = 1.5;
LWd3 = (LW03-LWmin3)/(xMax2);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_04.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps
    
    subplot(2,2,1)
    geschw_dx1 = 0.5;
    geschw_x1 = geschw_dx1*Y(3,n)/vMax1;
    geschw_y1 = -0.3;
    if geschw_x1 > 0
        mark_geschw1 = '>';
    elseif geschw_x1 < 0
        mark_geschw1 = '<';
    else
        mark_geschw1 = '.';
    end

    geschw_dx2 = geschw_dx1;
    geschw_x2 = geschw_dx2*Y(4,n)/vMax2;
    geschw_y2 = geschw_y1;
    if geschw_x2 > 0
        mark_geschw2 = '>';
    elseif geschw_x2 < 0
        mark_geschw2 = '<';
    else
        mark_geschw2 = '.';
    end

    kraft_f_dx1 = 0.5;
    kraft_f_x1 = kraft_f_dx1*kraft1(1,n)/kraftMax1;
    kraft_f_y1 = geschw_y1-0.3;
    if kraft_f_x1 > 0
        mark_f_kraft1 = '>';
    elseif kraft_f_x1 < 0
        mark_f_kraft1 = '<';
    else
        mark_f_kraft1 = '.';
    end

    kraft_f_dx2 = kraft_f_dx1;
    kraft_f_x2 = kraft_f_dx2*kraft2(1,n)/kraftMax2;
    kraft_f_y2 = kraft_f_y1;
    if kraft_f_x2 > 0
        mark_f_kraft2 = '>';
    elseif kraft_f_x2 < 0
        mark_f_kraft2 = '<';
    else
        mark_f_kraft2 = '.';
    end

    p11 = plot(Y(1,1),0,'o','MarkerSize',10,'LineWidth',2,'Color','k');                                     % Anf. 1
    hold on
    plot(Y(2,1)+d,0,'o','MarkerSize',10,'LineWidth',2,'Color','k')                                          % Anf. 2
    p15 = plot([-L01 Y(1,n)],[0 0],'LineWidth',LW01-Y(1,n)*LWd1,'Color','#77AC30');                         % Feder 1
    plot([Y(1,n) Y(2,n)+d],[0 0],'LineWidth',LW02-(Y(2,n)-Y(1,n))*LWd2,'Color','#77AC30')                   % Feder 2
    plot([Y(2,n)+d d+L03],[0 0],'LineWidth',LW03+Y(2,n)*LWd3,'Color','#77AC30')                             % Feder 3
    plot([-L01 -L01],[-10 10],'LineWidth',5,'Color','k')                                                    % Wand 1
    plot([d+L03 d+L03],[-10 10],'LineWidth',5,'Color','k')                                                  % Wand 2

    plot([Y(1,n) Y(1,n)],[0 kraft_f_y1],'--','Color',[0.7 0.7 0.7])
    plot(Y(1,n)-geschw_dx1,geschw_y1,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+geschw_dx1,geschw_y1,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-geschw_dx1 Y(1,n)+geschw_dx1],[geschw_y1 geschw_y1],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p13 = plot([Y(1,n) Y(1,n)+geschw_x1],[geschw_y1 geschw_y1],'LineWidth',3,'Color','#D95319');            % Geschw. 1
    plot(Y(1,n)+geschw_x1,geschw_y1,'Marker',mark_geschw1,'LineWidth',3,'Color','#D95319')                  % Geschw. 1

    plot(Y(1,n)-kraft_f_dx1,kraft_f_y1,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+kraft_f_dx1,kraft_f_y1,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-kraft_f_dx1 Y(1,n)+kraft_f_dx1],[kraft_f_y1 kraft_f_y1],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p14 = plot([Y(1,n) Y(1,n)+kraft_f_x1],[kraft_f_y1 kraft_f_y1],'LineWidth',3,'Color','#7E2F8E');         % Kraft_Feder 1
    plot(Y(1,n)+kraft_f_x1,kraft_f_y1,'Marker',mark_f_kraft1,'LineWidth',3,'Color','#7E2F8E')               % Kraft_Feder 1

    plot([Y(2,n)+d Y(2,n)+d],[0 geschw_y2],'--','Color',[0.7 0.7 0.7])
    plot(Y(2,n)+d-geschw_dx2,geschw_y2,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(2,n)+d+geschw_dx2,geschw_y2,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(2,n)+d-geschw_dx2 Y(2,n)+d+geschw_dx2],[geschw_y2 geschw_y2],'LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(2,n)+d Y(2,n)+d+geschw_x2],[geschw_y2 geschw_y2],'LineWidth',3,'Color','#D95319')               % Geschw. 2
    plot(Y(2,n)+d+geschw_x2,geschw_y2,'Marker',mark_geschw2,'LineWidth',3,'Color','#D95319')                % Geschw. 2

    plot(Y(2,n)+d-kraft_f_dx2,kraft_f_y2,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(2,n)+d+kraft_f_dx2,kraft_f_y2,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(2,n)+d-kraft_f_dx2 Y(2,n)+d+kraft_f_dx2],[kraft_f_y2 kraft_f_y2],'LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(2,n)+d Y(2,n)+d+kraft_f_x2],[kraft_f_y2 kraft_f_y2],'LineWidth',3,'Color','#7E2F8E')            % Kraft_Feder 2
    plot(Y(2,n)+d+kraft_f_x2,kraft_f_y2,'Marker',mark_f_kraft2,'LineWidth',3,'Color','#7E2F8E')             % Kraft_Feder 2

    p12 = plot(Y(1,n),0,'.','MarkerSize',30,'Color','k');                                                   % Gewicht 1
    plot(Y(2,n)+d,0,'.','MarkerSize',30,'Color','k')                                                        % Gewicht 2
    xlim([-L01 d+L03])
    xticks(-fix(L01):1:fix(d+L03))
    ylim([-1 1.5])
    yticks([])
    legend([p11 p12 p13 p14 p15],{'Anfang','Masse','Geschw','F_{Feder}','Feder'},'location','NorthEast')
    text(0.5,0.95,{['m_1=',num2str(m1),'kg     k_1=',num2str(k1),'N/m     x_{10}=',num2str(x1),'m     v_{10}=',num2str(v1),'m/s'],['m_2=',num2str(m2),'kg     k_2=',num2str(k2),'N/m     x_{20}=',num2str(x2),'m     v_{20}=',num2str(v2),'m/s'],['                  k_3=',num2str(k3),'N/m']})
    title('Animation','FontSize',16,'FontWeight','normal')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,2,2)
    p25 = plot([0 0],[0 Ep3(1,n)],'LineWidth',150,'Color',[0 0.2470 0.5410]);
    hold on
    p24 = plot([0 0],[Ep3(1,n) Ep3(1,n)+Ek2(1,n)],'LineWidth',150,'Color',[0.8500 0.3250 0.0980]);
    p23 = plot([0 0],[Ep3(1,n)+Ek2(1,n) Ep3(1,n)+Ek2(1,n)+Ep2(1,n)],'LineWidth',150,'Color',[0 0.4470 0.7410]);
    p22 = plot([0 0],[Ep3(1,n)+Ek2(1,n)+Ep2(1,n) Ep3(1,n)+Ek2(1,n)+Ep2(1,n)+Ek1(1,n)],'LineWidth',150,'Color',[1.0 0.5250 0.2980]);
    p21 = plot([0 0],[Ep3(1,n)+Ek2(1,n)+Ep2(1,n)+Ek1(1,n) Ep3(1,n)+Ek2(1,n)+Ep2(1,n)+Ek1(1,n)+Ep1(1,n)],'LineWidth',150,'Color',[0.2 0.6470 0.9410]);
    plot([-0.5 0.5],[EMax EMax],'--k')
    xlim([-0.4 0.4])
    xticks([])
    ylim([0 fix(EMax*12)/10])
    legend([p21 p22 p23 p24 p25],{'E_{Pot1}','E_{Kin1}','E_{Pot2}','E_{Kin2}','E_{Pot3}'},'Location','southeast')
    text(-0.45,EMax+0.04,'E_{Summe}')
    title('Energie','FontSize',16,'FontWeight','normal')
    ylabel('E [J]')
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,2,[3 4])
    t=0:h:h*(n-1);
    p31 = plot(t,Y(1,1:n),'Color','#0072BD');
    hold on
    p32 = plot(t,Y(2,1:n),'--','Color','#0072BD');
    p33 = plot(t,Y(3,1:n),'Color','#D95319');
    p34 = plot(t,Y(4,1:n),'--','Color','#D95319');
    p35 = plot(t,kraft1(1,1:n),'Color','#7E2F8E');
    p36 = plot(t,kraft2(1,1:n),'--','Color','#7E2F8E');
%    plot(t,x1_a(t),'x','MarkerSize',5,'Color','k')
%    plot(t,x2_a(t),'x','MarkerSize',5,'Color','k')
    xlim([0 h*steps])
    ylim([fix(-yMax*2) fix(yMax*3)])
    legend([p31 p32 p33 p34 p35 p36],{'x_1','x_2','v_1','v_2','F_{Feder 1}','F_{Feder 2}'},'location','NorthEast')
    text(0.7,fix(yMax*2),'A_1 =')
    text(1.4,fix(yMax*2),num2str(A1))
    text(2.1,fix(yMax*2),'A_2 =')
    text(2.8,fix(yMax*2),num2str(A2))
    text(3.5,fix(yMax*2),'\omega_n^2 =')
    text(4.2,fix(yMax*2),num2str(omega_n.^2))
    title('x-t / v-t / F_{Feder}-t','FontSize',16,'FontWeight','normal')
    xlabel('t [s]')
    ylabel('x [m],  v [m/s],  F_{Feder} [N]')
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
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
function [dy] = abl(y)
global m1;
global m2;
global k1;
global k2;
global k3
x1 = y(1,1);
x2 = y(2,1);
v1 = y(3,1);
v2 = y(4,1);

a1 = -(k1+k2)/m1*x1+k2/m1*x2;
a2 = k2/m2*x1-(k2+k3)/m2*x2;
dy = [v1;v2;a1;a2];
end

%% Runge-Kutta
function [y_neu] = ruku(funk,y,h)
k1 = funk(y);
k2 = funk(y+h/2*k1);
k3 = funk(y+h/2*k2);
k4 = funk(y+h*k3);
y_neu = y+h*(k1/6+k2/3+k3/3+k4/6);
end

%% analytische Funktion
function [A1,A2,C,omega_n] = anafunc(y0)
global m1;
global m2;
global k1;
global k2;
global k3
omega_n = sqrt(roots([1 -(k1+k2)/m1-(k2+k3)/m2 (k1*k2+k1*k3+k2*k3)/(m1*m2)]));
a1_1 = 1;
a1_2 = (-omega_n(1,1)^2*m1+k1+k2)/k2;
a2_1 = 1;
a2_2 = (-omega_n(2,1)^2*m1+k1+k2)/k2;
A1 = [a1_1;a1_2];
A2 = [a2_1;a2_2];
C = [a1_1 0 a2_1 0;a1_2 0 a2_2 0;0 a1_1*omega_n(1,1) 0 a2_1*omega_n(2,1);0 a1_2*omega_n(1,1) 0 a2_2*omega_n(2,1)]\y0;
end