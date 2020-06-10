%% Harmonischer Oszillator
%  horizontal
%  2 x Punktmasse
%  2 x Feder
clear all
clc

global m1;
global m2;
global k1;
global k2;

vid = 1;    % als Video zu speichern: vid = 1

steps = 220;
h = 0.1;    % Zeitschritt
m1 = 2;     % Gewicht 1
m2 = 1;     % Gewicht 2
k1 = 2;     % Federkoeff. 1
k2 = 1;     % Federkoeff. 2

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

[A1,A2,C,omega_n] = analitisch(y0);
x1_a = @(t) A1(1,1)*(C(1,1)*cos(omega_n(1,1)*t)+C(2,1)*sin(omega_n(1,1)*t))+A2(1,1)*(C(3,1)*cos(omega_n(2,1)*t)+C(4,1)*sin(omega_n(2,1)*t));
x2_a = @(t) A1(2,1)*(C(1,1)*cos(omega_n(1,1)*t)+C(2,1)*sin(omega_n(1,1)*t))+A2(2,1)*(C(3,1)*cos(omega_n(2,1)*t)+C(4,1)*sin(omega_n(2,1)*t));;

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
kraft2 = -(Y(2,:)-Y(1,:))*k2;
kraftMax2 = max(abs(kraft2));

yMax = max(max(xMax,vMax),max(kraftMax1,kraftMax2));

LW01 = 3;   % Linewidth 1
LWmin1 = 1.5;
LWd1 = (LW01-LWmin1)/(xMax1);

LW02 = 3;   % Linewidth 2
LWmin2 = 1.5;
LWd2 = (LW02-LWmin2)/(xMax2);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_03.avi');
    Video.FrameRate = 13;
    open(Video)
end

for n = 1:1:steps
    
    subplot(2,1,1)
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

    p11 = plot(Y(1,1),0,'o','MarkerSize',10,'LineWidth',2,'Color','k');                             % Anf. 1
    hold on
    p12 = plot(Y(2,1)+d,0,'o','MarkerSize',10,'LineWidth',2,'Color','k');                                   % Anf. 1
    p15 = plot([-L01 Y(1,n)],[0 0],'LineWidth',LW01-Y(1,n)*LWd1,'Color','#77AC30');                         % Feder 1
    p16 = plot([Y(1,n) Y(2,n)+d],[0 0],'LineWidth',LW02-(Y(2,n)-Y(1,n))*LWd2,'Color','#77AC30');            % Feder 2
    plot([-L01 -L01],[-10 10],'LineWidth',5,'Color','k')                                                    % Wand

    plot([Y(1,n) Y(1,n)],[0 kraft_f_y1],'--','Color',[0.7 0.7 0.7])
    plot(Y(1,n)-geschw_dx1,geschw_y1,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+geschw_dx1,geschw_y1,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-geschw_dx1 Y(1,n)+geschw_dx1],[geschw_y1 geschw_y1],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p17 = plot([Y(1,n) Y(1,n)+geschw_x1],[geschw_y1 geschw_y1],'LineWidth',3,'Color','#D95319');            % Geschw. 1
    plot(Y(1,n)+geschw_x1,geschw_y1,'Marker',mark_geschw1,'LineWidth',3,'Color','#D95319');                 % Geschw. 1

    plot(Y(1,n)-kraft_f_dx1,kraft_f_y1,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+kraft_f_dx1,kraft_f_y1,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-kraft_f_dx1 Y(1,n)+kraft_f_dx1],[kraft_f_y1 kraft_f_y1],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p18 = plot([Y(1,n) Y(1,n)+kraft_f_x1],[kraft_f_y1 kraft_f_y1],'LineWidth',3,'Color','#7E2F8E');         % Kraft_Feder 1
    plot(Y(1,n)+kraft_f_x1,kraft_f_y1,'Marker',mark_f_kraft1,'LineWidth',3,'Color','#7E2F8E');              % Kraft_Feder 1

    plot([Y(2,n)+d Y(2,n)+d],[0 geschw_y2],'--','Color',[0.7 0.7 0.7])
    plot(Y(2,n)+d-geschw_dx2,geschw_y2,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(2,n)+d+geschw_dx2,geschw_y2,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(2,n)+d-geschw_dx2 Y(2,n)+d+geschw_dx2],[geschw_y2 geschw_y2],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p19 = plot([Y(2,n)+d Y(2,n)+d+geschw_x2],[geschw_y2 geschw_y2],'LineWidth',3,'Color','#D95319');        % Geschw. 1
    plot(Y(2,n)+d+geschw_x2,geschw_y2,'Marker',mark_geschw2,'LineWidth',3,'Color','#D95319');               % Geschw. 1

    plot(Y(2,n)+d-kraft_f_dx2,kraft_f_y2,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(2,n)+d+kraft_f_dx2,kraft_f_y2,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(2,n)+d-kraft_f_dx2 Y(2,n)+d+kraft_f_dx2],[kraft_f_y2 kraft_f_y2],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p110 = plot([Y(2,n)+d Y(2,n)+d+kraft_f_x2],[kraft_f_y2 kraft_f_y2],'LineWidth',3,'Color','#7E2F8E');    % Kraft_Feder 2
    plot(Y(2,n)+d+kraft_f_x2,kraft_f_y2,'Marker',mark_f_kraft2,'LineWidth',3,'Color','#7E2F8E');            % Kraft_Feder 2

    p13 = plot(Y(1,n),0,'.','MarkerSize',30,'Color','k');                                                   % Gewicht 1
    p14 = plot(Y(2,n)+d,0,'.','MarkerSize',30,'Color','k');                                                 % Gewicht 2
    xlim([-L01 d+L02])
    xticks(-fix(L01):1:fix(d+L02))
    ylim([-1 1.5])
    yticks([])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    legend([p11 p13 p17 p18 p15],{'Anfang','Masse','Geschw','F_{Feder}','Feder'},'location','NorthEast')
%    title({['x_0=',num2str(Y(1,1)),',   v_0 = ',num2str(Y(2,1 ))],['m = ',num2str(m),',   k = ',num2str(k)]},'FontSize',15,'FontWeight','normal')
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,1,2)
    t=0:h:h*(n-1);
    p21 = plot(t,Y(1,1:n),'Color','#0072BD');
    hold on
    p22 = plot(t,Y(2,1:n),'--','Color','#0072BD');
    p23 = plot(t,Y(3,1:n),'Color','#D95319');
    p24 = plot(t,Y(4,1:n),'--','Color','#D95319');
    p25 = plot(t,kraft1(1,1:n),'Color','#7E2F8E');
    p26 = plot(t,kraft2(1,1:n),'--','Color','#7E2F8E');
%    plot(t,x1_a(t),'x','MarkerSize',5,'Color','k')
%    plot(t,x2_a(t),'x','MarkerSize',5,'Color','k')
    legend([p21 p22 p23 p24 p25 p26],{'x_1','x_2','v_1','v_2','F_{Feder 1}','F_{Feder 2}'},'location','NorthEast')
    text(2,fix(yMax*2),'A_1 = ')
    text(2.6,fix(yMax*2),num2str(A1))
    text(4,fix(yMax*2),'A_2 = ')
    text(4.6,fix(yMax*2),num2str(A2))
    xlim([0 h*steps])
    ylim([fix(-yMax*2) fix(yMax*3)])
    xlabel('t')
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
x1 = y(1,1);
x2 = y(2,1);
v1 = y(3,1);
v2 = y(4,1);

a1 = -(k1+k2)/m1*x1+k2/m1*x2;
a2 = k2/m2*x1-k2/m2*x2;
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

%% analitische Funktion
function [A1,A2,C,omega_n] = analitisch(y0)
global m1;
global m2;
global k1;
global k2;
omega_n = sqrt(roots([1 -(k1+k2)/m1-k2/m2 k1*k2/(m1*m2)]));
a1_1 = 1;
a1_2 = (-omega_n(1,1)^2*m1+k1+k2)/k2;
a2_1 = 1;
a2_2 = (-omega_n(2,1)^2*m1+k1+k2)/k2;
A1 = [a1_1;a1_2];
A2 = [a2_1;a2_2];
C = [a1_1 0 a2_1 0;a1_2 0 a2_2 0;0 a1_1*omega_n(1,1) 0 a2_1*omega_n(2,1);0 a1_2*omega_n(1,1) 0 a2_2*omega_n(2,1)]\y0;
end