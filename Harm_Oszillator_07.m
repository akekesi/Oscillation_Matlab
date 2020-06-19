%% Harmonischer Oszillator
% 2 x Feder
% 1 x Pendel
% 2 x Masse
clear all
clc

global g;
global m1;
global m2;
global k1;
global k2;
global L;

vid = 1;        % als Video zu speichern: vid = 1

steps = 200;
h = 0.05;       % Zeitschritt
g = 9.81;       % Gravkonst.
m1 = 1;         % Gewicht 1
m2 = 1;         % Gewicht 2
k1 = 10;        % Federkoeff. 1
k2 = k1;        % Federkoeff. 2
L = 1;          % Laenge 1

x = 1;          % Position_0
w = 45/180*pi;  % Winkel_0
dx = 0;         % Geschw_0
dw = 0;         % Winkelgeschw_0

y0 = [x;w;dx;dw];
Y = zeros(length(y0),steps);
E = zeros(5,1);

for n = 1:1:steps
    if n == 1
        y = y0;
    else
        y = ruku(@abl,y,h); 
    end
    Y(:,n) = y;
    E(:,n) = energie(y);
end

%% Plot
% Hilfswerte fuer Plot
xmax = max(Y(1,:));
xmin = min(Y(1,:));
xMax = max(xmax,abs(xmin));
L0 = fix(xMax*20)/10;

wmax = max(Y(2,:));
wmin = min(Y(2,:));
wMax = max(wmax,abs(wmin));

dxmax = max(Y(3,:));
dxmin = min(Y(3,:));
dxMax = max(dxmax,abs(dxmin));

dwmax = max(Y(4,:));
dwmin = min(Y(4,:));
dwMax = max(dwmax,abs(dwmin));

ymax = max(max(xmax,wmax),max(dxmax,dwmax));
ymin = min(min(xmin,wmin),min(dxmin,dwmin));

EMax = E(1,1)+E(2,1)+E(3,1)+E(4,1)+E(5,1);

LW0 = 3;   % Linewidth 1
LWmin = 1.5;
LWd = (LW0-LWmin)/(xMax);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_07.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps
    x1 = Y(1,n);
    y1 = 0;
    x2 = x1+L*sin(Y(2,n));
    y2 = -L*cos(Y(2,n));

    subplot(2,2,1)
    p13 = plot([-L0 x1],[y1 y1],'LineWidth',LW0-LWd*x1,'Color','#77AC30');  % Feder 1
    hold on
    plot([x1 L0],[y1 y1],'LineWidth',LW0+LWd*x1,'Color','#77AC30')          % Feder 2
    plot([-L0 -L0],[-10 10],'LineWidth',10,'Color','k')                     % Wand
    plot([L0 L0],[-10 10],'LineWidth',10,'Color','k')                       % Wand
    p11 = plot(x1,y1,'.','MarkerSize',35,'Color','k');                      % Masse 1
    p12 = plot(x2,y2,'.','MarkerSize',35,'Color','k');                      % Masse 2
    plot([x1 x2],[y1 y2],'k')
    xlim([-L0 L0])
    ylim([-L*1.5 0.5])
    legend([p11 p13],{'Masse','Feder'},'location','SouthEast')
    title('Animation','FontSize',15,'FontWeight','normal')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,2,2)
    plot([-10 10],[EMax EMax],'--k')
    hold on
    p21 = plot([1 1],[0 E(1,n)],'LineWidth',10,'Color',[0.4660 0.6740 0.1880]);
    plot([1 2],[E(1,n) E(1,n)],'--k')
    p22 = plot([2 2],[E(1,n) E(1,n)+E(2,n)],'LineWidth',10,'Color',[0.6660 0.8740 0.3880]);
    plot([2 3],[E(1,n)+E(2,n) E(1,n)+E(2,n)],'--k')
    p23 = plot([3 3],[E(1,n)+E(2,n) E(1,n)+E(2,n)+E(3,n)],'LineWidth',10,'Color','#7E2F8E');
    plot([3 4],[E(1,n)+E(2,n)+E(3,n) E(1,n)+E(2,n)+E(3,n)],'--k')
    p24 = plot([4 4],[E(1,n)+E(2,n)+E(3,n) E(1,n)+E(2,n)+E(3,n)+E(4,n)],'LineWidth',10,'Color',[0 0.4470 0.7410]);
    plot([4 5],[E(1,n)+E(2,n)+E(3,n)+E(4,n) E(1,n)+E(2,n)+E(3,n)+E(4,n)],'--k')
    p2f = plot([5 5],[E(1,n)+E(2,n)+E(3,n)+E(4,n)+E(5,n) EMax],'LineWidth',25,'Color',[1 0 0]);
    p25 = plot([5 5],[E(1,n)+E(2,n)+E(3,n)+E(4,n) E(1,n)+E(2,n)+E(3,n)+E(4,n)+E(5,n)],'LineWidth',10,'Color',[0.2 0.6470 0.9410]);
    xlim([0 6.5])
    xticks([])
    ylim([0 EMax*1.5])
    legend([p21 p22 p23 p24 p25 p2f],{'E_{Pot-FL}','E_{Pot-FR}','E_{Pot}','E_{Kin-O}','E_{Kin-U}','Fehler'},'Location','southeast')
    text(0.5,EMax*1.05,'E_{Summe}')
    title('Energie','FontSize',16,'FontWeight','normal')
    ylabel('E [J]')
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
    drawnow
    hold off

    subplot(2,2,[3 4])
    t = 0:h:h*(n-1);
    plot(t,Y(1,1:n),'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
    hold on
    plot(t,Y(2,1:n),'LineWidth',2,'Color',[1.0 0.5250 0.2980])
    plot(t,Y(3,1:n),'LineWidth',2,'Color',[0 0.4470 0.7410])
    plot(t,Y(4,1:n),'LineWidth',2,'Color',[0.2 0.6470 0.9410])
    xlim([0 steps*h])
    ylim([1.2*ymin 1.2*ymax])
    legend({'$x$','$\alpha$','$v$','$\beta$'},'Interpreter', 'latex','location','NorthEast')
    title('$x-t \,\,/\,\, \alpha-t \,\,/\,\, v-t \,\,/\,\, \beta-t$','Interpreter', 'latex','FontSize',16,'FontWeight','normal')
    xlabel('$t \, [s]$','Interpreter', 'latex')
    ylabel('$x \, [m], \,\, \alpha \, [rad], \,\, v \, [m/s], \,\, \beta \, [rad/s]$','Interpreter', 'latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
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
global k1;
global k2;
global L;
x = y(1,1);
w = y(2,1);
dx = y(3,1);
dw = y(4,1);

rechteseite_x = m2*L*dw^2*sin(w)-(k1+k2)*x;
rechteseite_w = -m2*g*L*sin(w);
ddxw = [m1+m2 m2*L*cos(w);m2*L*cos(w) m2*L^2]\[rechteseite_x;rechteseite_w];
dy= [dx;dw;ddxw(1,1);ddxw(2,1)];
end

%% Runge-Kutta
function y_neu = ruku(funk,y,h)
k1 = funk(y);
k2 = funk(y+h/2*k1);
k3 = funk(y+h/2*k2);
k4 = funk(y+h*k3);
y_neu = y+h*(k1/6+k2/3+k3/3+k4/6);
end

%% Energie
function [E] = energie(y)
global g;
global m1;
global m2;
global k1;
global k2;
global L;
x = y(1,1);
w = y(2,1);
dx = y(3,1);
dw = y(4,1);

Ep1 = 1/2*k1*x^2;
Ep2 = 1/2*k2*x^2;
Ep3 = m2*g*(L-L*cos(w));
Ek1 = 1/2*m1*dx^2;
Ek2 = 1/2*m2*((dx+L*dw*cos(w))^2+(L*dw*sin(x))^2);
E = [Ep1;Ep2;Ep3;Ek1;Ek2];
end