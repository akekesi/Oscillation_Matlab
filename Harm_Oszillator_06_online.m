%% Harmonischer Oszillator
%  Doppelpendel
%  mit Vergleich:
%  www.myphysicslab.com/pendulum/double-pendulum-en.html
%  www.youtube.com/watch?v=lSTI26zAex8

clear all
clc

global g;
global m1;
global m2;
global L1;
global L2;

vid = 0;        % als Video zu speichern: vid = 1

steps = 100;
h = 0.05;       % Zeitschritt
g = 9.81;       % Gravkonst.
m1 = 1;         % Gewicht 1
m2 = 1;         % Gewicht 2
L1 = 1;         % Laenge 1
L2 = 0.5;       % Laenge 2

w1 = 105/180*pi;    % Winkel_0 1
w2 = 180/180*pi;    % Winkel_0 2
dw1 = 0;            % Winkelgeschw_0 1
dw2 = 0;            % Winkelgeschw_0 2

y0 = [w1;w2;dw1;dw2];
Y = zeros(length(y0),steps);
Yo1 = Y;
Yo2 = Y;
E = zeros(4,1);
Eo1 = E;
Eo2 = E;

%% Berechnung
for n = 1:1:steps
    if n == 1
        y = y0;
        yo1 = y0;
        yo2 = y0;
    else
        y = ruku(@abl,y,h);
        yo1 = ruku(@abl_online1,yo1,h);
        yo2 = ruku(@abl_online2,yo2,h);
    end
    Y(:,n) = y;
    Yo1(:,n) = yo1;
    Yo2(:,n) = yo2;    
    E(:,n) = energie(y);
    Eo1(:,n) = energie(yo1);
    Eo2(:,n) = energie(yo2);
end

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_06_online.avi');
    Video.FrameRate = 13;
    open(Video)
end

for n = 1:1:steps
    x1 = L1*sin(Y(1,n));
    x2 = L1*sin(Y(1,n))+L2*sin(Y(2,n));
    y1 = -L1*cos(Y(1,n));
    y2 = -L1*cos(Y(1,n))-L2*cos(Y(2,n));

    x1o1 = L1*sin(Yo1(1,n));
    x2o1 = L1*sin(Yo1(1,n))+L2*sin(Yo1(2,n));
    y1o1 = -L1*cos(Yo1(1,n));
    y2o1 = -L1*cos(Yo1(1,n))-L2*cos(Yo1(2,n));

    x1o2 = L1*sin(Yo2(1,n));
    x2o2 = L1*sin(Yo2(1,n))+L2*sin(Yo2(2,n));
    y1o2 = -L1*cos(Yo2(1,n));
    y2o2 = -L1*cos(Yo2(1,n))-L2*cos(Yo2(2,n));

    subplot(1,2,1)
    plot(0,0,'x','MarkerSize',15,'LineWidth',2,'Color','k')
    hold on
    pruku = plot(x1,y1,'.','MarkerSize',35,'Color','k');
    plot(x2,y2,'.','MarkerSize',35,'Color','k')
    plot([0 x1],[0 y1],'Linewidth',3,'Color','k')
    plot([x1 x2],[y1 y2],'Linewidth',3,'Color','k')

    po1 = plot(x1o1,y1o1,'.','MarkerSize',25,'Color','r');
    plot(x2o1,y2o1,'.','MarkerSize',25,'Color','r')
    plot([0 x1o1],[0 y1o1],'Linewidth',2,'Color','r')
    plot([x1o1 x2o1],[y1o1 y2o1],'Linewidth',2,'Color','r')

    po2 = plot(x1o2,y1o2,'.','MarkerSize',15,'Color','g');
    plot(x2o2,y2o2,'.','MarkerSize',15,'Color','g')
    plot([0 x1o2],[0 y1o2],'Linewidth',1,'Color','g')
    plot([x1o2 x2o2],[y1o2 y2o2],'Linewidth',1,'Color','g')

    legend([pruku po1 po2],{'RK4','online-1','online-1'})

    xlim([-2.5 2.5])
    ylim([-3 3])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(1,2,2)
    plot([-10 10],[sum(E(:,1)) sum(E(:,1))],'--k')
    hold on
    p21 = plot([1 1],[0 E(1,n)],'LineWidth',15,'Color','#0072BD');
    plot([1 2],[E(1,n) E(1,n)],'--k')
    p22 = plot([2 2],[E(1,n) E(1,n)+E(2,n)],'LineWidth',15,'Color','#D95319');
    plot([2 3],[E(1,n)+E(2,n) E(1,n)+E(2,n)],'--k')
    p23 = plot([3 3],[E(1,n)+E(2,n) E(1,n)+E(2,n)+E(3,n)],'LineWidth',15,'Color','#EDB120');
    plot([3 4],[E(1,n)+E(2,n)+E(3,n) E(1,n)+E(2,n)+E(3,n)],'--k')
    p24 = plot([4 4],[E(1,n)+E(2,n)+E(3,n) E(1,n)+E(2,n)+E(3,n)+E(4,n)],'LineWidth',15,'Color','#7E2F8E');
    xlim([0 5])
    xticks([])
    ylim([1.2*(min(E(1,:))+min(E(2,:))) 1.3*sum(E(:,1))])
    legend([p21 p22 p23 p24],{'E_{Pot1}','E_{Pot2}','E_{Kin1}','E_{Kin2}'})

    text(0.5,1.17*sum(E(:,1)),['\Sigma E_0: ',num2str(sum(E(:,1)))])
    text(0.5,1.07*sum(E(:,1)),['\Sigma E_n: ',num2str(sum(E(:,n)))])
    text(0.5,0.7*sum(E(:,1)),{'\Sigma E_0-\Sigma E_n',num2str(sum(E(:,1))-sum(E(:,n))),num2str(sum(Eo1(:,1))-sum(Eo1(:,n))),num2str(sum(Eo2(:,1))-sum(Eo2(:,n)))})

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
global L1;
global L2;
w1 = y(1,1);
w2 = y(2,1);
dw1 = y(3,1);
dw2 = y(4,1);

rechteseite1 = -m2*L1*L2*dw2^2*sin(w1-w2)-(m1+m2)*g*L1*sin(w1);
rechteseite2 = m2*L1*L2*dw1^2*sin(w1-w2)-m2*g*L2*sin(w2);
ddw = [(m1+m2)*L1^2 m2*L1*L2*cos(w1-w2);m2*L1*L2*cos(w1-w2) m2*L2^2]\[rechteseite1;rechteseite2];
dy = [dw1;dw2;ddw(1,1);ddw(2,1)];
end

%% Online Version 1.
function dy = abl_online1(y)
global g;
global m1;
global m2;
global L1;
global L2;
w1 = y(1,1);
w2 = y(2,1);
dw1 = y(3,1);
dw2 = y(4,1);

% www.myphysicslab.com/pendulum/double-pendulum-en.html
ddw1 = 1/(L1*(2*m1+m2-m2*cos(2*w1-2*w2)))*(-g*(2*m1+m2)*sin(w1)-m2*g*sin(w1-2*w2)-2*sin(w1-w2)*m2*(dw2^2*L2+dw1^2*L1*cos(w1-w2)));
ddw2 = 1/(L2*(2*m1+m2-m2*cos(2*w1-2*w2)))*(2*sin(w1-w2)*(dw1^2*L1*(m1+m2)+g*(m1+m2)*cos(w1)+dw2^2*L2*m2*cos(w1-w2)));

dy = [dw1;dw2;ddw1;ddw2];
end

%% Online Version 2.
function dy = abl_online2(y)
global g;
global m1;
global m2;
global L1;
global L2;
w1 = y(1,1);
w2 = y(2,1);
dw1 = y(3,1);
dw2 = y(4,1);

% www.youtube.com/watch?v=lSTI26zAex8
M = 1+m1/m2;
ddw1 = 1/(L1*(M-cos(w1-w2)^2))*(g*(sin(w2)*cos(w1-w2)-M*sin(w1))-sin(w1-w2)*(L2*dw2^2+L1*dw1^2*cos(w1-w2)));
ddw2 = 1/(L2*(M-cos(w1-w2)^2))*(M*g*(sin(w1)*cos(w1-w2)-sin(w2))+sin(w1-w2)*(M*L1*dw1^2+L2*dw2^2*cos(w1-w2)));

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

%% Energie
function [E] = energie(y)
global g;
global m1;
global m2;
global L1;
global L2;
w1 = y(1,1);
w2 = y(2,1);
dw1 = y(3,1);
dw2 = y(4,1);

Ep1 = m1*g*(-L1*cos(w1));
Ep2 = m2*g*(-L1*cos(w1)-L2*cos(w2));
Ek1 = 1/2*m1*(L1*dw1)^2;
Ek2 = 1/2*m2*((L1*dw1*cos(w1)+L2*dw2*cos(w2))^2+(L1*dw1*sin(w1)+L2*dw2*sin(w2))^2);
E = [Ep1;Ep2;Ek1;Ek2];
end