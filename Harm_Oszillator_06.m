%% Harmonischer Oszillator
%  Doppelpendel
clear all
clc

global g;
global m1;
global m2;
global L1;
global L2;

vid = 0;        % als Video zu speichern: vid = 1

steps = 200;
h = 0.05;       % Zeitschritt
g = 9.81;       % Gravkonst.
m1 = 1;         % Gewicht 1
m2 = 1;         % Gewicht 2
L1 = 1;         % Laenge 1
L2 = 1;         % Laenge 2

w1 = 105/180*pi;    % Winkel_0 1
w2 = 180/180*pi;    % Winkel_0 2
dw1 = 0;            % Winkelgeschw_0 1
dw2 = 0;            % Winkelgeschw_0 2

y0 = [w1;w2;dw1;dw2];
Y = zeros(length(y0),steps);
E = zeros(4,1);

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
EMax = E(1,1)+E(2,1)+E(3,1)+E(4,1);
w1max = max(Y(1,:));
w1min = min(Y(1,:));

w2max = max(Y(2,:));
w2min = min(Y(2,:));

wMax = max(w1max,w2max);
wMin = min(w1min,w2min);

dw1max = max(Y(3,:));
dw1min = min(Y(3,:));
%dw1Max = max(dw1max,abs(dw1min));
dw2max = max(Y(4,:));
dw2min = min(Y(4,:));
%dw2Max = max(dw2max,abs(dw2min));
dwMax = max(dw1max,dw2max);
dwMin = min(dw1min,dw2min);

Max = max(wMax,dwMax);
Min = min(wMin,dwMin);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_06.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps
    x1 = L1*sin(Y(1,n));
    x2 = L1*sin(Y(1,n))+L2*sin(Y(2,n));
    y1 = -L1*cos(Y(1,n));
    y2 = -L1*cos(Y(1,n))-L2*cos(Y(2,n));

    subplot(2,3,1)
    plot(0,0,'x','MarkerSize',15,'LineWidth',2,'Color','k')
    hold on
    plot([0 x1],[0 y1],'Color','#0072BD')
    plot([x1 x2],[y1 y2],'Color','#D95319')

    if n > 6
        for m = 0:1:5
            C = m/20;
            plot(L1*sin(Y(1,n-1-m:n-m)),-L1*cos(Y(1,n-1-m:n-m)),'LineWidth',2,'Color',[0+C 0.4470+C 0.7410+C])
            plot(L1*sin(Y(1,n-1-m:n-m))+L2*sin(Y(2,n-1-m:n-m)),-L1*cos(Y(1,n-1-m:n-m))-L2*cos(Y(2,n-1-m:n-m)),'LineWidth',2,'Color',[0.8500+C/10 0.3250+C 0.0980+C])
        end
    elseif n > 2 && n <=6
        for m = 0:1:n-2
            C = m/20;    
            plot(L1*sin(Y(1,n-1-m:n-m)),-L1*cos(Y(1,n-1-m:n-m)),'LineWidth',2,'Color',[0+C 0.4470+C 0.7410+C])
            plot(L1*sin(Y(1,n-1-m:n-m))+L2*sin(Y(2,n-1-m:n-m)),-L1*cos(Y(1,n-1-m:n-m))-L2*cos(Y(2,n-1-m:n-m)),'LineWidth',2,'Color',[0.8500+C/10 0.3250+C 0.0980+C])
        end
    end

    p11 = plot(x1,y1,'.','MarkerSize',35,'Color','#0072BD');
    p12 = plot(x2,y2,'.','MarkerSize',35,'Color','#D95319');

    xlim([-(L1+L2)*1.2 (L1+L2)*1.2])
    ylim([-(L1+L2)*1.2 (L1+L2)*1.2])
    legend([p11 p12],{'Masse-1','Masse-2'},'location','SouthEast')
    title('Animation','FontSize',15,'FontWeight','normal')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,3,2)
    p21 = plot([0 0],[0 E(2,n)],'LineWidth',150,'Color',[0.8500 0.3250 0.0980],'LineWidth',15);
    hold on
    p22 = plot([0 0],[E(2,n) E(2,n)+E(4,n)],'LineWidth',150,'Color',[1.0 0.5250 0.2980],'LineWidth',15);
    p23 = plot([0 0],[E(2,n)+E(4,n) E(2,n)+E(4,n)+E(1,n)],'LineWidth',150,'Color',[0 0.4470 0.7410],'LineWidth',15);
    p24 = plot([0 0],[E(2,n)+E(4,n)+E(1,n) E(2,n)+E(4,n)+E(1,n)+E(3,n)],'LineWidth',150,'Color',[0.2 0.6470 0.9410],'LineWidth',15);
    plot([-EMax EMax],[EMax EMax],'--k')
    xlim([-1 1])
    xticks([])
    ylim([0 fix(EMax*1.2)])
    legend([p24 p23 p22 p21],{'E_{Kin1}','E_{Pot1}','E_{Kin2}','E_{Pot2}'},'Location','SouthEast')
    text(-0.7*EMax,EMax*1.05,'\Sigma E')
    title('Energie','FontSize',16,'FontWeight','normal')
    ylabel('E [J]')
    grid on
    grid minor
%    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,3,[3 6])
    plot(Y(1,n),Y(2,n),'.','MarkerSize',15,'Color','k')
    hold on
    plot(Y(1,1:n),Y(2,1:n),'k')
    xlim([w1min*1.2 w1max*1.2])
    ylim([w2min*1.2 w2max*1.2])
    title('Phasenportait','FontSize',16,'FontWeight','normal')
    xlabel('\alpha_1')
    ylabel('\alpha_2')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(2,3,[4 5])
    t=0:h:h*(n-1);
    p41 = plot(t,Y(1,1:n),'Color','#0072BD');
    hold on
    p42 = plot(t,Y(3,1:n),'--','Color','#0072BD');
    p43 = plot(t,Y(2,1:n),'Color','#D95319');
    p44 = plot(t,Y(4,1:n),'--','Color','#D95319');
    xlim([0 steps*h])
    ylim([Min*1.2 Max*1.2])
    legend([p41 p42 p43 p44],{'\alpha_1','\beta_1','\alpha_2','\beta_2'},'location','NorthEast')
    title('$\alpha-t \,\,/\,\, \beta-t$','Interpreter', 'latex','FontSize',16,'FontWeight','normal')
    xlabel('$t \, [s]$','Interpreter', 'latex')
    ylabel('$\alpha \, [rad], \,\, \beta \, [rad/s]$','Interpreter', 'latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
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
w1 = y(1,1);
w2 = y(2,1);
dw1 = y(3,1);
dw2 = y(4,1);

rechteseite1 = -m2*L1*L2*dw2^2*sin(w1-w2)-(m1+m2)*g*L1*sin(w1);
rechteseite2 = m2*L1*L2*dw1^2*sin(w1-w2)-m2*g*L2*sin(w2);
ddw = [(m1+m2)*L1^2 m2*L1*L2*cos(w1-w2);m2*L1*L2*cos(w1-w2) m2*L2^2]\[rechteseite1;rechteseite2];
dy = [dw1;dw2;ddw(1,1);ddw(2,1)];
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

Ep1 = m1*g*(L1+L2-L1*cos(w1));
Ep2 = m2*g*(L1+L2-L1*cos(w1)-L2*cos(w2));
Ek1 = 1/2*m1*(L1*dw1)^2;
Ek2 = 1/2*m2*((L1*dw1*cos(w1)+L2*dw2*cos(w2))^2+(L1*dw1*sin(w1)+L2*dw2*sin(w2))^2);
E = [Ep1;Ep2;Ek1;Ek2];
end