%% Harmonischer Oszillator
%  horizontal
%  1 x Punktmasse
%  1 x Feder
%  1 x Daempfung
clear all
clc

global m;
global k;
global c;

vid = 0;    % als Video zu speichern: vid = 1

steps = 150;
h = 0.1;    % Zeitschritt
m = 1;      % Gewicht
k = 1;      % Federkoeff.

%c = 0;                  % keine Daempf. (zeta = 0)
c = 2*sqrt(m*k)*0.25;    % schwache Daempf. (0 < zeta < 1)
%c = 2*sqrt(m*k);        % krit. Daempf. (zeta = 1)
%c = 2*sqrt(m*k)*1.2;    % starke Daempf. (zeta > 1)
%c = ???;                % starke & schnellste Daempf. (zeta > 1 & v0 = x0*Lambda2)

x = 0;      % Anfangspos.
v = 5;      % Anfangsgeschw.

y0 = [x;v];
Y = zeros(length(y0),steps);

% Analytische Loesung
omega_n = sqrt(k/m);                    % csillapitatlan rendszer sajat korfrekvenciaja [rad/s]
zeta = c/(2*m*omega_n);
omega_d = omega_n*sqrt(1-zeta^2);       % csillap. rendsz. sajat koerfrekv. [rad/s]

if c == 0                               % keine Daempfung
    C1 = x;
    C2 = v/omega_n;
    A = sqrt(x^2+v^2/omega_n^2);
    epsilon = atan(omega_n*x/v);
    Tn = 2*pi/omega_n;                  % periodusido [s]
    fn = 1/Tn;                          % csillapitatlan szabad rezges sajat korfrekvenciaja [Hz = 1/s]
    x_a = @(t) C1*cos(omega_n*t)+C2*sin(omega_n*t);
    x_b = @(t) A*sin(omega_n*t+epsilon);
    Fall = 'keine Daempf.';
elseif 0 < zeta && zeta < 1             % schwache Daempfung
    C1 = x;
    C2 = (v+zeta*omega_n*x)/omega_d;
    A = sqrt(C1^2+C2^2);
    epsilon = atan(C1/C2);
    Td = 2*pi/omega_d;                   % periodusido [s]
    fd = 1/Td;                            % rezges sajat korfrekv. [Hz = 1/s]
    x_a = @(t) exp(-zeta*omega_n*t).*(C1*cos(omega_d*t)+C2*sin(omega_d*t));
    x_b = @(t) A*exp(-zeta*omega_n*t).*sin(omega_d*t+epsilon);
    x_huelle = @(t) A*exp(-zeta*omega_n*t);
    Fall = '0 < \zeta < 1';
elseif zeta == 1                        % krit. Daempfung (Aperiodischer Grenzfall)
    C1 = x;
    C2 = v+omega_n*x;
    x_a = @(t) exp(-omega_n*t).*(C1+C2*t);
    Fall = '\zeta = 1';
elseif 1 < zeta                         % starke Daempfung
    Lambda1 = -zeta*omega_n+omega_n*sqrt(zeta^2-1);
    Lambda2 = -zeta*omega_n-omega_n*sqrt(zeta^2-1);
    C1 = (v-x*Lambda2)/(Lambda1-Lambda2);
    C2 = (v-x*Lambda1)/(Lambda2-Lambda1);
    x_a = @(t) C1*exp(Lambda1*t)+C2*exp(Lambda2*t);
    if v-Lambda2*x < 10^(-12)
        Fall = '1 < \zeta & v_0 = x_0\cdot\lambda_2 \rightarrow schnellste';
    else
        Fall = '1 < \zeta & v_0 \neq x_0\cdot\lambda_2 \rightarrow nicht schnellste';
    end
end

% Numerische Loesung
for n = 1:1:steps
    if n == 1
        y = y0;
    else
        y = ruku(@abl,y,h);
    end
    Y(:,n) = y;
end

%% Plot
% Hilfswerte fuer Plot
xmax = max(Y(1,:));
xmin = min(Y(1,:));
xMax = max(xmax,abs(xmin));
L0 = fix(15*xMax)/10;
if L0 < 1
    L0 = 1;
end

vmax = max(Y(2,:));
vmin = min(Y(2,:));
vMax = max(vmax,abs(vmin));

ymax = fix(15*max(xmax,vmax))/10;
ymin = fix(15*min(xmin,vmin))/10;

Max = max(ymax,abs(ymin));

LW0 = 3;    % Linewidth
LWmin = 1.5;
LWd = (LW0-LWmin)/(xMax);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_01.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps

    sgtitle({'Harmonischer Oszillator','Horizontal','Masse:    1x','Feder:     1x','Daempf.: 1x'},'FontSize',22)

    subplot(1,3,1)
    geschw_dx = 1.0;
    geschw_x = geschw_dx*Y(2,n)/vMax;
    geschw_y = -0.9;
    if geschw_x > 0
        mark_geschw = '>';
    elseif geschw_x < 0
        mark_geschw = '<';
    else
        mark_geschw = '.';
    end

    kraft_f_dx = 1.0;
    kraft_f_x = kraft_f_dx*(-1)*Y(1,n)*k/(xMax*k);
    kraft_f_y = geschw_y-0.5;
    if kraft_f_x > 0
        mark_f_kraft = '>';
    elseif kraft_f_x < 0
        mark_f_kraft = '<';
    else
        mark_f_kraft = '.';
    end

    kraft_d_dx = 1.0;
    kraft_d_x = kraft_d_dx*(-1)*Y(2,n)*c/(vMax*c);
    kraft_d_y = kraft_f_y-0.5;
    if kraft_d_x > 0
        mark_d_kraft = '>';
    elseif kraft_d_x < 0
        mark_d_kraft = '<';
    else
        mark_d_kraft = '.';
    end

    p11 = plot(Y(1,1),0,'o','MarkerSize',10,'LineWidth',2,'Color','k');                         % Anf.
    hold on
    p13 = plot([-L0 Y(1,n)],[0 0],'LineWidth',LW0-(Y(1,n))*LWd,'Color','#77AC30');              % Feder
    plot([-L0 -L0],[-Max Max],'LineWidth',5,'Color','k')                                        % Wand
    plot([Y(1,n) Y(1,n)],[0 kraft_d_y],'--','Color',[0.7 0.7 0.7])
    plot(Y(1,n)-geschw_dx,geschw_y,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+geschw_dx,geschw_y,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-geschw_dx Y(1,n)+geschw_dx],[geschw_y geschw_y],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p14 = plot([Y(1,n) Y(1,n)+geschw_x],[geschw_y geschw_y],'LineWidth',3,'Color','#D95319');   % Geschw.
    plot(Y(1,n)+geschw_x,geschw_y,'Marker',mark_geschw,'LineWidth',3,'Color','#D95319');        % Geschw.
    plot(Y(1,n)-kraft_f_dx,kraft_f_y,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+kraft_f_dx,kraft_f_y,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-kraft_f_dx Y(1,n)+kraft_f_dx],[kraft_f_y kraft_f_y],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p15 = plot([Y(1,n) Y(1,n)+kraft_f_x],[kraft_f_y kraft_f_y],'LineWidth',3,'Color','#7E2F8E');      % Kraft_Feder
    plot(Y(1,n)+kraft_f_x,kraft_f_y,'Marker',mark_f_kraft,'LineWidth',3,'Color','#7E2F8E');           % Kraft_Feder
    plot(Y(1,n)-kraft_d_dx,kraft_d_y,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+kraft_d_dx,kraft_d_y,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-kraft_d_dx Y(1,n)+kraft_d_dx],[kraft_d_y kraft_d_y],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p16 = plot([Y(1,n) Y(1,n)+kraft_d_x],[kraft_d_y kraft_d_y],'LineWidth',3,'Color','#A2142F');      % Kraft_Daempf.
    plot(Y(1,n)+kraft_d_x,kraft_d_y,'Marker',mark_d_kraft,'LineWidth',3,'Color','#A2142F');           % Kraft_Daempf.
    p12 = plot(Y(1,n),0,'.','MarkerSize',30,'Color','k');                                              % Gewicht
    xlim([-L0 L0])
    xticks(-fix(L0):1:fix(L0))
    ylim([-3 4])
    yticks([])
    legend([p11 p12 p13 p14 p15 p16],{'Anfang','Masse','Feder','Geschw.','Kraft-F','Kraft-D'},'location','NorthEast')
    title('Animation','FontSize',16,'FontWeight','normal')
    text(-L0*0.9,3.0,{['$ m = $',num2str(m)','$ kg $'],['$ \, k \, = $',num2str(k),'$ \frac{N}{m} $'],['$ \, c \, = $',num2str(c),'$ \frac{Ns}{m} $']},'Interpreter', 'latex')
    text(-L0*0.9+1.75,3.0,{['$ x_0 = $',num2str(Y(1,1)),'$ m $'],['$ v_0 = $',num2str(Y(2,1)),'$ \frac{m}{s} $']},'Interpreter', 'latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(1,3,2)
    t=0:h:h*(n-1);
    p21 = plot(t,Y(1,1:n),'Color','#0072BD','LineWidth',2);
    hold on
    p22 = plot(t,Y(2,1:n),'Color','#D95319','LineWidth',2);
    xlim([0 h*steps])
    ylim([-Max Max])
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
    if c == 0
        p23 = plot(t,x_a(t),'x','MarkerSize',5,'Color','k');
        p24 = plot(t,x_b(t),'o','MarkerSize',5,'Color','k');
        legend([p21 p22 p23 p24],{'x (RK-4)','v (RK-4)','x (ana. a)','x (ana. b)'},'location','NorthEast')
        text(h*steps/10,0.8*Max,['$T_n = \frac{2\pi}{\omega_n} = 2\pi\sqrt{\frac{m}{k}} =$',num2str((double(int64(Tn*100)))/100),'$s$'],'Interpreter', 'latex')
        title(Fall,'FontSize',16,'FontWeight','normal')
    elseif 0 < zeta && zeta < 1
        p23 = plot(t,x_a(t),'x','MarkerSize',5,'Color','k');
        p24 = plot(t,x_b(t),'o','MarkerSize',5,'Color','k');
        p25 = plot(t,x_huelle(t),'k');
        plot(t,-x_huelle(t),'k')
        legend([p21 p22 p23 p24 p25],{'x (RK-4)','v (RK-4)','x (ana. a)','x (ana. b)','Huelle'},'location','NorthEast')
        text(h*steps/10,0.8*Max,['$T_d = \frac{2\pi}{\omega_d} =$',num2str((double(int64(Td*100)))/100),'$s$'],'Interpreter', 'latex')
        title(Fall,'FontSize',16,'FontWeight','normal')
    elseif zeta == 1
        p23 = plot(t,x_a(t),'x','MarkerSize',5,'Color','k');
        legend([p21 p22 p23],{'x (RK-4)','v (RK-4)','x (ana. a)'},'location','NorthEast')
        title(Fall,'FontSize',16,'FontWeight','normal')
    elseif 1 < zeta
        p23 = plot(t,x_a(t),'x','MarkerSize',5,'Color','k');
        legend([p21 p22 p23],{'x (RK-4)','v (RK-4)','x (ana. a)'},'location','NorthEast')
        title(Fall,'FontSize',16,'FontWeight','normal')
    end
    daspect([1 1 1])
    drawnow
    hold off

    subplot(1,3,3)
    p31 = plot(Y(1,1),Y(2,1),'o','MarkerSize',10,'LineWidth',2,'Color','k');
    hold on
    p32 = plot(Y(1,1:n),Y(2,1:n),'k','LineWidth',2);
    plot(Y(1,n),Y(2,n),'.','Color','k','Markersize',30)
    xlim([-Max Max])
    ylim([-Max Max])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    legend([p31 p32],{'Anfang','Phasenkurve'},'location','NorthEast')
    title('Phasenportrait','FontSize',16,'FontWeight','normal')
    xlabel('x')
    ylabel('v')
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
function dy = abl(y)
global m;
global k;
global c;
x = y(1,1);
v = y(2,1);
a = -k/m*x-c/m*v;
dy = [v;a];
end

%% Runge-Kutta
function y_neu = ruku(funk,y,h)
k1 = funk(y);
k2 = funk(y+h/2*k1);
k3 = funk(y+h/2*k2);
k4 = funk(y+h*k3);
y_neu = y+h*(k1/6+k2/3+k3/3+k4/6);
end