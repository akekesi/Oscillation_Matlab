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

c1 = 0;                 % keine Daempf. (zeta = 0)
c2 = 2*sqrt(m*k)*0.25;  % schwache Daempf. (0 < zeta < 1)
c3 = 2*sqrt(m*k);       % krit. Daempf. (zeta = 1)
c4 = 2*sqrt(m*k)*4;     % starke Daempf. (zeta > 1)
%c5 = 0;                 % ??? starke & schnellste Daempf. (zeta > 1 & v0 = x0*Lambda2) ???
C = [c1;c2;c3;c4];
Lc = length(C);

x = 0;      % Anfangspos.
v = 5;      % Anfangsgeschw.
y0 = [x;v];

%% Loesungen
%  numerische Loesung
t = 0:h:h*(steps-1);
Y = zeros(length(y0),steps,Lc);
Fall = cell(Lc,1);
for L = 1:1:Lc
    c = C(L,1);
    Y(:,:,L) = num_lsg(y0,steps,h);
    [~,~,Fall_,~] = anafunc(m,k,c,x,v);
    Fall{L,1} = Fall_{2,1};
end

%% Plot
% Hilfswerte fuer Plot
maxvx = max(Y,[],[2 3]);
minvx = min(Y,[],[2 3]);
xmax = maxvx(1,1);
xmin = minvx(1,1);
xMax = max(xmax,abs(xmin));
L0 = fix(15*xMax)/10;
if L0 < 1
    L0 = 1;
end

vmax = maxvx(2,1);
vmin = minvx(2,1);
vMax = max(vmax,abs(vmin));

ymax = fix(15*max(xmax,vmax))/10;
ymin = fix(15*min(xmin,vmin))/10;

Max = max(ymax,abs(ymin));

LW0 = 3;    % Linewidth
LWmin = 1.5;
LWd = (LW0-LWmin)/(xMax);

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_02.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps

    for L = 1:1:Lc

        subplot(Lc,3,3*L-2)

        p11 = plot(Y(1,1,L),0,'o','MarkerSize',10,'LineWidth',2,'Color','k');               % Anf.
        hold on
        p13 = plot([-L0 Y(1,n,L)],[0 0],'LineWidth',LW0-(Y(1,n,L))*LWd,'Color','#77AC30');  % Feder
        plot([-L0 -L0],[-Max Max],'LineWidth',5,'Color','k')                                % Wand
        p12 = plot(Y(1,n,L),0,'.','MarkerSize',30,'Color','k');                             % Gewicht
        xlim([-L0 L0])
        xticks(-fix(L0):1:fix(L0))
        ylim([-1 1])
        yticks([])
        title(Fall{L,1},'FontSize',16,'FontWeight','normal')
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        grid on
        grid minor
        daspect([1 1 1])
        drawnow
        hold off

        subplot(Lc,3,3*L-1)
        p21 = plot(t(1:n),Y(1,1:n,L),'Color','#0072BD','LineWidth',1);
        hold on
        p22 = plot(t(1:n),Y(2,1:n,L),'Color','#D95319','LineWidth',1);
        xlim([0 h*steps])
        ylim([-Max Max])
        ax = gca;
        ax.XAxisLocation = 'origin';
        grid on
        grid minor
        daspect([1 1 1])
        drawnow
        hold off

        subplot(Lc,3,3*L)
        p31 = plot(Y(1,1,L),Y(2,1,L),'o','MarkerSize',5,'LineWidth',1,'Color','k');
        hold on
        p32 = plot(Y(1,1:n,L),Y(2,1:n,L),'k','LineWidth',1);
        plot(Y(1,n,L),Y(2,n,L),'.','Color','k','Markersize',15)
        xlim([-Max Max])
        ylim([-Max Max])
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
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
end

if vid == 1
    close(Video)
end

%% Numerische Loesung
function [Y] = num_lsg(y0,steps,h) 
Y = zeros(length(y0),steps);
for n = 1:1:steps
    if n == 1
        y = y0;
    else
        y = ruku(@abl,y,h);
    end
    Y(:,n) = y;
end
end

%% Abl.
function [dy] = abl(y)
global m;
global k;
global c;
x = y(1,1);
v = y(2,1);
a = -k/m*x-c/m*v;
dy = [v;a];
end

%% Runge-Kutta
function [y_neu] = ruku(funk,y,h)
k1 = funk(y);
k2 = funk(y+h/2*k1);
k3 = funk(y+h/2*k2);
k4 = funk(y+h*k3);
y_neu = y+h*(k1/6+k2/3+k3/3+k4/6);
end

%% Analytische Loesung
function [returnfunc,returntext,Fall,T] = anafunc(m,k,c,x,v)
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
    returnfunc = cell(2,1);
    returntext = cell(2,1);
    Fall = cell(2,1);
    returnfunc{1,1} = @(t) C1*cos(omega_n*t)+C2*sin(omega_n*t);
    returntext{1,1} = 'x (ana. a)';
    returnfunc{2,1} = @(t) A*sin(omega_n*t+epsilon);
    returntext{2,1} = 'x (ana. b)';
    Fall{1,1} = 1;
    Fall{2,1} = 'keine Daempf.';
    T = Tn;
elseif 0 < zeta && zeta < 1             % schwache Daempfung
    C1 = x;
    C2 = (v+zeta*omega_n*x)/omega_d;
    A = sqrt(C1^2+C2^2);
    epsilon = atan(C1/C2);
    Td = 2*pi/omega_d;                   % periodusido [s]
    fd = 1/Td;                           % rezges sajat korfrekv. [Hz = 1/s]
    returnfunc = cell(3,1);
    returntext = cell(3,1);
    Fall = cell(2,1);
    returnfunc{1,1}= @(t) exp(-zeta*omega_n*t).*(C1*cos(omega_d*t)+C2*sin(omega_d*t));
    returntext{1,1} = 'x (ana. a)';
    returnfunc{2,1} = @(t) A*exp(-zeta*omega_n*t).*sin(omega_d*t+epsilon);
    returntext{2,1} = 'x (ana. b)';
    returnfunc{3,1} = @(t) A*exp(-zeta*omega_n*t);
    returntext{3,1} = 'Huelle';
    Fall{1,1} = 2;
    Fall{2,1} = '0 < \zeta < 1';
    T = Td;
elseif zeta == 1                        % krit. Daempfung (Aperiodischer Grenzfall)
    C1 = x;
    C2 = v+omega_n*x;
    returnfunc = cell(1,1);
    returntext = cell(1,1);
    Fall = cell(2,1);
    returnfunc{1,1} = @(t) exp(-omega_n*t).*(C1+C2*t);
    returntext{1,1} = 'x (ana. a)';
    Fall{1,1} = 3;
    Fall{2,1} = '\zeta = 1';
    T = 0;
elseif 1 < zeta                         % starke Daempfung
    Lambda1 = -zeta*omega_n+omega_n*sqrt(zeta^2-1);
    Lambda2 = -zeta*omega_n-omega_n*sqrt(zeta^2-1);
    C1 = (v-x*Lambda2)/(Lambda1-Lambda2);
    C2 = (v-x*Lambda1)/(Lambda2-Lambda1);
    returnfunc = cell(1,1);
    returntext = cell(1,1);
    Fall = cell(2,1);
    returnfunc{1,1} = @(t) C1*exp(Lambda1*t)+C2*exp(Lambda2*t);
    returntext{1,1} = 'x (ana. a)';
    Fall{1,1} = 4;
    if v-Lambda2*x < 10^(-12)
        Fall{2,1} = '1 < \zeta & v_0 = x_0\cdot\lambda_2 \rightarrow schnellste';
    else
        Fall{2,1} = '1 < \zeta & v_0 \neq x_0\cdot\lambda_2 \rightarrow nicht schnellste';
    end
    T = 0;
end
end