%% Harmonischer Oszillator
%  hizontal
%  1 x Punktmasse
%  1 x Feder
clear all
clc

global m;
global k;

vid = 0;    % als Video zu speichern: vid = 1

steps = 230;
h = 0.1;    % Zeitschritt
m = 1;      % Gewicht
k = 1;      % Federkoeff.

x = 1;      % Anfangspos.
v = 0;      % Anfangsgeschw.

y0 = [x;v];
Y = zeros(length(y0),steps);

% Analytische Loesung
omega_n = sqrt(k/m);    % csillapitatlan rendszer sajat korfrekvenciaja [rad/s]
C1 = x;
C2 = v/omega_n;
A = sqrt(x^2+v^2/omega_n^2);
epsilon = atan(omega_n*x/v);
Tn = 2*pi/omega_n;      % Periodusido [s]
fn = 1/Tn;              % csillapitatlan szabad rezges sajat korfrekvenciaja [Hz = 1/s]

x_a = @(t) C1*cos(omega_n*t)+C2*sin(omega_n*t);
x_b = @(t) A*sin(omega_n*t+epsilon);

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

nT = zeros(1,fix(h*steps/Tn));  % Periodenzeiten
snT = 0;

% Video
if vid == 1
    Video = VideoWriter('Harm_Oszillator_00_inf_b.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for n = 1:1:steps

    sgtitle({'Harmonische Oszillator','Horizontal','Masse: 1x','Feder:  1x'},'FontSize',22)

    subplot(1,3,1)
    geschw_dx = 0.3;
    geschw_x = geschw_dx*Y(2,n)/vMax;
    geschw_y = -0.4;
    if geschw_x > 0
        mark_geschw = '>';
    elseif geschw_x < 0
        mark_geschw = '<';
    else
        mark_geschw = '.';
    end

    kraft_dx = 0.3;
    kraft_x = kraft_dx*(-1)*Y(1,n)*k/(xMax*k);
    kraft_y = geschw_y-0.2;
    if kraft_x > 0
        mark_kraft = '>';
    elseif kraft_x < 0
        mark_kraft = '<';
    else
        mark_kraft = '.';
    end

    p11 = plot(Y(1,1),0,'o','MarkerSize',10,'LineWidth',2,'Color','k');                         % Anf.
    hold on
    p13 = plot([-L0 Y(1,n)],[0 0],'LineWidth',LW0-(Y(1,n))*LWd,'Color','#77AC30');              % Feder
    plot([-L0 -L0],[1 -1],'LineWidth',5,'Color','k')                                            % Wand
    plot([Y(1,n) Y(1,n)],[0 kraft_y],'--','Color',[0.7 0.7 0.7])
    plot(Y(1,n)-geschw_dx,geschw_y,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+geschw_dx,geschw_y,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-geschw_dx Y(1,n)+geschw_dx],[geschw_y geschw_y],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p14 = plot([Y(1,n) Y(1,n)+geschw_x],[geschw_y geschw_y],'LineWidth',3,'Color','#D95319');   % Geschw.
    plot(Y(1,n)+geschw_x,geschw_y,'Marker',mark_geschw,'LineWidth',3,'Color','#D95319');        % Geschw.
    plot(Y(1,n)-kraft_dx,kraft_y,'Marker','<','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot(Y(1,n)+kraft_dx,kraft_y,'Marker','>','LineWidth',3,'Color',[0.7 0.7 0.7])
    plot([Y(1,n)-kraft_dx Y(1,n)+kraft_dx],[kraft_y kraft_y],'LineWidth',3,'Color',[0.7 0.7 0.7])
    p15 = plot([Y(1,n) Y(1,n)+kraft_x],[kraft_y kraft_y],'LineWidth',3,'Color','#7E2F8E');      % Kraft
    plot(Y(1,n)+kraft_x,kraft_y,'Marker',mark_kraft,'LineWidth',3,'Color','#7E2F8E');           % Kraft
    p12 = plot(Y(1,n),0,'.','MarkerSize',30,'Color','k');                                       % Gewicht
    xlim([-L0 L0])
    ylim([-1 1])
    yticks(-1:1:1)
    legend([p11 p12 p13 p14 p15],{'Anfang','Masse','Feder','Geschw.','Kraft'},'location','NorthEast')
    title('Animation','FontSize',16,'FontWeight','normal')
    text(-L0*0.9,0.7,{['$ x_0 = $',num2str(Y(1,1)),'$ m $'],['$ v_0 = $',num2str(Y(2,1)),'$ \frac{m}{s} $']},'Interpreter', 'latex')
    text(-L0*0.9+0.6,0.7,{['$ m = $',num2str(m)','$ kg $'],['$ \,\,k = $',num2str(k),'$ \frac{N}{m} $']},'Interpreter', 'latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(1,3,2)
    x2lim = Tn*1.5;
    t = 0:h:h*(n-1);
    if t(n) < x2lim
        x2e = x2lim;
        x2a = x2e - x2lim;
    else
        x2e = t(n);
        x2a = x2e - x2lim;
    end
    p23 = plot(t,x_a(t),'x','MarkerSize',5,'Color','k');              % x_a
    hold on
    p24 = plot(t,x_b(t),'o','MarkerSize',5,'Color','k');              % x_b
    p21 = plot(t,Y(1,1:n),'Color','#0072BD','LineWidth',2);           % x
    p22 = plot(t,Y(2,1:n),'Color','#D95319','LineWidth',2);           % v
    text(x2a+0.3,1.6*ymax,['$T_n = \frac{2\pi}{\omega_n} = 2\pi\sqrt{\frac{m}{k}} = $',num2str((double(int64(Tn*100)))/100),'$ s$'],'Interpreter', 'latex')
    if t(n) > (1+snT)*Tn && t(n)-h < (1+snT)*Tn
        snT = snT + 1;
        nT(1,snT) = n;
    end
    for nn = 1:1:snT
        plot([t(nT(nn)) t(nT(nn))],[ymax 0.9*ymin],':','Color','k','LineWidth',1.5)   % Tn
        if t(nT(nn)) > x2a+0.5 && t(nT(nn)) < x2e-0.5
            text(t(nT(nn))+0.1,1.1*ymin,[num2str(nn),'xT_n'],'HorizontalAlignment','center')
        end
    end

    xlim([x2a x2e])
    ylim([2*ymin 2*ymax])
    yticks(int64(2*ymin-1):1:int64(2*ymax-1))
    title('x-t / v-t','FontSize',16,'FontWeight','normal')
    legend([p21 p22 p23 p24],{'x (RK-4)','v (RK-4)','x (ana. a)','x (ana. b)'},'location','NorthEast')
    xlabel('t')
    ylabel('x / v')
    ax = gca;
    ax.XAxisLocation = 'origin';
    grid on
    grid minor
    daspect([1 1 1])
    drawnow
    hold off

    subplot(1,3,3)
    p30 = plot(Y(1,1),Y(2,1),'o','MarkerSize',10,'LineWidth',2,'Color','k');
    hold on
    p31 = plot(Y(1,1:n),Y(2,1:n),'k','LineWidth',2);
    p32 = plot([0 Y(1,n)],[0 0],'Color','#0072BD','LineWidth',2);
    p34 = plot([0 0],[0 Y(2,n)],'Color','#D95319','LineWidth',2);
    plot([Y(1,n) Y(1,n)],[0 Y(2,n)],'k')
    plot([0 Y(1,n)],[Y(2,n) Y(2,n)],'k')
    plot([0 Y(1,n)],[0 Y(2,n)],'k')
    plot(Y(1,n),Y(2,n),'.','Color','k','Markersize',30)
    xlim([-Max Max])
    ylim([-Max Max])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    legend([p30 p31 p32 p34],{'Anfang','Phasenkurve','Position','Geschw.'},'location','NorthEast')
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
x = y(1,1);
v = y(2,1);
a = -k/m*x;
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