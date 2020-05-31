%% Tuerme von Hanoi
clear all
clc

%% Anzahl der Scheiben
N = 5;

%% Anfangsaufbau
x0 = zeros(N,3);
for n = 1:1:N
    x0(n,1) = n;
end
x = x0;
X = x;

%% Ergebnismatrix E
E = strings(N+1,1);
for n = 1:1:N
    E(n+1,1) = num2str(n);
end

%% Hilfmatrix A
% 1.Reihe: Daten von Scheibe-1
% 2.Reihe: Daten von Scheibe-1
% ...

% Spalte 1-2-3: Reihenfolge der Position von Scheibe (Abhaengigkeit: N-gerade/ungera)
% 1.Spalt: Position von Scheibe nach 1. Zug
% 2.Spalt: Position von Scheibe nach 2. Zug
% 3.Spalt: Position von Scheibe nach 3. Zug
% wiederholen
% 4.Spalt: Schritt wann 1.Zug dieser Scheibe durchgeführt wird
% 5.Spalt: Abstand zwischen 2 nachenander kommenden Züge
% 6.Spalt: Schritte --> Bestimmt die nächste Position (Spalte 1-2-3)
A = zeros(1,6);

for n = 1:1:N
    if mod(N,2) == 0
        if mod(n,2) == 1
            r = [2 3 1];
        else
            r = [3 2 1];
        end
    else
        if mod(n,2) == 1
            r = [3 2 1];
        else
            r = [2 3 1];
        end
    end
    A(n,1:3) = r;
    A(n,4) = 2^(n-1);
    A(n,5) = 2^n;
    A(n,6) = 1;
end

%% Hauptcode
s = 1;
while x(1,3) ~= 1
    for n = N:-1:1
        if mod(s-2^(n-1),2^n) == 0 % Welche Scheibe wird in diesem Schritt umgelegt?
            [m_alt,n_alt] = find(x == n); % Wo ist diese Scheibe?
            [m_neu,n_neu] = wofrei(x,A(n,A(n,6))); % Wo ist die neue Pos. dieser Scheibe?
            [x(m_alt,n_alt),x(m_neu,n_neu)] = swap(x(m_alt,n_alt),x(m_neu,n_neu)); % Durchführung des Zuges.
            E(1,s+1) = num2str(s); % Züge wurde im ersten Zeilen abgespeichert.
            E(n+1,s+1) = num2str(A(n,A(n,6))); % In welchem Zug, welche Scheibe wird auf welchem Stab umgelegt.
            A(n,6) = schritt(A(n,6),3); % Die naechste Pos. abspeichern.
            s = s+1; % Inkrementieren des Schrittes.
            X(:,:,s) = x;
        end
    end
end

disp(X)
disp(['Steps: ',num2str(s-1)])
%xlswrite('Tuerme_von_Hanoi_00.xls',E)

%% Plot
fig = figure('WindowState','maximized');
Video = VideoWriter('Tuerme_von_Hanoi_00.avi');
Video.FrameRate = 1; 
open(Video)

for m = 1:1:s
    for n = N:-1:1
        multi = 30;
        stab = 5;
        col_stab = 'k'; 
        txt = '';
        if X(n,1,m) == 0
            LW1 = stab;
            COL1 = col_stab;
            txt1 = txt;
        else
            LW1 = multi*X(n,1,m);
            COL1 = [1-X(n,1,m)/10 X(n,1,m)/10 X(n,1,m)/10];
            txt1 = num2str(X(n,1,m));
        end
        if X(n,2,m) == 0
            LW2 = stab;
            COL2 = col_stab;
            txt2 = txt;
        else
            LW2 = multi*X(n,2,m);
            COL2 = [1-X(n,2,m)/10 X(n,2,m)/10 X(n,2,m)/10];
            txt2 = num2str(X(n,2,m));
        end
        if X(n,3,m) == 0
            LW3 = stab;
            COL3 = col_stab;
            txt3 = txt;
        else
            LW3 = multi*X(n,3,m);
            COL3 = [1-X(n,3,m)/10 X(n,3,m)/10 X(n,3,m)/10];
            txt3 = num2str(X(n,3,m));
        end

        plot([1 1],[N-n N-n+1],'Color',COL1,'LineWidth',LW1)
        hold on
        plot([2 2],[N-n N-n+1],'Color',COL2,'LineWidth',LW2)
        plot([3 3],[N-n N-n+1],'Color',COL3,'LineWidth',LW3)

        text(1-0.01,N-n+0.5,txt1,'FontSize',15,'FontWeight','bold')
        text(2-0.01,N-n+0.5,txt2,'FontSize',15,'FontWeight','bold')
        text(3-0.01,N-n+0.5,txt3,'FontSize',15,'FontWeight','bold')

        text(3.5,N-0.5,['Steps: ',num2str(m-1)])
        xlim([0 4])
        xticks(1:1:3)
        ylim([0 N])
        yticks(0:1:N)
        grid on
        daspect([1 N 1])
    end
    drawnow
    hold off
    frame = getframe(fig);
    writeVideo(Video,frame);
end
close(Video)

%% swap
function [a,b] = swap(a,b)
tmp = a;
a = b;
b = tmp;
end

%% Wo frei?
function [m_neu,n_neu] = wofrei(x,spalt)
[m,~] = size(x);
n_neu = spalt;
for n = m:-1:1
    if x(n,spalt) == 0
        m_neu = n;
        break
    end
end
end

%% Schritt
function [neu] = schritt(zahl,max)
if zahl == max
    neu = 1;
else
    neu = zahl+1;
end
end