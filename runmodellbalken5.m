close all
global n homvars homsols E Q_0 rho A L I Q_1 Q_2
%E...E-Modul in N/m2
E=0.01;
%rho...Dichte
rho=0.01;
%A... Querschnitt
A=0.01;
%L... Länge der Welle
L=0.01;
%I... Flächenträgheitsmoment (Biegung)
I=0.01;
%Q_0... konstante Kraft
Q_0=0.01;
%Q_1... Kraft(t)
Q_1=0.01;
%Q_2... Kraft(u,u_dot)
Q_2=0.01;

%n...Anzahl der Moden
n = 3;
[w, phi] = EVP(n);

moden = phi(0.3);

%% Einfache Lösung der linearen DFG um ANfangswert zu bekommen
homvars.hom=0;
%sol0 ist die Menge der x und der y (x, y, phi, M, V) und der Parameter (P - Horizontale Kraft)
sol0.x=linspace(0,1,11);
sol0.y=zeros(5, size(sol0.x,2));
sol0.parameters=[P0];
opts=bvpset('RelTol', 1.0e-9);
homvars.yalt=zeros(6,1);
homvars.yalt(6)=P0;
sol1=bvp5c(@modellbalken5dgl, @modellbalken5bc, sol0, opts);
figure(1)
plot(sol1.y(1,:), sol1.y(2,:), 'b-')
xlabel("x (dimensionslos)")
ylabel("Durchbiegung (dimensionslos)")
hold off
%% Fortsetzungsmethode
% % homvars.hom=1;
%yalt ist der Lösungsvektor des ersten Schrittes auf der Linken Seite
% % homvars.yalt=[sol1.y(:,1);sol1.parameters];
%Note: In hom4dgl wird die Schrittweite von P geändert als auch die Anzahl
%der AUfteilungen des Balkens (Delta x...kann kleiner oder größer werden)
% % homsol=hom4dgl(@modellbalken5dgl,@modellbalken5bc,sol1,opts,6,0.1,1200);
%homsol(6,:) ist die horizontale Kraft
%homsol(4,:) ist das Biegemoment
%plot(homsol(6,:), homsol(4,:), 'b-');
%xlabel("Horizontale Kraft P (dimensionslos)")
%ylabel("Biegemoment")
savehomsol=homsols;


%% Suche von Umkehrpunkten und die Lösung dort erstellen!
% Bei den Umkehrpunkten ist die Jacobimatrix des Problemes singulär
% das bedeutet man findet keine Lösung für das Randwertproblem bei diesem
% Umkehrpunkt! Um jedoch eine Lösung zu finden muss das Gleichungssystem
% erweitert werden um die Gleichungen (2.7a + b+ c) aus dem Skript Ist wichtig für Stabilitätsuntersuchung,
% siehe Seite 26 im Skript!
% nk=size(homsol,2);
% nc=0;
% solc=[];
% ic=[];
% for i=2:nk-1
%     %Wann wird die Kraft wieder kleiner - Extremum von P gesucht?
%     if ((homsol(6,i+1)-homsol(6,i))*(homsol(6,i)-homsol(6,i-1)) <=0)
%         display(i)
%         solc0=savehomsol{i};
%         solc0.y(6:10,:)=0.1;
%         solc0.y(9,1)=1;
%     nc=nc+1;
%     ic(nc)=i;
%     Umkehrpunkt_solution = solc0.y(:,1)
%     homvars.yalt = [Umkehrpunkt_solution(1);Umkehrpunkt_solution(2);Umkehrpunkt_solution(3);Umkehrpunkt_solution(4);Umkehrpunkt_solution(5);... 
%         cos(Umkehrpunkt_solution(3)); sin(Umkehrpunkt_solution(3)); -Umkehrpunkt_solution(4); Umkehrpunkt_solution(5)*cos(Umkehrpunkt_solution(3))+solc0.parameters*sin(Umkehrpunkt_solution(3)); -q ;solc0.parameters]
%     solcrit=bvp4c(@modellbalken5singdgl, @modellbalken5singbc, solc0, opts);
%     solvec=[solcrit.y(:,1); solcrit.parameters];
%     solc=[solc,solvec];
%     solcritmat{nc}=solcrit;
%     end
% end
% figure(2)
% plot(homsol(6,:), homsol(4,:), 'b-', solc(11,:), solc(4,:), 'ro')
% xlabel("Horizontale Kraft P (dimensionslos)")
% ylabel("Biegemoment (dimensionslos)")
% hold off
% 
% figure(3)
% plot(solcrit.y(1,:), solcrit.y(2,:), 'b-')
% xlabel("x in meter")
% ylabel("Durchbiegung (dimensionslos)")
% title(sprintf("Durchbiegung bei P=" + solcrit.parameters) + " (dimensionslos)")
% %sqr = @(x) x.^2;
