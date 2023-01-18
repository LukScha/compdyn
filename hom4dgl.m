function [sol, lastsol] = hom4dgl(odefun, bcfun, solinit, bvpoptions, firstk, firststep, ...
    nsteps, homoptions)
%Variante der Fortsetzungsmethode für Differentialgleichungen
%mit Phasenbedingungen
% (y_neu - y_alt)^2 = d^2.
%
% Um dieses Verfahren zu verwenden, sollte in der Function, die die
% Randbedingungen beschreibt, die Struktur homvars mit folgenden
% Elementen global vereinbart werden:
%
% homvars.hom: Wird innerhalb der Methode auf 1 gesetzt. In diesem
% Fall muss die Homotopie-bedingung ausgewertet werden. Andere
% Werte können zB. verwendet werden, um einen Startpunkt zu
% berechnen.

% homvars.yalt:  Vektor der vorigen Lösung am linken Rand;
% besteht aus den unbekannten Funktionswerten y(a) und den freien
% Parametern.

% homvars.step: Schrittweite

% Siehe zB. duffing_bc.m, wie die Phasenbedingung spezifiziert
% wird.

global homvars;
global ghomsol; % Bisher berechnete Lösung. Nützlich, falls hom4dgl
% abbricht
global homsols; % Die diversen Lösungsdaten, nicht nur Anfangswerte.
% Kann viel Speicher kosten
  global homstartsols; % Startvektoren
  homstartsols={};
lastsolx=solinit;
sol = [];
ghomsol = [];
homsols = [];
nx=length(solinit.y(:,1));
npar=length(solinit.parameters);
if (firstk < 0 || firstk > nx+npar) % Index k out of range
    if (npar > 0)
        firstk = nx+1;
    else
        firstk = 1;
    end
end
k=firstk;
step=abs(firststep);
%if (nargin < 8) % Default-Optionen
hom.epsmin = 1.0e-4; % Minimale Schrittweite
hom.itaim = 8;       % Optimale # Newtonschritte
hom.gew =       1;       % Relatives Gewicht des Parameters
hom.iex =       1;       % Extrapolation
hom.alphamax = 1.1;
hom.fid =       1;
%else
if (nargin > 7)
    %  hom = homoptions;        % Vorhandene Optionen aus Parameterliste übernehmen
    optnames=fieldnames(homoptions);
for iopt=1:size(optnames,1)
	   optname=optnames{iopt,1};
if (isfield(hom, optname))
  hom.(optname)=homoptions.(optname);
        end
        %hom=homoptions;
    end
    if (hom.epsmin <= eps)
        hom.epsmin = 1.0e-4;
    end
    if (hom.itaim <= 0)
        hom.itaim = 8;
    end
    if (hom.iex) < 0
        hom.iex = 0;
    end
    if (hom.alphamax < 1)
        hom.alphamax = 1.1;
    end
    
end
    alpha=0;
    for ihom = 1:nsteps
    if (npar > 0)
        homvars.yalt=[solinit.y(:,1); solinit.parameters];
    else
        homvars.yalt=solinit.y(:,1);
    end
    if (ihom > 1)
        alpha = 1; % Keine Schrittweitenkontrolle möglich
        if (alpha > hom.alphamax)
            alpha = hom.alphamax;
        end
        step = step*alpha;
    end
    homvars.step=step;
    ifail = 0;
    for iinner=1:3
        solx.x = solinit.x;
        if (ihom > 1 & hom.iex == 1)
            ik=1;
            solx.y=zeros(size(solinit.y,1), size(solx.x,2));
            nkalt=length(solalt.x);
            for i=1:length(solx.x)
                for j=ik:nkalt
                    if (solalt.x(j)==solx.x(i))
                        yi=solalt.y(:,j);
                        ik=j+1;
                        break;
                    elseif (solalt.x(j) > solx.x(i))
                        if (isfield(solalt,'idata'))
                            yi = deval(solinit.x(i), solalt); % Alte Loesung auf neuem
                        % Gitterpunkt
                        else
                            yi=solalt.y(:,j);
                        end
                        break;
                    else
                        ik=j;
                    end
                end
                yni = solinit.y(:,i);
                solx.y(:, i) = yni*(1+alpha) - alpha*yi;  % Neuer Startwert
            end
            solx.parameters = solinit.parameters*(1+alpha) -  ...
                alpha*solalt.parameters;
        else
            solx = solinit;
            if (firstk <= nx)
                solx.y(firstk, 1) = solx.y(firstk, 1) + firststep;
            else
                solx.parameters(firstk-nx) = solx.parameters(firstk-nx) + ...
                    firststep;
            end
        end
        homvars.step=step;
        fprintf(hom.fid, 'ihom: %i, step: %g, alpha: %g, iinner: %d, start:\n', ihom, homvars.step, alpha, iinner);
        fprintf(hom.fid, ' %15.7g %15.7g %15.7g %15.7g %15.7g\n', solx.y(:,1));
        fprintf(hom.fid, ' %15.7g %15.7g %15.7g %15.7g %15.7g\n', solx.parameters);
%        fprintf(hom.fid, '\nhomvars.yalt:\n');
%        fprintf(hom.fid, '%15.7g ', homvars.yalt); fprintf(hom.fid, '\n');
        homstartsols=[homstartsols, solx];
        [solx, sol_flag, iter] = solver(odefun, bcfun, solx, bvpoptions);
        if (sol_flag == 0)  % Loesung gefunden
            lastsolx = solx;
            solalt = solinit;
            if (ihom==1)
                solalt.solver=solx.solver;
            end
            solinit = solx;
            solvec=[solx.y(:,1); solx.parameters];
            fprintf(hom.fid, ' %15.7g %15.7g %15.7g %15.7g %15.7g\n', solvec);
            fprintf(hom.fid, '\n');
            sol = [sol, solvec];
            ghomsol = [ghomsol, solvec];
            homsols{ihom} = solx;
            if (ihom < nsteps)
                if (hom.iex == 0)
                    break;        % Beende inneren Loop
                end
                
                
            end            % ihom < nsteps
            break;         % continue outer loop
        else             % sol_flag ~= 0: Fehlerbehandlung
            fprintf(hom.fid, 'flag = %d iinner = %d \n', sol_flag, iinner);
            if (iinner < 4 & step > 5*hom.epsmin)
                step = 0.2*step; % Verkürze die Schrittweite
                alpha = 0.2*alpha;
                continue;
            else
                if (nargout > 1)
                    flag = -1;
                end
                return;
            end
        end % Sol_flag
    end
end % Outer loop
if (nargout > 1)
    lastsol = lastsolx;
end


% Innere Funktion: Interface zu Randwertlöser
function [x, flag, iter] = solver(odefun, bcfun, x0, bvpoptions)
Nmax=500;
RelTol=1.0e-5;
if (nargin > 3)
    if (isfield(bvpoptions, 'Nmax'))
        Nmax=bvpoptions.Nmax;
    end
    if (isfield(bvpoptions, 'RelTol'))
        RelTol=bvpoptions.RelTol;
    end
end

options=bvpset(bvpoptions, 'RelTol', RelTol, 'NMax', Nmax, 'Stats', 1);
bvperr=0;
try
    x = bvp5c(odefun, bcfun, x0, options);
catch
    bvperr=1;
    x=x0;
end
if (bvperr == 0)
    ya=x.y(:,1);
    nswpoint=0;
    for i=1:size(x.x,2)-1
        if (x.x(i)==x.x(i+1)) % Switching Point
            nswpoint=nswpoint+1;
            ya=[ya, x.y(:,i+1)];
            yb(:,nswpoint)=x.y(:,i);
        end
    end
    yb(:,nswpoint+1)=x.y(:,end);
    bcres=bcfun(ya, yb, x.parameters);
    %  fprintf(1, 'Mesh size: %d BC Residual: ', size(x.y, 2));
    %  fprintf(1, ' %14.7g', bcres); fprintf(1, '\n');
    if norm(bcres) > 2.0e-3
        flag = -2;
    else
        flag = 0; % Diese beiden Werte werden in bvp4c leider nicht
        % gesetzt.
    end
else
    flag=-1;
end
%  plot (x.x, x.y(1,:))
iter = 1; % Dadurch ist die Fehler- und Schrittweitenkontrolle
% unwirksam.
