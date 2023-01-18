function [w] = modellbalken5bc(ya,yb,parameters)
%modellbalken5bc Randbedingungen für stetiges Homotopiebeispiel
%   mit Fortsetzungsmethode
global q homvars;
if (homvars.hom==0)
    bchom=parameters(1)-homvars.yalt(6);
else
    yneu=[ya; parameters];
    bchom=norm(yneu-homvars.yalt)-homvars.step;
end
w=[ya(1); ... % x(0)=0
    ya(2); ... % y(0)=0
    ya(3); ... % phi(0)=0
    yb(2); ... % y(1)=0
    yb(4); ... % M(1)=0
    bchom];
end

