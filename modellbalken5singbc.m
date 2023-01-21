function [w] = modellbalken5singbc(ya,yb,parameters)
%modellbalken5bc Randbedingungen für stetiges Homotopiebeispiel
%   mit Fortsetzungsmethode HALLO
global q homvars;
if (homvars.hom==0)
    bchom=parameters(1)-homvars.yalt(6);
else
    yneu=[ya; parameters];
    bchom=norm(yneu-homvars.yalt)-homvars.step;
end
w=[ya(1); ... % x(0)=0
    ya(2); ... % y(0)=0
    ya(3);... % phi(0)=0
    ya(6)-1;... % x'(0)=cos(ya(3)) = 1
    ya(7)-0;... % y'(0)=sin(ya(3)) = 0
    ya(8)+ya(4);... % phi'(0)=-ya(4)
    yb(2); ... % y(1)=0
    yb(4); ... % M(1)=0
    yb(7)-sin(yb(3)); ... % y'(1)=sin(yb(3))
    yb(8)+yb(4); ... % phi'(1)=-yb(4)
    bchom;];
end
