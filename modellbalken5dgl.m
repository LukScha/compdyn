function [dy] = modellbalken5dgl(s, y, parameters)
%modellbalken5dgl: Differentialgleichungen für stetigen Modellbalken
%   Detailed explanation goes here
global q;
p=parameters(1);
phi=y(3);
dy=[cos(phi); sin(phi); -y(4); y(5)*cos(phi)+p*sin(phi); -q];
end

