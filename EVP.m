function [w, phi] = EVP(n)
%Funktion um Eigenfrequenzen (Kreisfrequenzen) und Eigenmoden zu berechnen
%Hier handelt es sich um das linearisierte und das ungedämpfte System
global E Q_0 rho A L I n

i = linspace(1,n, n);
w = sqrt(  (E*I/rho*A) .* (i .* pi/L).^4 + (Q_0/rho*A) .* (i .* pi/L).^2  );
phi = @(x) sin(i .* pi*x);

end