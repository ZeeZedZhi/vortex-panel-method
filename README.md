# Vortex-Panel-Method
The Vortex Panel Method as seen in Kuethe and Chow's Foundations of Aerodynamics: Basics of Aerodynamic Design 5th Edition implemented in MATLAB.

## NACA Four-Digit Airfoil Coordinates Example
Calculation of the x and y coordinates ```X``` and ```Y``` of a NACA four-digit airfoil's profile and the y coordinates of its camber line ```yc```.
```
digits = 2412;
n = 48;
[X, Y, yc] = naca4digit(digits, n, true);
```
```digits``` specifies the four-digit airfoil. ```n``` specifies the number of panels. ```true``` indicates that the trailing edge should be smoothly closed.

<!---
close all
plot(X, Y)
hold on
domain = linspace(0, 1);
plot(domain, yc(domain));
legend("Profile", "Camber Line")
title(sprintf("NACA %d Contour", digits))
axis equal
xlim([0, 1])
ylim([-.16, .16])
-->

## Coefficients and Center of Pressure Example
Calculation of the coefficients of pressure ```Cp``` along with the lift ```Cls```, drag ```Cds```, and moment about the leading edge ```Cmles``` coefficients and the centers of pressure ```Xcps``` of a general simple closed loop.
```
alpha = 8*pi/180;
Cp = vortex_panel_method(X, Y, alpha);
[Cls, Cds, Cmles, Xcps] = get_coefficients(X, Y, alpha);
```

```alpha``` is the angle of attack in radians. ```X``` and ```Y``` are respectively the x and y coordinates of the loop's nodes.

<!---
resolutions = [12, 48, 120];
alpha = 8*pi/180;
digits = 2412;

close all
hold on
for n = resolutions
    [X, Y] = naca4digit(digits, n, true);
    Cp = vortex_panel_method(X, Y, alpha);
    
    XC = (X(1:end-1) + X(2:end)) / 2;
    if n ~= resolutions(end)
        scatter(XC, Cp, "filled")
    end
end

plot(XC(1:n/2), Cp(1:n/2), "Color", "#EDB120")
plot(XC(n/2+1:end), Cp(n/2+1:end), "Color", "#EDB120")
hold off

xlabel("x")
ylabel("Cp")
title(sprintf("NACA %d Coefficient of Pressures", digits))
legend(arrayfun(@num2str, resolutions, 'UniformOutput', 0))
set(gca, 'ydir', 'reverse')
-->

<!---
alphas = (-8:8)*pi/180;
alphaCount = length(alphas);
Cls = zeros(alphaCount, 1);
CM4s = zeros(alphaCount, 1);
Xcps = zeros(alphaCount, 1);

for j = 1:alphaCount
    [Cls(j), ~, ~, Xcps(j)] = get_coefficients(X, Y, alphas(j));
end
CM4s = -Cls .* (Xcps - .25);

alphas = alphas'*180/pi;
clfit = fit(alphas, Cls, "poly1")
cm4fit = fit(alphas, CM4s, "poly1")

alphaZLC = -clfit.p2 / clfit.p1;

xcpfit = fit(alphas, Xcps, sprintf("p1/(x-%f) + p2", alphaZLC))

dcmdcl = cm4fit.p1 / clfit.p1;
xac = .25 - dcmdcl;

macfit = fit(alphas, -Cls.*(Xcps - xac), "poly1")

alphaZLC, dcmdcl, xac

close all

scatter(alphas, Cls, "filled")
hold on
scatter(alphas, CM4s, "filled")
scatter(alphas, Xcps, "filled")
xlabel("Angle of Attack (°)")

domain = linspace(alphas(1), alphas(end));
plot(domain, polyval(coeffvalues(clfit), domain), "Color", "#0072BD")
plot(domain, polyval(coeffvalues(cm4fit), domain), "Color", "#D95319")

domain = linspace(alphas(1), alphaZLC-.5);
plot(domain, xcpfit.p1 ./ (domain - alphaZLC) + xcpfit.p2, "Color", "#EDB120")
domain = linspace(alphaZLC+.15, alphas(end));
plot(domain, xcpfit.p1 ./ (domain - alphaZLC) + xcpfit.p2, "Color", "#EDB120")
legend("cl", "cmc/4", "xcp")
-->
