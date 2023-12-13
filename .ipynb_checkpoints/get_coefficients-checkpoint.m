function [cl, cd, cmle, xcp] = get_coefficients(X, Y, alpha)
% get_coefficients The 2D coefficients of lift, drag, and moment
%   about the leading edge and center of pressure
%   get_coefficients(cl, cd, cmle, xcp) returns the coefficients of 
%   lift cl, drag cd, and moment about the leading edge cmle 
%   and center of pressure xcp given the coordinates of panels 
%   indicated with vectors X and Y. alpha is the angle of attack 
%   given in radians.   

    Cp = vortex_panel_method(X, Y, alpha);
    XC = (X(1:end-1) + X(2:end)) / 2;
    YC = (Y(1:end-1) + Y(2:end)) / 2;

    n = length(Cp);
    dyUP = diff(Y(n/2+1:end))';
    dxUP = diff(X(n/2+1:end))';
    
    dyDOWN = diff(Y(n/2+1:-1:1))';
    dxDOWN = diff(X(n/2+1:-1:1))';
    
    dydxUP = dyUP ./ dxUP;
    dydxDOWN = dyDOWN ./ dxDOWN;
    
    Clx = Cp(n/2:-1:1) - Cp(n/2+1:end);
    Cdx = Cp(n/2+1:end) .* dydxUP - Cp(n/2:-1:1) .* dydxDOWN; 
    
    CnCa = [trapz(XC(n/2+1:end), Clx); trapz(XC(n/2+1:end), Cdx)];
    ClCd = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*CnCa;
    cl = ClCd(1);
    cd = ClCd(2);
    
    Cmle1 = trapz(XC(n/2+1:end), -Clx .* XC(n/2+1:end)');
    Cmle2 = trapz(XC(n/2+1:end), Cp(n/2+1:end) .* YC(n/2+1:end)' .* dydxUP);
    Cmle3 = -trapz(XC(n/2+1:end), Cp(n/2:-1:1) .* YC(1:n/2)' .* dydxDOWN);
    CmleXcp = (Cmle1 + Cmle2 + Cmle3) .* [1; -1 / cl];
    cmle = CmleXcp(1);
    xcp = CmleXcp(2);
end