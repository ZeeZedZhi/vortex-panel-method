function [Cp, gammas] = vortex_panel_method(X, Y, alpha)
% VORTEX_PANEL_METHOD The coefficients of pressure and gammas at the
%   midpoint of panels
%   VORTEX_PANEL_METHOD(X, Y, alpha) returns the coefficients of 
%   pressure and gammas at the midpoints of panels whose nodes'
%   coordinates indicated with vectors X and Y. alpha is the angle
%   of attack given in radians.

    n = length(X) - 1;
    XC = (X(1:end-1) + X(2:end)) / 2;
    YC = (Y(1:end-1) + Y(2:end)) / 2;
    thetas = atan2(Y(2:end) - Y(1:end-1), X(2:end) - X(1:end-1));
    S = sqrt((X(2:end) - X(1:end-1)).^2 + (Y(2:end) - Y(1:end-1)).^2);
    
    A = zeros(n);
    B = zeros(n);
    C = zeros(n);
    D = zeros(n);
    E = zeros(n);
    F = zeros(n);
    G = zeros(n);
    P = zeros(n);
    Q = zeros(n);
    
    for i = 1:n
        for j = 1:n
            A(i, j) = (X(j)-XC(i))*cos(thetas(j)) + (Y(j)-YC(i))*sin(thetas(j));
            B(i, j) = (X(j)-XC(i))^2 + (Y(j)-YC(i))^2;
            C(i, j) = sin(thetas(i) - thetas(j));
            D(i, j) = cos(thetas(i) - thetas(j));
            E(i, j) = (Y(j)-YC(i))*cos(thetas(j)) - (X(j)-XC(i))*sin(thetas(j));
            F(i, j) = log(1. + (S(j) + 2*A(i, j))*S(j) / B(i, j));
            G(i, j) = atan2(E(i, j)*S(j), B(i, j) + A(i, j)*S(j));
            P(i, j) = (XC(i)-X(j))*sin(thetas(i) - 2*thetas(j)) + ... 
                (YC(i)-Y(j))*cos(thetas(i) - 2*thetas(j));
            Q(i, j) = (XC(i)-X(j))*cos(thetas(i) - 2*thetas(j)) - ... 
                (YC(i)-Y(j))*sin(thetas(i) - 2*thetas(j));
        end
    end
    
    Cn2 = D + ((.5)*Q.*F - (A.*C + D.*E).*G)./S;
    Cn1 = (.5).*D .*F + C.*G - Cn2;
    
    An = zeros(n+1);
    for i = 1:n
        An(i, 1) = Cn1(i, 1);
        An(i, n+1) = Cn2(i, n);
        for j = 2:n
            An(i, j) = Cn1(i, j) + Cn2(i, j-1);
        end
    end
    An(n+1, 1) = 1;
    An(n+1, n+1) = 1;
    
    RHS = [sin(thetas - alpha), 0]';
    
    gammas = An \ RHS;
    
    Ct2 = C + (.5*P.*F + (A.*D - C.*E).*G)./S;
    Ct1 = (.5)* C.*F - D.*G - Ct2;
    for i = 1:n
        Ct1(i, i) = pi / 2;
        Ct2(i, i) = pi / 2;
    end
    
    At = zeros(n, n+1);
    for i = 1:n
        At(i, 1) = Ct1(i, 1);
        At(i, n+1) = Ct2(i, n);
        for j = 2:n
            At(i, j) = Ct1(i, j) + Ct2(i, j-1);
        end
    end
    
    V = cos(thetas - alpha)' + At*gammas;
    Cp = 1 - V.^2;

end