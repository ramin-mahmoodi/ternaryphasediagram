%% INITIALIZE MATLAB
clear;
clc;
close all;
format long;

%% DASHBOARD
Z = [0.5 1e-20 0.5]; % Component Mole Fraction
P = input('Enter Pressure value (MPa): '); % MPa
T = input('Enter Temperature value (K): '); % Kelvin
R = 8.31446261815324e-3; % L.MPa/K/mol 
Curve_i = [];
Curve_j = [];
Curve_k = [];

% Critical Parameters and Acentric Factor
Pc = [4.599 3.796 2.110]; % MPa
Tc = [190.56 425.12 617.7]; % Kelvin
Pr = P./Pc;
Tr = T./Tc;
omega = [0.0115 0.2002 0.4923];

%% 1. Calculating Ki using Wilson eq
K = (1./Pr).*exp(5.37.*(1+omega).*(1-(1./Tr)));
x=0;
y=1;
while abs(x(1)-y(1))>0.01
%% 2. Calculate x & y & nV
err1 = inf;
iter = 0;
while err1 > 1e-12
    % calculating nV, x, y
    nV = 0.5; % guess
    err = inf;
    while err>1e-6
        F = sum(Z.*(K-1)./(1 + nV.*(K-1)));
        dF = -sum((Z.*(K-1).^2)./((1+nV.*(K-1)).^2));
        nV_old = nV;
        nV = nV - F/dF;
        err = abs(nV-nV_old);
    end
    x = Z./(1+(K-1).*nV);
    y = K.*Z./(1+(K-1).*nV);
    
    %% 3. PR EOS For Vapor phase
    %%%%%%%%%%%% Vapor %%%%%%%%%%%%%
    m = 0.37464 + 1.54226.*omega - 0.26992.*omega.^2;
    a_vector = .45724.*(R.^2).*(Tc.^2)./Pc.*(1+ m.*(1-sqrt(Tr))).^2;
    b_vector = .07780*R*Tc./Pc;
    
    % % Calculating a & b using Mixing Rule % %
    % Binary Interaction Parameters;
    BIP     =                [0         0.02     0.04;
                              0.02      0           0;
                              0.04      0           0];
            
    b = sum(y.*b_vector);     
    a = 0;
    Sum = 0;
    for i = 1:length(Z)
        for j = 1:length(Z)
            Sum = y(i)*y(j)*((a_vector(i)*a_vector(j))^0.5)*(1-BIP(i,j));
            a = Sum + a;
        end
    end
    
    %% 4. Calculating component Fugacity in gas phase 
    A = a.*P./((R.*T).^2);
    B = b.*P./(R.*T);
    
    r = roots([1 B-1 A-2*B-3*B^2 -(A*B-B^2-B^3)]);
    r = r(imag(r)==0);
    Zg=min(r);
    
    phi_G = 0;
    for i = 1:3
        sigma = 0; Bi = b_vector(i)*P/R/T;
        for j = 1:3
            sigma =  sigma + y(j)*(1-BIP(i,j))*(a_vector(i)*a_vector(j))^0.5;
        end
        phi_G(i) = exp(Bi./B*(Zg-1)-log(Zg-B)+A./(B.*(-2*sqrt(2))).*(Bi./B-2/a*sigma)*log((Zg+(1-sqrt(2)).*B)./(Zg+(1+sqrt(2)).*B)));
    end
    
    fg = phi_G.*y.*P;
    
    %% 5. PR EOS For Liquid phase
    %%%%%%%%%%%%%% Liquid %%%%%%%%%%%%%
    
    bl = sum(x.*b_vector);
    
    % a for liquid phase
    al = 0;        
    Sum = 0;
    for i = 1:length(Z)
        for j = 1:length(Z)
            Sum = x(i)*x(j)*((a_vector(i)*a_vector(j))^0.5)*(1-BIP(i,j));
            al = Sum + al;
        end
    end 
    

    %% 6. Calculating component Fugacity in liquid phase 
    A = al.*P./((R.*T).^2);
    B = bl.*P./(R.*T);
    r = roots([1 B-1 A-2*B-3*B^2 -(A*B-B^2-B^3)]);
    r = r(imag(r)==0);
    Zl=max(r);
    
    phi_L = 0;
    for i = 1:3
        sigma = 0; Bi = b_vector(i)*P/(R*T);
        for j = 1:3
            sigma =  sigma + x(j)*(1-BIP(i,j))*(a_vector(i)*a_vector(j))^0.5;
        end
        phi_L(i) = exp(Bi./B*(Zl-1)-log(Zl-B)+A./(B.*(-2*sqrt(2))).* ...
            (Bi./B-2/al*sigma)*log((Zl+(1-sqrt(2)).*B)./(Zl+(1+sqrt(2)).*B)));
    end
    
    fl = phi_L.*x.*P;
    
    %% 7. Determining Error Criteria and Calculating the New K
    err1 = sum((1 - fl./fg).^2);
    K = K.*(fl./fg);
    iter = iter + 1;
end

%% Table of Results
sprintf(' The results after %d iterations are: ', iter)
disp(['gas phase molar amount (nV): ', num2str(nV)]);
disp(['liquid phase molar amount (nL): ', num2str(1-nV)]);

VarNames = {'GasComposition', 'LiquidComposition'};
names = {'C1', 'nC4', 'nC10'};
results = table(y', x');
results.Properties.RowNames = names;
results.Properties.VariableNames = VarNames;
disp(results);

Z(1)=(x(1)+y(1))/2;
Z(2)=Z(2)+0.005;
Z(3)=1-(Z(1)+Z(2));
Curve_i(:, end+1) = x(1);
Curve_j(:, end+1) = x(2);
Curve_k(:, end+1) = x(3);
Curve_i(:, end+1) = y(1);
Curve_j(:, end+1) = y(2);
Curve_k(:, end+1) = y(3);
end
ternplot(Curve_k, Curve_j, Curve_i)
ternlabel('C10', 'nC4', 'C1')