function [g,h] = constraints_3truss(x_con,x_dis)

%counter=0;

x0(1:3) = x_con;
x0(4:6) = x_dis;
% Material Variables
material = x0(4:6);

% Possible material properties for each element
for ii = 1:3
    if material(ii) == 1      % Aluminum
        E(ii) = 68.9e9;         % [Pa]
        rho(ii) = 2700;         % [kg/m^3]      
        sigma_y(ii) = 55.2e6;   % [Pa]
    elseif material(ii) == 2  % Titanium
        E(ii) = 116e9;          % [Pa]
        rho(ii) = 4500;         % [kg/m^3]
        sigma_y(ii) = 140e6;    % [Pa]
    elseif material(ii) == 3  % Steel
        E(ii) = 205e9;          % [Pa]
        rho(ii) = 7872;         % [kg/m^3]
        sigma_y(ii) = 285e6;    % [Pa]
    elseif material(ii) == 4  % Nickel
        E(ii) = 207e9;          % [Pa]
        rho(ii) = 8800;         % [kg/m^3]
        sigma_y(ii) = 59.0e6;   % [Pa]
    end
end

% Area Variables [m^2]
A = x0(1:3)/1e4;

% Lengths of Elements
L(1) = sqrt(1.2^2+1.2^2);	% length of element 1 [m]
L(2) = 1.2;               % length of element 2 [m]
L(3) = sqrt(1.2^2+1.2^2);   % length of element 3 [m]

% Stress Values
[sigma,DL] = stress_3truss(A,E);
%[sigma,DL,counter] = stressHW3(A,E,counter);

% Constraints
%const = abs(sigma) ./ sigma_y - 1;

g(1) = ((abs(sigma(1)))/sigma_y(1)) - 1;
g(2) = ((abs(sigma(2)))/sigma_y(2)) - 1;
g(3) = ((abs(sigma(3)))/sigma_y(3)) - 1;

g = [g(1) g(2) g(3)];

h = [];
end