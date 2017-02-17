function f = objfunc_initial_3truss(x0,flag,count)

%count=0;

% Material Variables
material = x0(4:6);

% Possible material properties for each element
for k = 1:3
    if material(k) == 1      % Aluminum
        E(k) = 68.9e9;         % [Pa]
        rho(k) = 2700;         % [kg/m^3]      
        sigma_y(k) = 55.2e6;   % [Pa]
    elseif material(k) == 2  % Titanium
        E(k) = 116e9;          % [Pa]
        rho(k) = 4500;         % [kg/m^3]
        sigma_y(k) = 140e6;    % [Pa]
    elseif material(k) == 3  % Steel
        E(k) = 205e9;          % [Pa]
        rho(k) = 7872;         % [kg/m^3]
        sigma_y(k) = 285e6;    % [Pa]
    elseif material(k) == 4  % Nickel
        E(k) = 207e9;          % [Pa]
        rho(k) = 8800;         % [kg/m^3]
        sigma_y(k) = 59.0e6;   % [Pa]
    end
end

% Area Variables [m^2]
A = x0(1:3)/1e4;

% Lengths of Elements
L(1) = sqrt(1.2^2+1.2^2);	% length of element 1 [m]
L(2) = 1.2;               % length of element 2 [m]
L(3) = sqrt(1.2^2+1.2^2);   % length of element 3 [m]

% Objective Functions
W = sum(rho.*A.*L);

% Stress Values
[sigma,dL] = stress_3truss(A,E);
%[sigma,DL,counter] = stressHW3(A,E,counter);


% Constraints
%g = abs(sigma) ./ sigma_y - 1;

%no penalty in hybridga
% %step-linear penalty function
% P = 0;
% for i = 1:length(g)
%     if g(i) >= 0
%         P = P + 1e6*(1+g(i));
%     end
% end

count = count + 1; %count of no. of times objfunc is used
%count = count+counter+1;

% Fitness Functions
phi1 = W; %weight
phi2 = dL; %displacement
f = [phi1 phi2];