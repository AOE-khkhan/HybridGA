function f = objfunc_initial_10bar(x0,flag,count)

if(length(x0)~=20)
   error('Vector needs 20 inputs arguments');
end

rhoInit=[0.10 0.284 0.318 0.162]; % Note: some people use 0.101 for the density of aluminium
Einit=[1e7 2.8e7 3e7 1.65e7];

i = 1; 
for counter = 1:2:20 
   rho(i) =  rhoInit(x0(counter));
   A(i) = x0(counter+1);
   E(i) =  Einit(x0(counter));
   i=i+1;
end

W = calculateWeight(A,rho);
displacement = tenBarTruss(A,E);

%[c] = constraintsTenBarTruss(x);

dL = sqrt(displacement(1).^2 + displacement(2).^2); %resultant displacement

count = count + 1; %count of no. of times objfunc is used
%count = count+counter+1;

% Fitness Functions
phi1 = W; %weight
phi2 = dL; %displacement

f = [phi1 phi2];