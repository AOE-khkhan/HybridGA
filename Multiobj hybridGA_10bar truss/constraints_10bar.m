function [g, h] = constraints_10bar(x_con,x_dis)
%counter=0;

%converting x_con, x_dis to x0
i = 1; 
for counter = 1:2:20
   x0(counter) = x_dis(i);  
   x0(counter+1) = x_con(i);
   i=i+1;
end

if(length(x0)~=20)
    error('Vector needs 20 inputs arguments');
end

%Sets the maximum displacement
displacementMax = 2;
%Sets the maximum stresses
stress_max=[25 20 13.3 73.3]*10^3;

Einit=[1e7 2.8e7 3e7 1.65e7];

i = 1; 
for counter = 1:2:20
   E(i) =  Einit(x0(counter));  
   A(i) = x0(counter+1);
   i=i+1;
end

[displacement,stress] = tenBarTruss(A,E);

for i=1:10
    c(i) =  abs(stress(i))./stress_max(x0(2*i-1))-1;
end

c(11) = max(abs(displacement))/displacementMax-1;
g = c(11);

%ceq = [];
h = [];
end