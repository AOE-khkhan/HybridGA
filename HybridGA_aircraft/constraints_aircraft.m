function [g, h] = constraints_aircraft(x_con,x_dis)
%counter=0;

%calling FLOPS and getting Outputs
Inputs.seats = 162;
Inputs.GW = 174200;
Inputs.DESRNG = 2940;

Filename = 'aircraft_cons';

output = analyze_discrete(x_con,x_dis);

FLOPSInputGen(x_con,output,Inputs,Filename)
cmmndline = ['flops < ' Filename '.in> ' Filename '.out'];
[s,w] = dos(cmmndline); 
if s==1; 
   disp(w); 
end

Outputs = ReadFLOPSOutput(Filename);

%constraints
if isnan(Outputs.TD);
    TD = 10000; %high takeoff distance as penalty
else TD = Outputs.TD;
end

if isnan(Outputs.LD);
    LD = 10000; %high landing distance as penalty
else LD = Outputs.LD;
end
    
g(1) = TD/8500 - 1; %take-off distance
g(2) = LD/6500 - 1; %landing distance

h = [];
end