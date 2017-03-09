function f = objfunc_aircraft(x_con,x_dis,Filename)

%calling FLOPS and getting Outputs
Inputs.seats = 162;
Inputs.GW = 174200;
Inputs.DESRNG = 2940;

output = analyze_discrete(x_con,x_dis);

FLOPSInputGen(x_con,output,Inputs,Filename)
cmmndline = ['flops < ' Filename '.in> ' Filename '.out'];
[s,w] = dos(cmmndline); 
if s==1; 
   disp(w); 
end

[Outputs,nan_count,ac_count] = ReadFLOPSOutput(Filename);

FUEL = 55000;
NOX = 650;
FARE = 1e3;
if nan_count==0 && ac_count == 4
    %multi-objective functions
    FUEL = Outputs.FUEL;
    NOX = Outputs.NOX;
    FARE = Outputs.FARE;
else
%     fprintf('\n%s\n','Mission failed!')
end
% Fitness Functions
phi1 = FUEL;
phi2 = FARE;
%phi3 = FARE;

%phi4 = TOC;

f = [phi1 phi2]; %1. fuel burn vs NOx
%f = [phi1 phi3]; %2. fuel burn vs ticket price
%f = [phi2 phi3]; %3. NOx vs ticket price


%f = [phi1 phi4]; %check