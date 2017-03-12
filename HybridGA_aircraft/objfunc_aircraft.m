function f = objfunc_aircraft(x_con_hat,x_dis,lb_con,ub_con,Filename)
if size(x_con_hat,1)>1
    x_con_hat = x_con_hat';
end
x_con = lb_con + x_con_hat.*(ub_con-lb_con);

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
%NOX = 650;
TOC = 1e5;
if nan_count==0 && ac_count == 4
    %multi-objective functions
    FUEL = Outputs.FUEL;
    TOC = Outputs.TOC;
    %NOX = Outputs.NOX;
else
%     fprintf('\n%s\n','Mission failed!')
end
% Fitness Functions
phi1 = FUEL;
%phi2 = NOX;
phi3 = TOC;

f = [phi1 phi3]; %1. Fuel burn vs TOC
%f = [phi1 phi2]; %1. fuel burn vs NOx
%f = [phi2 phi3]; %3. NOx vs TOC