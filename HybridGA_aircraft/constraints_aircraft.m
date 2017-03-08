function [g, h] = constraints_aircraft(x_con,x_dis, Filename)
%counter=0;

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

[Outputs, nan_count, ac_count] = ReadFLOPSOutput(Filename);

TD = 10000;
LD = 10000;
if nan_count == 0 && ac_count == 4
	TD = Outputs.TD;
    LD = Outputs.LD;
else
%     fprintf('\n%s\n','Mission failed!')
end
g(1) = TD/8500 - 1; %take-off distance
g(2) = LD/7000 - 1; %landing distance

h = [];
end