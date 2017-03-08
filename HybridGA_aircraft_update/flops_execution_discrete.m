clc;clear all
Inputs.seats = 162;
Inputs.GW = 174200;
Inputs.DESRNG = 2940;
x_con = [9.4; 0.159; 0.1338; 1345.5; 25; 24200];
Filename = 'AC-X1';

x_dis = [0; 0; 0; 0; 1; 1];
output = analyze_discrete(x_con,x_dis);
%output.TRUW

FLOPSInputGen(x_con,output,Inputs,Filename)
cmmndline = ['flops < ' Filename '.in> ' Filename '.out'];
[s,w] = dos(cmmndline); 
if s==1; 
   disp(w); 
end

Outputs = ReadFLOPSOutput(Filename);

%fprintf('\n%s%f\n','TOC: ',Outputs.TOC)
%fprintf('\n%s%f\n','BH: ',Outputs.BH)
fprintf('\n%s%f\n','NOX: ',Outputs.NOX)
fprintf('\n%s%f\n','FUEL: ',Outputs.FUEL)
fprintf('\n%s%f\n','FARE: ',Outputs.FARE)