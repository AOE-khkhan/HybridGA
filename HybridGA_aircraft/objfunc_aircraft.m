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

[Outputs, nan_count] = ReadFLOPSOutput(Filename);

if nan_count==0
    %multi-objective functions
    if isnan(Outputs.FUEL);
        FUEL = 55000; %high fuel burn as penalty
    else
        FUEL = Outputs.FUEL;
    end

    if isnan(Outputs.NOX);
        NOX = 450; %high NOX value as penalty
    else
        NOX = Outputs.NOX;
    end
else
    FUEL = 55000;
    NOX = 450;
end
 
% if isnan(Outputs.FARE)
%     FARE = 450; %high fare as penalty
% else FARE = Outputs.FARE;
% end

% %check
% if isnan(Outputs.TOC);
%     TOC = 50000; %high fare as penalty
% else TOC = Outputs.TOC;
% end

% Fitness Functions
phi1 = FUEL; 
phi2 = NOX;
%phi3 = FARE;

%phi4 = TOC;

f = [phi1 phi2]; %1. fuel burn vs NOx
%f = [phi1 phi3]; %2. fuel burn vs ticket price
%f = [phi2 phi3]; %3. NOx vs ticket price


%f = [phi1 phi4]; %check
