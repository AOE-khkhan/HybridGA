clc;clear;close all;
Timestart = tic;
Start_time=clock;
format longg;

Filename_org = 'OutputFiles/AC_X';
options = goptions([]);
ND0=[];xopt=[];xopt0=[];
count=0;

%Sequence of the design variables
%Continuous: x_con
    %1. Aspect ratio AR (5 bits res 0.129)
    %2. Taper ratio TR (5 bits res 0.0064)
    %3. Thickness to chord ratio TCA (5 bits res 0.0026)
    %4. Wing area[ft^2] S (5 bits res=16.13)
    %5. Wing sweep at 25% chord[deg] Sweep (5 bits res 1.29)
    %6. Thrust per engine[lbs] Thrust (5 bits res 322.581)
    %7. By-Pass Ratio (BPR) (5 bits res 0.3225)
    %8. Turbine Inlet Temperature[R] (TIT) (5 bits res 24.839)
    %9. Operating Pressure Ratio (OPR) (5 bits res 0.645)
    %10.Fan Pressure Ratio (FPR) (5 bits res 0.0129)
    bits_con = [5 5 5 5 5 5 5 5 5 5]; %Bits for continuous variable
    % Make sure bounds are are always row vectors
    lb_con = [8 0.3 0.09 1000 0 20000 0 0 0 0]; %Lower bounds for the continuous design variables
    ub_con = [12 0.5 0.17 1500 40 30000 10 770 20 0.4]; %Upper bound for the continuous design variables
    
%Discrete: x_dis 
    % 1. Composites: x_dis(1:4)
    %    1 - 0/1 for wing
    %    2 - 0/1 for fuselage
    %    3 - 0/1 for nacelle
    %    4 - 0/1 for tail
    % 2. Position and number of engine/s: x_dis(5)
    %    1 - 2 wing 
    %    2 - 2 fuselage 
    %    3 - 2 wing + 1 fuselage
    %    4 - 3 fuselage
    %    5 - 4 wing 
    %    6 - 2 wing + 2 fuselage
    %    7 - 1 fuselage
    %    8 - 4 wing + 1 fuselage
    % 3. Laminar flow technologies: x_dis(6)
    %    1 - No laminar flow technologies
    %    2 - NLF wing
    %    3 - HLFC wing
    %    4 - HLFC nacelle
    %    5 - HLFC tail
    %    6 - HLFC wing + nacelle
    %    7 - HLFC wing + tail 
    %    8 - HLFC tail + nacelle
    %    9 - HLFC wing+tail+nacelle
    %   10 - NLF wing + HLFC tail
    %   11 - NLF wing + HLFC nacelle
    %   12 - NLF wing + HLFC tail + HLFC nacelle
    % 4. Engine Technologies
    %    1 - DDF
    %    2 - GTF
    %    3 - CRTF
    %    4 - OR
    bits_dis = [1 1 1 1 3 4 2];
    lb_dis = [0 0 0 0 1 1 1];
    ub_dis = [1 1 1 1 8 16 4];
    dis_cases = [2 2 2 2 8 12 4]; %ENTER NO. OF CASES BEING INPUT FOR EACH DISCRETE VARIABLE
    
    % Initial population for GA
    lb  = [lb_con, lb_dis];
    ub = [ub_con, ub_dis];
    x0 = zeros(48,17);
    x0(1,:) = [9.4;0.35;0.1338;1345.5;25;24200;5;250;10;0.05;0;0;0;0;1;1;1];
    for ii=2:48
        for jj = 1:size(x0,2)
            if jj<=10
                x0(ii,jj) = lb(jj) + rand*(ub(jj)-lb(jj));
            else
                x0(ii,jj) = randi([lb_dis(jj-10), ub_dis(jj-10)]);
            end
        end
    end
num_con = length(lb_con);
num_dis = length(lb_dis);
   
ToTFunc = 0;
 
%fmin = [8000, 2]; fmax=[100000, 1000];       %1: Fuel burn VS NOX
%fmin = [8000, 10000]; fmax=[100000, 1000000];  %2: Fuel burn VS TOC
fmin=[2, 10000]; fmax=[1000, 1000000];          %3: NOX VS TOC

flag = 1;

[xopt,ND,funeval,lgen,fmin_hist] = hybridGA_aircraft(x0,ND0,xopt,options,...
    lb_con,ub_con,bits_con,lb_dis,ub_dis,bits_dis,fmin,flag,count,Filename_org,dis_cases)

ToTFunc=ToTFunc+funeval;

cd OutputFiles
save('Results.mat','xopt','ND','funeval','lgen','fmin_hist')
cd ..

Start_time
End_time=clock
Run_time = toc(Timestart)

% Pareto front plot
figure
plot(ND(:,1),ND(:,2),'k.');
xlabel('f_1')
ylabel('f_2')
