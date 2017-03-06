clc;
clear all;
close all;
Timestart = tic;
Start_time=clock;
format longg;
matlabpool 8
options = goptions([]);
%x0=ones(48,12); %population X variables - CHECK FOR EVERY CHANGE IN DESIGN VARIABLES/POP. SIZE
ND0=[];xopt=[];xopt0=[];
count=0;

%Initial population for GA
for ii=1:48
    x0(ii,:) = [9.4;0.159;0.1338;1345.5;25;24200;0;0;0;0;1;1];
end

%Sequence of the design variables
%Continuous: x_con
    %1. Aspect ratio AR (5 bits res 0.129)
    %2. Taper ratio TR (5 bits res 0.0129)
    %3. Thickness to chord ratio TCA (5 bits res 0.0026)
    %4. Wing area[ft^2] S (5 bits res=16.13)
    %5. Wing sweep at 25% chord[deg] Sweep (5 bits res 1.29)
    %6. Thrust per engine[lbs] Thrust (5 bits res 322.581)
    bits_con = [5 5 5 5 5 5]; %Bits for continuous variable
    lb_con = [8 0.1 0.09 1000 0 20000]; %Lower bounds for the continuous design variables
    ub_con = [12 0.5 0.17 1500 40 30000]; %Upper bound for the continuous design variables
    
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
    %%%%will be implemented later:
    %    9 - HLFC wing+tail+nacelle
    %   10 - NLF wing + HLFC tail
    %   11 - NLF wing + HLFC nacelle
    %   12 - NLF wing + HLFC tail + HLFC nacelle
    bits_dis = [1 1 1 1 3 3]; %change 3bits to 4bits
    lb_dis = [0 0 0 0 1 1];
    ub_dis = [1 1 1 1 8 8]; %change 8 to 16 

num_con = length(lb_con);
num_dis = length(lb_dis);
   
ToTFunc = 0;
 
%pair=3; %Select the objective pair

fmin = [8000, 2]; fmax=[100000, 500];       %1: Fuel burn VS NOX
% fmin = [8000, 100]; fmax=[100000, 1000];  %2: Fuel burn VS TicketPrice
% fmin=[2, 100]; fmax=[500, 1000];          %3: NOX VS TicketPrice

flag = 1;

[xopt,ND,funeval,lgen,fmin_hist] = hybridGA_aircraft(x0,ND0,xopt,options,lb_con,ub_con,bits_con,lb_dis,ub_dis,bits_dis,fmin,flag,count)

ToTFunc=ToTFunc+funeval;

% %Write results in text file
% Filename='output_aircraft';
% WriteResults(ToTFunc,xopt,ND,Filename)

%Write results in text file
Filename='output_aircraft';
WriteResults(ToTFunc,xopt,ND,Filename)
A = ToTFunc;
B = xopt;
C = ND;
fid=fopen('output_10bar.txt','w');
fprintf(fid, '%f\n',A);
fprintf(fid, '%f\n',B);
fprintf(fid, '%f\n',C);
fclose(fid);

%Sending an email upon completion
% string='Lam Max_gen=50, pop_size=50, obj_pair=3';
% sendemail(string)
Start_time
End_time=clock
Run_time = toc(Timestart)

%pareto front plot
figure
plot(ND(:,1),ND(:,2),'k.');
xlabel('f_1')
ylabel('f_2')

matlabpool close