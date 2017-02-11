clc;
clear all;
close all;
Timestart = tic;
Start_time=clock;
format longg;
matlabpool 8
options = goptions([]);
x0=[];ND0=[];xopt0=[];
load Model
%Generating initial population for GA
for ii=1:(8*6)
    AR_w=8+(4*rand);TR_w=0.1+(.5-.1)*rand;t2c_w=0.09+(0.17-0.09)*rand;
    S_w=1000+(1500-1000)*rand;LamdaDeg_w=0+40*rand;
    BPR=10*rand; TET=500*rand;OPR=20*rand;
    FPR=0+0.1*rand;thrust=20000+(30000-20000)*rand; 
    x_con=[AR_w,TR_w,t2c_w,S_w,LamdaDeg_w,BPR,TET,OPR,FPR,thrust]; 
    x_dis=[randi([2,4]),randi([0,1]),randi([0,1]),randi([3,8]),randi([1,4]),...
    randi([0,1]),randi([0,1]),randi([0,1]),randi([0,1])];
    x0(ii,:)=[x_con x_dis];
end
% x0(1,:)=[9.4 0.201 0.1338 1345.5 25 0.6 500 0 0 24200 2 0 0 3 1 0 0 0 0]; %B737-800 typish input
%% -------------------------------------------------------------------------

%Sequence of the design variables

%Continuous
    %1. Aspect ratio AR (5 bits res 0.129)
    %2. Taper ratio TR (5 bits res 0.0129)
    %3. Thickness to chord ratio TCA (5 bits res 0.0026)
    %4. Wing area[ft^2] S (5 bits res=16.13)
    %5. Wing sweep at 25% chord[deg] Sweep (5 bits res 1.29)
    %6. Bypass ratio BPR (5 bits res 0.323)
    %7. Turbine entry temparature[R] TIT (5 bits 16.13)
    %8. Overall pressure ratio OPR (5 bits res=0.645)
    %9. Fan pressure ratio FPR (5 bits res 0.003226)
    %10. Thrust per engine[lbs] Thrust (5 bits res 322.581)
    bits_con=[5*ones(1,9) 10]; %Bits for continuous variable
    lb_con=[8 0.1 0.09 1000 0 0 0 0 0 20000]; %Lower bounds for the continuous design variables
    ub_con=[12 0.5 0.17 1500 40 10 500 20 0.1 30000]; %Upper bound for the continuous design variables
    
%Discrete 
    %A. Laminar flow technology (1-15) (4 bits eliminate 0000)
    %   1.No Laminar flow 2.NLF wing 3.NLF nascelle 4.NLF tail 5.NLF wing+nascelle
    %   6.NLF wing+tail 7.NLF tail+nascelle 8.HLFC wing+tail+nascelle 
    %   9.HLFC wing 10. HLFC nascelle 11. HLFC tail 12.HLFC wing+nascelle
    %   13..HLFC wing+tail 14.HLFC tail+nascelle 15.HLFC wing+tail+nascelle
    %   
    %B. Engine position (2-7) (2-7) (3 bits eliminate 0000, 1111)
    %   2. 2 on wing  3. 2 on fuselage(3) 4. 2 on wing,1 on fuselage
    %   5. 3 on fuselage(5) 6. 4 on wing(6) 7. 2 on wing, 2 on fuselage(7)
    %C. Engine type (1-4)
    %   1. DDF 2. GTF 3. CRF 4. Open rotors
    %D. Composite
    %   1. Wing: i.Yes(1) ii. No(0)
    %   2. Fuselage: i.Yes(1) ii. No(0)
    %   3. Nascelle: i.Yes(1) ii. No(0)
    %   4. Tail: i.Yes(1) ii. No(0)
    bits_dis=[2 1 1 3 2 1 1 1 1];
    lb_dis=[1 0 0 1 1 0 0 0 0];
    ub_dis=[4 1 1 8 4 1 1 1 1];

num_con=length(lb_con);
num_dis=length(lb_dis);
%--------------------------------------------------------------------------   
%% ------------------------------------------------------------------------    
ToTFunc=0;
 
pair=3; %Select the objective pair

%1. Wf VS TP fmin=[8000,100]; fmax=[100000,1000];
%2. NOX VS Wf fmin=[2, 8000]; fmax=[200, 100000];
%3. TP VS NOX fmin=[100,2]; fmax=[1000,200];

% pair=3; %Choose the objective pair
% fmin=[100,2]; fmax=[1000,200];

flag=1;
[xopt,ND,funeval,lgen,fmin_hist]= MOMDNP_bin_lam_NB_PAR(x0,ND0,xopt0,options,lb_con,ub_con,bits_con,lb_dis,ub_dis,bits_dis,flag,[],pair,Model)

ToTFunc=ToTFunc+funeval;
% Constraint Check
[len_x,len_y]=size(xopt);
for ii=1:len_x;
    x_con=xopt(ii,1:num_con);
    x_dis=xopt(ii,num_con+1:end);
    g(ii,:)=constraints(x_con,x_dis,Model);
    max_g(ii)=max(g(ii,:));
end
max_g

%Write results in text file
Filename='output3';
WriteResults(ToTFunc,xopt,ND,Filename)

%Sending an email upon completion
string='Lam Max_gen=50, pop_size=50, obj_pair=3';
sendemail(string)
Start_time
End_time=clock
Run_time = toc(Timestart)
matlabpool close
