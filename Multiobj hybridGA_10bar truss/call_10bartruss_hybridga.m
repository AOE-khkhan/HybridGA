clc;
clear all;
close all;
Timestart = tic;
Start_time=clock;
format longg;
matlabpool 8
options = goptions([]);
x0=ones(48,20); %population X variables - CHECK FOR EVERY CHANGE IN DESIGN VARIABLES/POP. SIZE
ND0=[];xopt0=[];xopt=[];
count=0;
fmin=[1e-6,1e-6];


for i = 1:2:19
    vlb(i) = 1; 
    vlb(i+1) = 0.1;
    vub(i) = 4; 
    vub(i+1) = 40;
end
bits =[2 8 2 8 2 8 2 8 2 8 2 8 2 8 2 8 2 8 2 8];

i = 1; 
for counter1 = 1:2:19
   lb_dis(i) = vlb(counter1); %Discrete lower bound 
   lb_con(i) = vlb(counter1 + 1); %Continuous upper bound
   ub_dis(i) = vub(counter1); %Discrete upper bound 
   ub_con(i) = vub(counter1 + 1);  %Continuous upper bound
   i=i+1;
end

i = 1;
for counter2 = 1:2:19
    bits_dis(i) = bits(counter2); %Bits for discrete variables
    bits_con(i) = bits(counter2 + 1); %Bits for continuous variables
    i=i+1;
end

% %Continuous
% bits_con=[8 8 8 8 8 8 8 8 8 8]; %Bits for continuous variables
% lb_con = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %Lower bounds for the continuous
% ub_con=[40 40 40 40 40 40 40 40 40 40]; %Upper bound for the continuous
%     
% %Discrete    
% bits_dis=[2 2 2 2 2 2 2 2 2 2];
% lb_dis=[1 1 1 1 1 1 1 1 1 1];
% ub_dis=[4 4 4 4 4 4 4 4 4 4];

num_con=length(lb_con);
num_dis=length(lb_dis);

ToTFunc=0;
flag=1;

%objfunc file called in the hybridGA code itself
[xopt,ND,funeval,lgen,fmin_hist]= hybridGA_10bar(x0,ND0,xopt,options,lb_con,ub_con,bits_con,lb_dis,ub_dis,bits_dis,fmin,flag,count)

ToTFunc=ToTFunc+funeval;

% %Write results in text file
% Filename='output_3bar';
% WriteResults(ToTFunc,xopt,ND,Filename)
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