clc;
clear all;
close all;
Timestart = tic;
Start_time=clock;
format longg;
matlabpool 8
options = goptions([]);
x0=ones(48,6); %population X variables
ND0=[];xopt0=[];xopt=[];
count=0;
fmin=[1e-6,1e-6];

%Continuous
bits_con=[8 8 8]; %Bits for continuous variable
lb_con=[1e-10 1e-10 1e-10]; %Lower bounds for the continuous design variables
ub_con=[5 5 5]; %Upper bound for the continuous design variables
    
%Discrete    
bits_dis=[2 2 2];
lb_dis=[1 1 1];
ub_dis=[4 4 4];

num_con=length(lb_con);
num_dis=length(lb_dis);

ToTFunc=0;
flag=1;

%objfunc file called in the hybridGA code itself
[xopt,ND,funeval,lgen,fmin_hist]= hybridGA_3truss(x0,ND0,xopt,options,lb_con,ub_con,bits_con,lb_dis,ub_dis,bits_dis,fmin,flag,count)

ToTFunc=ToTFunc+funeval;

% %Write results in text file
% Filename='output_3bar';
% WriteResults(ToTFunc,xopt,ND,Filename)
A = ToTFunc;
B = xopt;
C = ND;
fid=fopen('output_3bar.txt','w');
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

figure
plot(ND(:,1),ND(:,2),'k.');
xlabel('f_1')
ylabel('f_2')

matlabpool close