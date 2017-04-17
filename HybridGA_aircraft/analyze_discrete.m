%% Analyzing discrete technologies
function [output] = analyze_discrete(x_con,x_dis)

% Inputs: x_dis
% 1. Composites: x_dis(1:4)
%    1 - 0/1 for wing
%    2 - 0/1 for fuselage
%    3 - 0/1 for nacelle
%    4 - 0/1 for tail

% 2. Position and number of engine/s x_dis(5)

%    1 - 2 wing 
%    2 - 2 fuselage 
%    3 - 2 wing + 1 fuselage
%    4 - 3 fuselage
%    5 - 4 wing 
%    6 - 2 wing + 2 fuselage
%    7 - 1 fuselage
%    8 - 4 wing + 1 fuselage
% 3. Laminar flow technologies x_des(6)

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

% 4. Engine Technologies: x_dis(7)
%    1 - DDF
%    2 - GTF
%    3 - CRTF
%    4 - OR
%% Composites
xcomp = x_dis(1:4);
FCOMP = 0.0; FRFU = 1.0; FRNA = 1.0; FRHT = 1.0; FRVT = 1.0;
if xcomp(1) == 1
    FCOMP = 1.0;
end
if xcomp(2) == 1
    FRFU = 0.82;
end
if xcomp(3) == 1
    FRNA = 0.8;
end
if xcomp(4) == 1
    FRHT = 0.75;
    FRVT = 0.75;
end

%% Engine position and number
xeng = x_dis(5);
switch xeng
    case 1
        NEW = 2;
        NEF = 0;
    case 2
        NEW = 0;
        NEF = 2;
    case 3
        NEW = 2;
        NEF = 1;
    case 4
        NEW = 0;
        NEF = 3;
    case 5
        NEW = 4;
        NEF = 0;
    case 6
        NEW = 2;
        NEF = 2;
    case 7
        NEW = 0;
        NEF = 1;
    case 8
        NEW = 4;
        NEF = 1;
end
   
%% Laminar flow technologies
xltech = x_dis(6);
XLLAM = 1.0;
TRUW = 0.0;TRLW =0.0;
FMWING = 1.0;FOWING = 1.0;WAC = 1.0;
TRUN = 0.0;TRLN = 0.0;
FMNAC = 1.0;FONAC = 1.0;
TRUH = 0.0;TRLH = 0.0;
TRUV = 0.0;TRLV = 0.0;
FMTAIL = 1.0;FOAC = 1.0;

switch xltech
    case 1
        XLLAM = 0.0;
    case 2 %NLF wing
        flag = 1;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.02;
        FOWING = 1.02;
    case 3 %HLFC wing
        flag = 2; 
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.5;
        FOWING = 1.5;
        WAC = 2.5;
    case 4 %HLFC nacelle
        TRUN = 60;
        TRLN = 60;
        FMNAC = 1.5;
        FONAC = 1.5;
        WAC = 1.1667;
    case 5 %HLFC tail
        TRUH = 50; TRLH = 50;
        TRUV = 50; TRLV = 50;
        FMTAIL = 1.5;
        FOAC = 1.5;
        WAC = 1.2;
    case 6
        flag = 2;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.5;
        FOWING = 1.5;
        WAC = 2.6667;
        TRUN = 60;
        TRLN = 60;
        FMNAC = 1.5;
        FONAC = 1.5;
    case 7
        flag = 2;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.5;
        FOWING = 1.5;
        WAC = 2.7;
        TRUH = 50; TRLH = 50;
        TRUV = 50; TRLV = 50;
        FMTAIL = 1.5;
        FOAC = 1.5;
    case 8
        TRUN = 60;
        TRLN = 60;
        FMNAC = 1.5;
        FONAC = 1.5;
        WAC = 1.3667;
        TRUH = 50; TRLH = 50;
        TRUV = 50; TRLV = 50;
        FMTAIL = 1.5;
        FOAC = 1.5;
    case 9
        flag = 2;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.5;
        FOWING = 1.5;
        WAC = 2.8667;
        TRUH = 50; TRLH = 50;
        TRUV = 50; TRLV = 50;
        FMTAIL = 1.5;
        FOAC = 1.5;
        TRUN = 60;
        TRLN = 60;
        FMNAC = 1.5;
        FONAC = 1.5;
    case 10
        flag = 1;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.02;
        FOWING = 1.02; 
        TRUH = 50; TRLH = 50;
        TRUV = 50; TRLV = 50;
        FMTAIL = 1.5;
        FOAC = 1.5;
        WAC = 1.2;
    case 11
        flag = 1;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.02;
        FOWING = 1.02;
        TRUN = 60;
        TRLN = 60;
        FMNAC = 1.5;
        FONAC = 1.5;
        WAC = 1.1667;
    case 12
        flag = 1;
        [TRUW] = laminarflow_code(x_con,flag);
        FMWING = 1.02;
        FOWING = 1.02;
        TRUN = 60;
        TRLN = 60;
        FMNAC = 1.5;
        FONAC = 1.5;
        WAC = 1.3667;
        TRUH = 50; TRLH = 50;
        TRUV = 50; TRLV = 50;
        FMTAIL = 1.5;
        FOAC = 1.5;
end

%% Engine Technologies
engtech = x_dis(7);
switch engtech
    case 1 %DDF
        BPRDES = 5 + x_con(7)/2;
        TETDES = 3010 + x_con(8);
        OPRDES = 35 + x_con(9);
        FPRDES = 1.6 + x_con(10);
        WENG = 0.80;
    case 2 %GTF
        BPRDES = 10 + x_con(7);
        TETDES = 3010 + x_con(8);
        OPRDES = 30 + x_con(9);
        FPRDES = 1.5 + x_con(10);
        WENG = 0.85;
    case 3 %CRTF/CRR
        BPRDES = 15 + x_con(7);
        TETDES = 3010 + x_con(8);
        OPRDES = 30 + x_con(9);
        FPRDES = 1.4 + x_con(10);
        WENG = 0.90;
    case 4 %OR
        BPRDES = 25 + x_con(7);
        TETDES = 3010 + x_con(8);
        OPRDES = 25 + x_con(9);
        FPRDES = 1.3 + x_con(10);
        WENG = 0.75;
end
%% Collect the information and store it in an object
output.FCOMP = FCOMP;output.FRFU = FRFU;
output.FRNA = FRNA;output.FRHT = FRHT;output.FRVT = FRVT;

output.NEW = NEW;output.NEF = NEF; output.NPOD = NEW+NEF;

output.XLLAM = XLLAM;
output.TRUW = TRUW;output.TRLW =TRLW; 

output.FMWING = FMWING;output.FOWING = FOWING;output.WAC = WAC;
output.TRUN = TRUN;output.TRLN = TRLN;
output.FMNAC = FMNAC;output.FONAC = FONAC;
output.TRUH = TRUH;output.TRLH = TRLH;
output.TRUV = TRUV;output.TRLV = TRLV;
output.FMTAIL = FMTAIL;output.FOAC = FOAC;


=======
output.BPRDES = BPRDES; output.TETDES = TETDES;
output.OPRDES = OPRDES; output.FPRDES = FPRDES;
output.WENG = WENG;
>>>>>>> d756537166b88afba233ed14c3b962140d28f3f8

function [TRUW] = laminarflow_code(x_con,flag)

if flag == 1% NLF
    S = 28e6;
elseif flag ==2 % HLFC
    S = 50e6; 
end
swp_25 = x_con(5);

%rho=0.2847;
mu=1.12e-5;
V=0.787*316.3742;
del_beta=2; %degree

b=sqrt(x_con(1)*x_con(4));
c_root=(2*x_con(4)/(b*(1+x_con(2))));
c_tip=x_con(2)*c_root;
beta_LE=tand(swp_25)+((1-x_con(2))/(x_con(1)*(1+x_con(2))));

%Re1=-((180/pi*atan((180/pi*atan(slp_LE))-13)/(13*1.6))*1.6)
Re=-atan((beta_LE + del_beta - 13)/(13*1.6))*1.6 + S;

li=(Re*mu/V)/0.3048; %in feet

l_lam = zeros(51,1);
l_lam(1) = min(li, c_root/2);

Aw=0;
res=(b/2)/50;
for i=2:51
    xi=(i-1)*res;
    y_LE=-beta_LE*xi;
    y_TE=((2/b)*(c_root-c_tip-(b/2*beta_LE))*xi)-c_root;
    ci=abs(y_LE-y_TE);
    l_lam(i) = min(ci/2,li);
    Ai=0.5*(l_lam(i-1)+l_lam(i))*res;
    Aw=Aw+Ai;
end

TRUW=2*Aw/x_con(4)*100;