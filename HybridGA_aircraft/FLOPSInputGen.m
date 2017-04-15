% FLOPS Input File Generator
% function [] = FLOPSInputGen(x_con,mission,Inputs,Filename)
function [] = FLOPSInputGen(x_con,output,Inputs,Filename)
%% Currently this is a 737-8 aircraft ish aircraft setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wt_pax = 165;
bag_pax = 30;
NPF = round(0.07*Inputs.seats);
if mod(NPF,2) == 1;NPF = NPF-1;end
NPT = Inputs.seats - NPF;

titleline = ['Input Deck ',Filename];
fid = fopen([Filename '.in'],'w');

fprintf(fid,'%s\r\n', titleline);
%% -------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $OPTION ');
fprintf(fid,'%s\r\n','  IOPT=1, IANAL=3, ICOST=1,');
fprintf(fid,'%s\r\n',' $END ');
%% -------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $WTIN ');
fprintf(fid,'%s%6.2f%s\r\n','  DGW=', Inputs.GW,', ');
fprintf(fid,'%s%6.2f%s\r\n','  VMMO=', 0.82,', ');
fprintf(fid,'%s%6.2f%s\r\n','  DIH=', 6.0,', ');
fprintf(fid,'%s%6.2f%s\r\n','  HYDPR=', 3000.0,', ');
fprintf(fid,'%s%6.2f%s\r\n','  WPAINT=', 0.033,', ');
fprintf(fid,'%s%6.2f%s\r\n','  XL=', 129.5,', ');
fprintf(fid,'%s%6.2f%s\r\n','  WF=', 12.33,', ');
fprintf(fid,'%s%6.2f%s\r\n','  DF=', 13.5,', ');
fprintf(fid,'%s%6.2f%s\r\n','  XLP=', 98.5,', ');
fprintf(fid,'%s%6.2f%s\r\n','  SHT=', 353.1,', ');
fprintf(fid,'%s%6.2f%s\r\n','  SWPHT=', 284.2,', ');
fprintf(fid,'%s%6.2f%s\r\n','  ARHT=', 353.1,', ');
fprintf(fid,'%s%6.2f%s\r\n','  TRHT=', 0.281,', ');
fprintf(fid,'%s%6.2f%s\r\n','  TCHT=', 0.09,', ');
fprintf(fid,'%s%6.2f%s\r\n','  SVT=', 284.2,', ');
fprintf(fid,'%s%6.2f%s\r\n','  SWPVT=', 39.4,', ');
fprintf(fid,'%s%6.2f%s\r\n','  ARVT=', 1.24,', ');
fprintf(fid,'%s%6.2f%s\r\n','  TRVT=', 0.386,', ');
fprintf(fid,'%s%6.2f%s\r\n','  TCVT=', 0.09,', ');
% fprintf(fid,'%s%d%s\r\n','  NEW=', Inputs.NEW,', ');
fprintf(fid,'%s%d%s\r\n','  NEW=', output.NEW,', ');
fprintf(fid,'%s%d%s\r\n','  NEF=', output.NEF,', ');
fprintf(fid,'%s%d%s\r\n','  NPF=', NPF,', ');
fprintf(fid,'%s%d%s\r\n','  NPT=', NPT,', ');
fprintf(fid,'%s%d%s\r\n','  WPPASS=', wt_pax,', ');
fprintf(fid,'%s%d%s\r\n','  BPP=', bag_pax,', ');
fprintf(fid,'%s%d%s\r\n','  CARGOF=', 5500,', ');
% fprintf(fid,'%s%d%s\r\n','  NFLCR=', Inputs.NumCrew,', ');
% fprintf(fid,'%s%6.2f%s\r\n','  XMLG=', Inputs.X_mainLG,', ');
% fprintf(fid,'%s%6.2f%s\r\n','  XNLG=', Inputs.X_noseLG,', ');

fprintf(fid,'%s\r\n','  WSRV=1.8, ');
%fprintf(fid,'%s%1.0f%s\r\n','  FULWMX=', Inputs.wing_fuelcap,', ');
fprintf(fid,'%s\r\n','  IFUFU=1, ');
fprintf(fid,'%s\r\n','  MLDWT=0,  WAPU=1.0,  WHYD=1.0, ');

% Discrete technology modeling
fprintf(fid,'%s%2.1f%s\r\n','  FCOMP=', output.FCOMP,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FRHT=', output.FRHT,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FRVT=', output.FRVT,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FRFU=', output.FRFU,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FRNA=', output.FRNA,', ');
fprintf(fid,'%s%2.1f%s\r\n','  WAC=', output.WAC,', ');

% Engine technology modeling
fprintf(fid,'%s%6.2f%s\r\n','  WINL=', 0.0,', ');
fprintf(fid,'%s%6.2f%s\r\n','  WNOZ=', 0.0,', ');
fprintf(fid,'%s%6.2f%s\r\n','  WENG=', output.WENG,', ');

fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $CONFIN ');

fprintf(fid,'%s%7.2f%s\r\n','  DESRNG=', Inputs.DESRNG,', ');
fprintf(fid,'%s%7.2f%s\r\n','  GW=', Inputs.GW,', ');
fprintf(fid,'%s%7.2f%s\r\n','  AR=', x_con(1),', '); 
fprintf(fid,'%s%7.4f%s\r\n','  TR=', x_con(2),', ');
fprintf(fid,'%s%7.4f%s\r\n','  TCA=', x_con(3),', ');
fprintf(fid,'%s%7.2f%s\r\n','  SW=', x_con(4),', ');
fprintf(fid,'%s%5.2f%s\r\n','  SWEEP=', x_con(5),', ');
fprintf(fid,'%s%7.2f%s\r\n','  THRUST=', x_con(6),', '); 
fprintf(fid,'%s%7.3f%s\r\n','  VCMN=', 0.787,', '); 
fprintf(fid,'%s%7.2f%s\r\n','  CH=', 41000.0,', ');
fprintf(fid,'%s%7.2f%s\r\n','  HTVC=', 2.84,', ');
fprintf(fid,'%s%7.2f%s\r\n','  VTVC=', 0.243,', ');
fprintf(fid,'%s\r\n','  OFG=1., OFF=0., OFC=0.,');

fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $AERIN ');

fprintf(fid,'%s%7.5f%s\r\n','  VAPPR=', 142.0,',');
fprintf(fid,'%s%7.5f%s\r\n','  AITEK=', 1.819,', ');
fprintf(fid,'%s\r\n','  E=0.93365, ');
fprintf(fid,'%s%2.1f%s\r\n','  XLLAM=', output.XLLAM,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRUW=', output.TRUW,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRUH=', output.TRUH,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRLH=', output.TRLH,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRUV=', output.TRUV,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRLV=', output.TRLV,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRUN=', output.TRUN,', ');
fprintf(fid,'%s%2.1f%s\r\n','  TRLN=', output.TRLN,', ');

fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $COSTIN ');
fprintf(fid,'%s%d%s\r\n','  ROI=', 7.0,', ');
fprintf(fid,'%s%d%s\r\n','  FARE=', 0.10,', ');
fprintf(fid,'%s%d%s%s\r\n','  DEVST=', 2010,'.',', ');
fprintf(fid,'%s%d%s%s\r\n','  PLMQT=', 2015,'.',', ');
fprintf(fid,'%s%d%s\r\n','  DYEAR=', 2017,', ');
fprintf(fid,'%s%d%s%s\r\n','  TEMP=', 3321,'.',', ');
fprintf(fid,'%s%2.1f%s\r\n','  FUELPR=', 2.2,', ');
fprintf(fid,'%s%d%s\r\n','  NPOD=', output.NPOD,', ');
% fprintf(fid,'%s%f%s\r\n','  FMCOMP=', output.FMCOMP,', ');
% fprintf(fid,'%s%f%s\r\n','  FOCOMP=', output.FOCOMP,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FMWING=', output.FMWING,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FOWING=', output.FOWING,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FMNAC=', output.FMNAC,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FONAC=', output.FONAC,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FMTAIL=', output.FMTAIL,', ');
fprintf(fid,'%s%2.1f%s\r\n','  FOAC=', output.FOAC,', ');

fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $ENGDIN ');
fprintf(fid,'%s\r\n','  IDLE=1, IGENEN=1, NOX=1,');
fprintf(fid,'%s\r\n','  MAXCR=1, NGPRT=0, ');
fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $ENGINE ');

fprintf(fid,'%s\r\n','  IENG=2, IPRINT=0, ');

% Engine technology modeling
fprintf(fid,'%s%7.5f%s\r\n','  BPRDES=', output.BPRDES,', ');
fprintf(fid,'%s%7.5f%s\r\n','  TETDES=', output.TETDES,', ');
fprintf(fid,'%s%7.5f%s\r\n','  OPRDES=', output.OPRDES,', ');
fprintf(fid,'%s%7.5f%s\r\n','  FPRDES=', output.FPRDES,', ');
    
% fprintf(fid,'%s%7.5f%s\r\n','  OPRDES=', 29.5,', ');
% fprintf(fid,'%s%7.5f%s\r\n','  FPRDES=', 1.67,', ');
% fprintf(fid,'%s%7.5f%s\r\n','  TETDES=', 2660.0,', ');
    
fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n',' $MISSIN ');

fprintf(fid,'%s\r\n','  IFLAG=2, ');
fprintf(fid,'%s\r\n','  IRW=1, ');
fprintf(fid,'%s\r\n','  ITTFF=1, ');
fprintf(fid,'%s\r\n','  TAKOTM=0.4, ');
fprintf(fid,'%s\r\n','  TAXOTM=10, ');
fprintf(fid,'%s\r\n','  TAXITM=10, ');
fprintf(fid,'%s\r\n','  FWF=-1, ');
fprintf(fid,'%s\r\n','  THOLD=0.05, ');
fprintf(fid,'%s\r\n','  RESRFU=0.05, ');
fprintf(fid,'%s\r\n',' $END ');

%% ------------------------------------------------------------------------
fprintf(fid,'%s\r\n','START ');
fprintf(fid,'%s\r\n','CLIMB ');
fprintf(fid,'%s\r\n','CRUISE ');
fprintf(fid,'%s\r\n','DESCENT ');
fprintf(fid,'%s\r\n','END ');

%% ------------------------------------------------------------------------
%% Rerun with the sized aircraft for off-design mission runs
% for jj = 1:length(mission.route)
%     fprintf(fid,'%s\r\n','$RERUN ');
%     fprintf(fid,'%s\r\n','mywts = 1, wsr = 0., twr = 0. , ');
%     fprintf(fid,'%s%7.2f%s\r\n','  desrng=', mission.route(jj),', ');
%     payload = mission.pax(jj)*(wt_pax + bag_pax);
%     fprintf(fid,'%s%7.2f%s\r\n','  paylod=', payload,', ');
%     fprintf(fid,'%s\r\n','$END ');
%     
%     fprintf(fid,'%s\r\n',' $MISSIN ');
%     fprintf(fid,'%s\r\n','  IFLAG=0, ');
%     fprintf(fid,'%s\r\n',' $END ');
%     
%     fprintf(fid,'%s\r\n','START ');
%     fprintf(fid,'%s\r\n','CLIMB ');
%     fprintf(fid,'%s\r\n','CRUISE ');
%     fprintf(fid,'%s\r\n','DESCENT ');
%     fprintf(fid,'%s\r\n','END ');
% end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






