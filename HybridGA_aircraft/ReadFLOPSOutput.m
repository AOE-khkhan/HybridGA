% FLOPS Output File Reader
function [Outputs,nan_count,acdata_count] = ReadFLOPSOutput(filename)

TOC = [];BH = [];TD = [];LD = [];SSFOR = [];AMFOR = [];FUEL = [];NOX = [];FARE = [];
fname = [filename '.out'];
fid = fopen(fname);
count = 0;
nlines = 0;
nan_count = 0;
acdata_count = 0;
%Calculating the total number of lines
while ~feof(fid) && nan_count==0
    tline = fgetl(fid);
    nlines = nlines + 1;
    if ~isempty(tline)
                %% Read NOx data
        if length(tline)>36 && strcmp(tline(1:36),'     TOTAL NITROGEN OXIDES EMISSIONS')==1% && count == 0
            acdata_count = acdata_count + 1;
            try
                [~,~,~,~,~,NOX,~] = textread(fname,'%s%s%s%s%s%f%s',1,'headerlines',nlines-1);
                %count = 1;
                if isnan(NOX) == 1
                    nan_count = nan_count + 1;
					break
                end
            catch
                nan_count=nan_count+1;
                break
            end
        end
        %% Read cost data and calculate Block hours
        %nan_count=0;
        if length(tline)>22 && strcmp(tline(1:23),' DIRECT OPERATING COSTS')==1
            acdata_count = acdata_count + 1;
            try
                [~,~,~,~,DOC] = textread(fname,'%s%s%s%s%f',1,'headerlines',nlines);
                [~,~,~,~,~,DOCperBH] = textread(fname,'%s%s%s%s%s%f',1,'headerlines',nlines+1);
                [~,~,~,~,~,IOC] = textread(fname,'%s%s%s%s%s%f',1,'headerlines',nlines+25);
                TOC = [TOC;DOC + IOC];
                BH = [BH;DOC/DOCperBH];
                if isnan(DOC) == 1 || isnan(IOC)==1 || isnan(DOCperBH)==1
                    nan_count = nan_count + 1;
					break
                end
            catch
                nan_count=nan_count+1;
                break
            end
        end
%         %% Read ticket price data
%         %nan_count=0;
%         if length(tline)>29 && strcmp(tline(1:29),' FOR AN ROI OF  7.000 PERCENT')==1% && count == 0
%             acdata_count = acdata_count + 1;
%             try
%                 [~,~,~,~,~,~,~,~,FARE,~] = textread(fname,'%s%s%s%s%f%s%s%s%f%s',1,'headerlines',nlines-1);
%                 if isnan(FARE)==1
% 					nan_count=nan_count+1;
% 					break
% 				end
%             catch
%                 nan_count=nan_count+1;
%                 break
%             end
%         end
        %% Read take-off, landing data
        %nan_count=0;
        if length(tline)>15 && strcmp(tline(1:15),'#OBJ/VAR/CONSTR')==1 && count == 0
            acdata_count = acdata_count + 1;
            try
                [FUEL,~,~,TD,LD,AMFOR,SSFOR] = textread(fname,'%f%f%f%f%f%f%f',1,'headerlines',nlines+2);
                count = 1;
				if isnan(FUEL) == 1 || isnan(TD)==1 || isnan(LD) == 1
					nan_count=nan_count+1;
					break
				end
            catch
                nan_count=nan_count+1;
                break
            end
        end
    end
end
fclose(fid);
Outputs.TOC = TOC;
Outputs.BH = BH;
Outputs.TD = TD;
Outputs.LD = LD;
Outputs.SS = SSFOR;
Outputs.tclimb = AMFOR;
Outputs.FUEL = FUEL;
Outputs.NOX = NOX;
Outputs.FARE = FARE;
fname1 = [filename,'.in'];
% delete(fname1)
% delete(fname)