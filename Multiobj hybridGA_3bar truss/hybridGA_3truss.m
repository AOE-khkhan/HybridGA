function [xopt,ND,funeval,lgen,fmin_hist] = hybridGA_3truss(x0,ND,xopt...
    ,options,lb_con,ub_con,bits_con,lb_dis,ub_dis,bits_dis,fmin,flag,count,~,~,~,~,~,~,~,~,~,~,~,~)
%GENETIC minimizes a fitness function using a simple genetic algorithm.
%
%	X=GENETIC('FUN',X0,OPTIONS,VLB,VUB) uses a simple (haploid) 
%       genetic algorithm to find a minimum of the fitness function 
%       FUN.  FUN can be a user-defined M-file: FUN.M, or it can be a 
%	string containing the function itself.  The user may define all
%       or part of an initial population X0. Any undefined individuals 
%	will be randomly generated between the lower and upper bounds
%	(VLB and VUB).  If X0 is an empty matrix, the entire initial
%	population will be randomly generated.  Use OPTIONS to specify 
%	flags, tolerances, and input parameters.  Type HELP GOPTIONS
%       for more information and default values.
%
%	X=GENETIC('FUN',X0,OPTIONS,VLB,VUB,BITS) allows the user to 
%	define the number of BITS used to code non-binary parameters
%	as binary strings.  Note: length(BITS) must equal length(VLB)
%	and length(VUB).  If BITS is not specified, as in the previous 
%	call, the algorithm assumes that the fitness function is 
%	operating on a binary population.
%
%	X=GENETIC('FUN',X0,OPTIONS,VLB,VUB,BITS,P1,P2,...) allows up 
%	to ten arguments, P1,P2,... to be passed directly to FUN.
%	F=FUN(X,P1,P2,...). If P1,P2,... are not defined, F=FUN(X).
%
%	[X,FOPT,STATS,NFIT,FGEN,LGEN,LFIT]=GENETIC(<ARGS>)
%          X       - design variables of best ever individual
%          FOPT    - fitness value of best ever individual
%          STATS   - [min mean max] fitness values for each generation
%          NFIT	 - number of fitness function evalations
%          FGEN    - first generation population
%          LGEN    - last generation population
%          LFIT    - last generation fitness
%
%       The algorithm implemented here is based on the book: Genetic
%       Algorithms in Search, Optimization, and Machine Learning,
%       David E. Goldberg, Addison-Wiley Publishing Company, Inc.,
%       1989.
%
%	Originally created on 1/10/93 by Andrew Potvin, Mathworks, Inc. 
%	Modified on 2/3/96 by Joel Grasmeyer.
%   Modified on 11/12/02 by Bill Crossley.
%   Modified on 7/20/04 by Bill Crossley.
%   Modified on 01/24/10 by Satadru Roy

%Combining the design variables (Continuous+Discrete)
%19 design variables
%   10 continuous
%   09 discrete
num_con=length(lb_con);
num_dis=length(lb_dis);
num_variables=num_con+num_dis;
%Combined continuous and discrete variables
bits=cat(2,bits_con,bits_dis);
vlb=cat(2,lb_con,lb_dis);
vub=cat(2,ub_con,ub_dis);
%% ------------------------------------------------------------------------
% Load input arguments and check for errors
if nargin<4,
    error('No population bounds given.')
    elseif (size(vlb,1)~=1) || (size(vub,1)~=1),
        % Remark: this will change if algorithm accomodates matrix variables
        error('VLB and VUB must be row vectors')
    elseif (size(vlb,2)~=size(vub,2)),
        error('VLB and VUB must have the same number of columns.')
    elseif (size(vub,2)~=size(x0,2)) && (size(x0,1)>0),
        error('X0 must all have the same number of columns as VLB and VUB.')
    elseif any(vlb>vub),
        error('Some lower bounds greater than upper bounds')
else
        x0_row = size(x0,1);
    for i=1:x0_row,
        if any(x0(x0_row,:)<vlb) || any(x0(x0_row,:)>vub),
            error('Some initial population not within bounds.')
        end % if initial pop not within bounds
    end % for initial pop
end % if nargin<4   

if nargin<6,
    bits = [];
    elseif (size(bits,1)~=1) || (size(bits,2)~=size(vlb,2)),
        % Remark: this will change if algorithm accomodates matrix variables
        error('BITS must have one row and length(VLB) columns')
    elseif any(bits~=round(bits)) || any(bits<1),
        error('BITS must be a vector of integers >0')
end % if nargin<6

% Determine all options
% Remark: add another options index for type of termination criterion
if size(options,1)>1,
    error('OPTIONS must be a row vector')
else
    % Use default options for those that were not passed in
    options = goptions(options);
end
PRINTING = options(1);
BSA = options(2);


% Since operators are tournament selection and uniform crossover and
% default coding is Gray / binary, set crossover rate to 0.50 and use
% population size and mutation rate based on Williams, E. A., and Crossley,
% W. A., "Empirically-derived population size and mutation rate guidelines
% for a genetic algorithm with uniform crossover," Soft Computing in
% Engineering Design and Manufacturing, 1998.  If user has entered values
% for these options, then user input values are used.
if options(11) == 0,
%     pop_size=4*sum(bits);
    n=6;
    pop_size = 8*n;
else
    pop_size = options(11);
end
if options(12) == 0,
    Pc = 0.5;
else
    Pc = options(12);
end
if options(13) == 0,
%     Pm = (sum(bits) + 1) / (2 * pop_size * sum(bits));
    Pm=0.005;
else
    Pm = options(13);
end
max_gen = options(14);
% Ensure valid options: e.q. Pc,Pm,pop_size,max_gen>0, Pc,Pm<1
if any([Pc Pm pop_size max_gen]<0) | any([Pc Pm]>1),
    error('Some Pc,Pm,pop_size,max_gen<0 or Pc,Pm>1')
end

% Encode fitness (cost) function if necessary
ENCODED = any(any(([vlb; vub; x0]~=0) & ([vlb; vub; x0]~=1))) |  ....
    ~isempty(bits);
if ENCODED,
    [fgen,lchrom] = encode(x0,vlb,vub,bits);
else
    fgen = x0;
    lchrom = size(vlb,2);
end

% Display warning if initial population size is odd
if rem(pop_size,2)==1,
    disp('Warning: Population size should be even.  Adding 1 to population.')
    pop_size = pop_size +1;
end

% Form random initial population if not enough supplied by user
if size(fgen,1)<pop_size,
    fgen = [fgen; (rand(pop_size-size(fgen,1),lchrom)<0.5)];
end

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
new_gen = fgen;
x_des=zeros(pop_size,size(vlb,2)); %%%%rows=pop_size;columns=no of lower bound=no of design variables 
funeval=0;
ter_crit=0;
ALL_obj_val=zeros(pop_size,2,max_gen+1);
ALL_x_des=zeros(pop_size,size(vlb,2),max_gen+1); %%%%x_des matrix for every gen
termi_crit=zeros(max_gen+1,1);
cons_check=10; %Termination criteria, number of gen to observe for no change
fmin_hist=zeros(max_gen+1,2);
frac_uto=0.9; %Fractional offset for the utopia point

%Goal attainment percentage
ck=ones(1,8*n);
dk=ones(1,8*n);

%Model1=Model;
% ND=[];
% xopt=[];
%x_con = x0(ii,1:num_con); x_dis = x0(ii,num_con+1:end);
%objective values for initial population
obj_val=zeros(pop_size,2);
for ii=1:pop_size                 %%%%initial obj val-both functions & continuous & discrete
    obj_val(ii,:)=objfunc_initial_3truss(x0(ii,1:num_variables),flag,count); 
    funeval=funeval+1;
end
% obj_val=[(fmin(1)+(fmax(1)-fmin(1))*rand(pop_size,1)),(fmin(2)+(fmax(2)-fmin(2))*rand(pop_size,1))];
fmin=frac_uto*[min(obj_val(:,1)),min(obj_val(:,2))]; %%%%obj_val(:,1) for function 1 & obj_val(:,2) for function 2
% Header display
if PRINTING>=1,
    if ENCODED,
        disp('Variable coding as binary chromosomes successful.')
        disp('')
        fgen = decode(fgen,vlb,vub,bits);
    end
end
    
% Set up main loop
for generation = 1:max_gen+1,
    old_gen = new_gen;
    % Decode binary strings if necessary
    if ENCODED,
        x_pop = decode(old_gen,vlb,vub,bits);
    else
        x_pop = old_gen;
    end
    fmin_Cgen=frac_uto*[min(obj_val(:,1)),min(obj_val(:,2))]; %%%%calculating func values to be optimized to(corner point in plot)                                                                                                        
    if fmin_Cgen(1)<fmin(1)                                                                                                                                          
        fmin(1,1)=fmin_Cgen(1); %%%%dynamic min point - min value replaced                                                                                                                                   
    end                                                                                                                                                              
    if fmin_Cgen(2)<fmin(2)                                                                                                                                          
        fmin(1,2)=fmin_Cgen(2); %%%%dynamic min point - min value replaced                                                                                                                                     
    end     
    fmin_hist(generation,:)=fmin %history of the fmin
    %Work by Satadru Roy
    %Hybrid Optimization-Multi Search Method
    %Developed by Satadru Roy on January 24,2010
    % as per the guidance found in the PhD dissertation of W.Hart UCSD,1994
    %Uses the technique of assigning SQP output to the closest GA fitness
    %value
    if generation>1  %%%%if generation=1, no goal point?
        coarse = (vub-vlb)./((2.^bits)-1);
        X_pop=x_pop;
        XX_pop=x_pop;
        Obj_val=obj_val;
        X_des=x_des;
        Goal=goal;
        cK=ck;
        dK=dk;
        parfor i=1:pop_size
%             Model=Model1;
%             t_st=tic;
            sprintf('%s%d%s%d','Generation:',generation,'  Individual:',i)
            %Keeping discrete variable constant for the SQP process
            xx_pop=XX_pop(i,1:num_con);
            xx_pop1=XX_pop(i,num_con+1:end);
            x0_con=xx_pop; %Initials for SQP from GA (Design variables for SQP)
            x_dis=xx_pop1; %Discrete variable sent as parameters. 
            aa=Obj_val(i,:);%;bb=Obj_val(i,2); %%%%obj values of ith gen for both functions 1&2 
%             goal=goalpointgen(fmin(1),fmax(2),fmax(1),fmin(2),aa(1),aa(2));
            gg=Goal(i,:); %%%%goal values of ith gen for both functions 1&2
            %New addition 12/11/11
            w1=cK(i)*abs(gg(1));
            w2=dK(i)*abs(gg(2));
            weight=[w1 w2];
            options=optimset('LargeScale','off','GradObj','off','Algorithm','active-set','MaxFunEvals',300,'UseParallel','always');
            [x_con,fval,~,exitflag,output]=fgoalattain(@(x_con) objfunc_3truss(x_con,x_dis,flag,count),x0_con,gg,weight,[],[],[],[],lb_con,ub_con,@(x_con) constraints_3truss(x_con,x_dis),options);
            funeval=funeval+output.funcCount;
            g_check=constraints_3truss(x_con,x_dis); %%%%check? should be constraints...
            if (exitflag>=0 && isempty(find(g_check>1e6))==1)
                Obj_val(i,:)=fval;  %%%%func value obtained from SQP set as next objective value 
                X_des(i,:)=[x_con,x_dis];
                %Applying Lamarckism 
                pos=1;    
                aaa=X_des(i,:);
                bbb=X_pop(i,:);
                ccc=old_gen(i,:);
                for j=1:num_con 
                    graycode=GCdevelop(bits(j));
                    for k=2:(2^bits(j))
                         fit_val_U=vlb(j)+((k-1)*coarse(j));
                         fit_val_L=vlb(j)+((k-2)*coarse(j));
                         if (aaa(j)<=fit_val_U && aaa(j)>=fit_val_L)
                             if ((aaa(j)-fit_val_L)/(fit_val_U-aaa(j))<=1)
                                 bbb(j)=fit_val_L;
                                 ccc(pos:pos+bits(j)-1)=graycode(k-1,2:end);
                                 pos=pos+bits(j);
                                 break;
                             else
                                 bbb(j)=fit_val_U;
                                 ccc(pos:pos+bits(j)-1)=graycode(k,2:end);
                                 pos=pos+bits(j);
                                 break;
                             end
                          end
                     end
                end
                X_des(i,:)=aaa;
                X_pop(i,:)=bbb;
                old_gen(i,:)=ccc;
                
            else
                Obj_val(i,:)=1e15; %%%%exitflag~=0; apply high value
            end
%             t_end=toc(t_st);
%             time_left=(((max_gen+1)*pop_size*t_end)-(generation*i*t_end))/3600;
%             sprintf('%s%f','Estimated time to complete (hrs): ',time_left)
        end
        obj_val=Obj_val;
        x_des=X_des;
        x_pop=X_pop;
        disp('Non-Dominated set')
        [ND,xopt,NC_count]=paretogen(ND,obj_val,xopt,x_des)
        ALL_obj_val(:,:,generation)=obj_val;
        ALL_x_des(:,:,generation)=x_des;
        obj_val
        termi_crit(generation)=NC_count;
    end
    %end hybrid optimization
    %%%%%%%%%% 
%   x_last=x_des(:,:,generation);        
% 
       
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
% Check for termination
% The default termination criterion is bit string affinity.  Also
% available are fitness tolerance across five generations and number of
% consecutive generations with same best fitness.  These can be used
% concurrently.

% if BSA > 0,  % bit string affinity criterion
%     if generation > 1
%         bitlocavg = mean(old_gen,1);
%         BSA_pop = 2 * mean(abs(bitlocavg - 0.5));
%         if BSA_pop >= BSA,
%             if PRINTING >=1
%                 fprintf('\n')
%                 disp('GA stopped based on bit string affinity value.')
%             end
%             lfit = obj_val;
%             if ENCODED,
%                 lgen = x_pop;
%             else
%                 lgen = old_gen;
%             end
%             return
%         end
%     end
% end
%    
 % Always save last generation.  This allows user to cancel and
 % 
  % restart with x0 = lgen
  % 
    lfit=obj_val;                                                                                                                                                        
    if ENCODED,                                                                                                         
        lgen = x_pop;                                                                                                                         
    else                                                                              
        lgen = old_gen;                                                              
    end                                                                                                                                         
                                                                                                                                                    
    save('GenerationHistory.mat','lgen','ALL_obj_val','ALL_x_des')                                                                       
                                                                                                      
                                                                                                                                                            
  % Check for termination
  % %                                                     
  % Termination Criteria: If no chnage in the set of non-dominated designs
  % are      %                                                            
  % found in the last "cons_check" generation
  % %                                                     
    if generation>cons_check                                                                                                                                    
        if (sum(termi_crit(generation-cons_check:generation))==cons_check)                                                                        
            fprintf('%s%d%s\n','Termination criteria satisfied: No change in the non dominated design') %                                                                         
            ter_crit=1;                                                                                                                                   
            break;                                                                                                                                               
        end                                                                                                                                                  
    end                                                                                                                                                         
        
       
    % Mutation
    new_gen = mutate(old_gen,Pm);
    % Tournament selection
    [new_gen, obj_val] = Nbranchtourney_new(new_gen,obj_val);
    
    %Generating the new goal points
    goal=goalpointgen_new(obj_val,fmin);
    % Crossover
    new_gen = uniformx(new_gen,Pc);
    
    % Always save last generation.  This allows user to cancel and
    % restart with x0 = lgen
    lfit=obj_val;
    if ENCODED,
        lgen = x_pop;
    else
        lgen = old_gen;
    end
    
    save LastGen.mat lgen    
    save OBJ_val.mat ALL_obj_val
    save X_DES.mat ALL_x_des
    
end % for max_gen

% Maximum number of generations reached without termination
lfit = obj_val;
if ENCODED,
    lgen = x_pop;
else
    lgen = old_gen;
end
if ter_crit==0,
    fprintf('\n')
    disp('Maximum number of generations reached!')
end


% end genetic


%% ------------------------------------------------------------------------
%Work by Satadru Roy February 20,2010
function [gen,lchrom,coarse] = encode(x,vlb,vub,bits)
lchrom = sum(bits);
coarse = (vub-vlb)./((2.^bits)-1);
[npop,x_col]=size(x);
gen=zeros(npop,lchrom);
for i=1:npop
    for j=1:x_col
        sl(i,j)=round(((x(i,j)-vlb(j))/coarse(j))+1);
    end
end
      
for i=1:npop
    pos=1;
    for j=1:x_col
        %Reading Graycode from the Graycode developer%  
        data_greycode=GCdevelop(bits(j));
        for k=1:bits(j)
            gen(i,pos)=data_greycode(sl(i,j),k+1);
            pos=pos+1;
        end   
    end
end

% end encode
%%%%%%%%%%

% Work by Satadru Roy February 20,2010
function [x,coarse]=decode(gen,vlb,vub,bits) 
coarse = (vub-vlb)./((2.^bits)-1);
[x_rowgen,x_colgen]=size(gen);
x=zeros(x_rowgen,length(bits));
for i=1:x_rowgen
    pos=1;
    for j=1:length(bits)
        %Reading Graycode from the Graycode developer%  
        data_greycode=GCdevelop(bits(j));
        for k=1:(2^bits(j))
            if (data_greycode(k,2:(bits(j)+1))==gen(i,pos:(pos+bits(j)-1)))
                pos=pos+bits(j);
                x(i,j)=vlb(j)+((data_greycode(k,1)-1)*coarse(j));
                break;
            end
        end
    end
end
%end decode
%%%%%%%%%%


function [new_gen,mutated] = mutate(old_gen,Pm)
%MUTATE Changes a gene of the OLD_GEN with probability Pm.
%	[NEW_GEN,MUTATED] = MUTATE(OLD_GEN,Pm) performs random
%       mutation on the population OLD_POP.  Each gene of each
%       individual of the population can mutate independently
%       with probability Pm.  Genes are assumed possess boolean
%       alleles.  MUTATED contains the indices of the mutated genes.
%
%	Copyright (c) 1993 by the MathWorks, Inc.
%	Andrew Potvin 1-10-93.

mutated = find(rand(size(old_gen))<Pm);
new_gen = old_gen;
new_gen(mutated) = 1-old_gen(mutated);

% end mutate

% Functions performs the N-Branch tournament selection as per the paper "Using
% the Two branch tournament Genetic Algorithm for multiobjective
% design"-W.A.Crossley, Andrea M.Cook, David W. Fanjoy and Vipperla
% Venkayya
%With added features to improve the spread of the Pareto frontier
%(W.A.Crossley, Satadru Roy)-version 3

function [new_gen, new_obj_val]=Nbranchtourney_new(old_gen,obj_val)

% [old_gen,obj_val]
% Initialize nselected vector and indices of old_gen
new_gen = [];
%nselected = zeros(size(old_gen,1),1);
%i_old_gen = 1:size(old_gen,1);

% Perform two "tournaments" to generate size(old_gen,1) new individuals

n=size(old_gen,1)/8; % pop_size=8*n
pot_I=zeros(4*n,size(old_gen,2),2);
obj_pot_I=zeros(4*n,size(obj_val,2),2);

% Selection phase I: Creating Level I Pots
for j = 1:2,
    count=1;
    % Shuffle the old generation and the corresponding fitness values
    [old_gen,i_shuffled] = shuffle(old_gen);
    obj_val = obj_val(i_shuffled,:);
    
    % Keep the best of each pair of individuals
    for index = 1:2:(size(old_gen,1)-1);
        [~,i_min] = min([obj_val(index,j);obj_val(index+1,j)]); %%%%finding min of 2 chosen individuals
        %selected = i_min + [0:2:size(old_gen,1)-2];
        pot_I(count,:,j) = old_gen((index-1)+i_min,:); %%%%adding individual to pool1 & then same process repeated for pool2 as j=1:2
        obj_pot_I(count,:,j)=obj_val((index-1)+i_min,:); %%%%adding obj value of above individual
        count=count+1;
    end
end
% [pot_I,obj_pot_I]

%Selection phase II: Creating Level II Pots
for j=1:2
    [pot_I(:,:,j),i_shuffled]=shuffle(pot_I(:,:,j));
    obj_pot_I(:,:,j)= obj_pot_I(i_shuffled,:,j);
end

pot_II(:,:,1)=pot_I(1:2*n,:,1); %%%%subpool1 - 2n from pool1 - phi1-phi1 strong
pot_II(:,:,2)=pot_I(1:2*n,:,2); %%%%subpool2 - 2n from pool2 - phi2-phi2 strong

obj_pot_II(:,:,1)=obj_pot_I(1:2*n,:,1); 
obj_pot_II(:,:,2)=obj_pot_I(1:2*n,:,2);

count=1;
for j=(2*n)+1:(4*n)
    pot_II(count,:,3)=pot_I(j,:,1); %%%%subpool3 - 2n from pool1 - phi1 strong 
    pot_II(count+1,:,3)=pot_I(j,:,2); %%%%subpool3 - 2n from pool2 - phi2 strong
    %%%%subpool 3 - phi1-phi2 strong
    obj_pot_II(count,:,3)=obj_pot_I(j,:,1);
    obj_pot_II(count+1,:,3)=obj_pot_I(j,:,2);
    
    count=count+2;
end
new_gen=[pot_II(1:2*n,:,1);pot_II(1:2*n,:,2);pot_II(1:count-1,:,3)]; %%%%final subpools 1,2,3

new_obj_val=[obj_pot_II(1:2*n,:,1);obj_pot_II(1:2*n,:,2);obj_pot_II(1:count-1,:,3)];
% [new_gen,new_obj_val]

% end tourney


function [new_gen,index] = shuffle(old_gen)
%SHUFFLE Randomly reorders OLD_GEN into NEW_GEN.
%	 [NEW_GEN,INDEX] = MATE(OLD_GEN) performs random reordering
%        on the indices of OLD_GEN to create NEW_GEN.
%	 INDEX is a vector containing the shuffled row indices of OLD_GEN.
%
%	 Created on 1/21/96 by Joel Grasmeyer

[~,index] = sort(rand(size(old_gen,1),1));
new_gen = old_gen(index,:);

% end shuffle


function [new_gen,sites] = uniformx(old_gen,Pc)
%UNIFORMX Creates a NEW_GEN from OLD_GEN using uniform crossover.
%	  [NEW_GEN,SITES] = UNIFORMX(OLD_GEN,Pc) performs uniform crossover
%         on consecutive pairs of OLD_GEN with probability Pc.
%	  SITES shows which bits experienced crossover.  1 indicates
%	  allele exchange, 0 indicates no allele exchange.  SITES has
%	  size(old_gen,1)/2 rows.
%
%  	  Created 1/20/96 by Joel Grasmeyer

new_gen = old_gen;
sites = rand(size(old_gen,1)/2,size(old_gen,2)) < Pc;
for i = 1:size(sites,1),
  new_gen([2*i-1 2*i],find(sites(i,:))) = old_gen([2*i 
2*i-1],find(sites(i,:)));

% %Genetically tailored new generation to account for discrete technologies for aircraft sizing
% forbid=[0 0 1];
% %Supress output 0 in Laminar flow tech
% if (sum(new_gen(2*i-1,56:57))==0)
%     position=randi([1,2]);
%     new_gen(2*i-1,55+position)=abs(new_gen(2*i-1,55+position)-1);
% end
% if (sum(new_gen(2*i,56:57))==0)
%     position=randi([1,2]);
%     new_gen(2*i,55+position)=abs(new_gen(2*i,55+position)-1);
% end
% 
% %Supress output 1 & 2 in engine position
% if (sum(new_gen(2*i-1,60:62))==0) 
%     position=randi([1,2]);
%     new_gen(2*i-1,59+position)=abs(new_gen(2*i-1,59+position)-1);
% end
% if (new_gen(2*i-1,60:62)==forbid)
%     position=randi([1,2]);
%     new_gen(2*i-1,59+position)=abs(new_gen(2*i-1,59+position)-1);
% end
% if (sum(new_gen(2*i,60:62))==0) 
%     position=randi([1,2]);
%     new_gen(2*i,59+position)=abs(new_gen(2*i,59+position)-1);
% end
% if (new_gen(2*i,60:62)==forbid)
%     position=randi([1,2]);
%     new_gen(2*i,59+position)=abs(new_gen(2*i,59+position)-1);
% end
    
end




% end uniformx

%Work by Satadru Roy Binary Code developer Aug 11,2011
%Work by Satadru Roy GrayCode developer February 20,2010
function graycode=GCdevelop(bits)
GCrow=2^bits;
GCcol=bits+1;
graycode=zeros(GCrow,GCcol); %Defining the dimensions of the graycode table
for i=1:GCrow
    graycode(i,1)=i; % Serial Number for the graycode table
end
for col=GCcol:-1:2
    row=1;
    while (row<=GCrow)
        if (row==1)
            flip=0;
            for k=1:(2^(GCcol-col))
                graycode(row,col)=flip;
                row=row+1;
            end
            flip=1;
        end
        
        for k=1:(2^(GCcol+1-col))
            graycode(row,col)=flip;
            row=row+1;
            if (row>GCrow)
                break;
            end
        end
        flip=flip+1;
        flip=mod(flip,2);
    end
end

%Developed by Satadru on Nov 20,2011
%Generates the goal points based on the location of the objective function on the objective space
function [goals]=goalpointgen_new(obj_val,fmin)

pop_size=size(obj_val,1);
n=pop_size/8;

for ii=1:pop_size
   
    if ii<=2*n %direction vector has slope 0 %%%%Is pop size divided into 2n,2n&4n subpools in sequence?
        goals(ii,1)=fmin(1); %%%%minimize in obj. 1
        goals(ii,2)=obj_val(ii,2); %%%%goal taken same as obj. value for obj. 2 
    elseif ii>2*n && ii<=4*n %direction vector has slope 90
        goals(ii,1)=obj_val(ii,1); %%%%minimize in obj. 2
        goals(ii,2)=fmin(2); %%%%goal taken same as obj. value for obj. 1
    else %%%%for subpool 3
        xA_in=fmin(1); %%%%minimize obj. 1 
        yA_in=obj_val(ii,2);
        xB_in=obj_val(ii,1);
        yB_in=fmin(2); %%%%minimize obj. 2 

        f1=(obj_val(ii,1)-fmin(1))/fmin(1); %%%%percentage attainment of goal values for obj. 1
        f2=(obj_val(ii,2)-fmin(2))/fmin(2); %%%%percentage attainment of goal values for obj. 2
        x0=0;y0=0;
        xA=0;yA=f2;
        xB=f1;yB=0;
        
        m=tan(pi/2-atan(f2/f1));
        %x0,y0- Coordinates of the ideal point (or utopian point)

%         m=tan(pi/2-atan(((f2-y0)/(yA-y0))*((xB-x0)/(f1-x0)))); %Slope of the pointer

        x_int=((y0-f2)/m)+(f1-x0); 
        y_int=(m*(x0-f1))+(f2-y0);

        dV=sqrt(f1^2+(f2-y_int)^2);
        dH=sqrt((f1-x_int)^2+f2^2);
        % dV=sqrt((f1-x0)^2+(f2-y_int)^2);
        % dH=sqrt((f1-x_int)^2+(f2-y0)^2);

        [~,ind]=min([dV,dH]); %Pick the closest intercept

        gp_norm=[x0 y_int+y0;x_int+x0 y0]; %%%%normalized goal point?
        gp=(gp_norm.*[fmin(1) fmin(2);fmin(1) fmin(2)])+[fmin(1) fmin(2);fmin(1) fmin(2)];
        goals(ii,:)=gp(ind,:);     
    end
end

%Pareto generator developed by Satadru Roy
function [ND,xopt,NC_count]=paretogen(ND,obj_val,xopt,x_des)
ND_prev=ND;NC_count=0;
ND_new=[];xopt_new=[];
pot=[ND;obj_val];
xpot=[xopt;x_des];
f1=0;
while size(pot,1)>0
    [~,f2max_in]=max(pot(:,2));
    if (pot(f2max_in,1)==min(pot(:,1)))
        if f1==0
            ND_new=[ND_new;pot(f2max_in,:)];
            xopt_new=[xopt_new;xpot(f2max_in,:)];
            f1=1;
        else
            f2=0;
            for kk=1:size(ND_new,1)
                if ((ND_new(kk,1)==pot(f2max_in,1)) && (ND_new(kk,2)==pot(f2max_in,2)))
                    f2=1;
                    break
                end      
            end
            if f2==0
                ND_new=[ND_new;pot(f2max_in,:)];
                xopt_new=[xopt_new;xpot(f2max_in,:)];
            end
        end
    end            
    if f2max_in==1
        pot_new=pot(f2max_in+1:end,:);
        xpot_new=xpot(f2max_in+1:end,:);
    elseif f2max_in==size(pot,1)
        pot_new=pot(1:end-1,:);
        xpot_new=xpot(1:end-1,:);
    else
        pot_new=[pot(1:f2max_in-1,:);pot(f2max_in+1:end,:)];
        xpot_new=[xpot(1:f2max_in-1,:);xpot(f2max_in+1:end,:)];
    end
    pot=pot_new;
    xpot=xpot_new;
end

ND=ND_new;
xopt=xopt_new;

if (size(ND,1)==size(ND_prev,1))                                                                                                                                   
    if ND==ND_prev                                                                                                                                                 
         NC_count=1;                                                                                                                                                 
    end                                                                                                                                                        
end   