%% ============ Two Primary Mutation Strategies and a Group of Secondary Ones Differential Evolution Algorithm (SMSDE) ============
% =========================================================================
%% Some part of this code is taken from UMOEA-II
%%==============================================================================
function [outcome,com_time,SR,avgFE,res_det,bestx,nums]= SMSDE_main(run,I_fno,dv,nums)

Par= Introd_Par(I_fno,dv);
iter=0;             %% current generation

filename = strcat(strcat('Fx_\CEC20_F',num2str(I_fno)),'_SMSDE_',num2str(Par.n),'.txt');
fp = fopen(filename,'a+');

%% =================== Define a random seed ===============================
%%-----Becase we ran experiments in parallel, we used "*run" to differentiate
%%-----among runs which started at the same time
RandStream.setGlobalStream (RandStream('mt19937ar','seed',run));

%% to record seeds for further validation, if needed
% seed_run=stream.Seed;

%% define variables
current_eval=0;             %% current fitness evaluations
PS1=Par.PopSize;            %% define PS1
PS2=100;
InitPop1=PS2;
%% ====================== Initalize x ==================================
x=repmat(Par.xmin,Par.PopSize,1)+repmat((Par.xmax-Par.xmin),Par.PopSize,1).*rand(Par.PopSize,Par.n);

%% calc. fit. and update FES
tic;
% fitx= cec22_func(x',I_fno); % 2022
fitx= cec20_func(x',I_fno); % 2020
current_eval =current_eval+Par.PopSize;
res_det= min(repmat(min(fitx),1,Par.PopSize), fitx); %% used to record the convergence

curdiv = 0.0;
for x1 =1 : Par.n
    midpoint(x1) = median(x(:,x1));
end
distobest = 1 : Par.PopSize;
for x1 = 1: Par.PopSize
    distobest (x1)= 0;
    for y = 1 : Par.n
        distobest(x1) = distobest(x1) + abs((x(x1,y) - midpoint(y))/(Par.xmax(y) - Par.xmin(y)));
    end
    distobest (x1) = distobest (x1) / Par.n;
    curdiv = curdiv + distobest (x1);
end
curdiv = curdiv / Par.PopSize;
fprintf(fp,'%d %e %e %e\r\n', current_eval, curdiv, mean(fitx-Par.f_optimal), min(fitx-Par.f_optimal));
disp(['nfes:' ,num2str(current_eval),'  curdiv:',num2str( curdiv), '  mean:',num2str(mean(fitx-Par.f_optimal)), '  min:',num2str(min(fitx-Par.f_optimal))] );

%% ====================== store the best ==================
[bestold, bes_l]=min(fitx);     bestx= x(bes_l,:);
%% ================== fill in for each  Algorithm ===================================


%% IMODE
EA_1= x(1:PS1,:);    EA_obj1= fitx(1:PS1);   EA_1old = x(randperm(PS1),:);
% Initialise op_3 score
op_3_list = [1,2,3,4];
init_score = ones(1,numel(op_3_list));
op_score = [op_3_list; init_score];
% initial position
op3_pos = 1;

%% CMA-ES
EA_2= x(PS1+1:size(x,1),:);    EA_obj2= fitx(PS1+1:size(x,1));
%% ================ define CMA-ES parameters ==============================
setting=[];bnd =[]; fitness = [];
[setting]= init_cma_par(setting,EA_2, Par.n, PS2);

%% ===== prob. of each DE operator
probDE1=1./Par.n_opr .* ones(1,Par.n_opr);
%% ===================== archive data ====================================
arch_rate=2.6;
archive.NP = arch_rate * PS1; % the maximum size of the archive
archive.pop = zeros(0, Par.n); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions
%% ==================== to adapt CR and F =================================
hist_pos=1;
memory_size=20*Par.n;
archive_f= ones(1,memory_size).*0.2;
archive_Cr= ones(1,memory_size).*0.2;
archive_T = ones(1,memory_size).*0.1;
archive_freq = ones(1, memory_size).*0.5;
%%
stop_con=0; avgFE=Par.Max_FES; InitPop=PS1; thrshold=1e-08;

cy=0;indx = 0; Probs=ones(1,2);
F = normrnd(0.5,0.15,1,PS1);
cr= normrnd(0.5,0.15,1,PS1);
%% main loop
xxx=0;

%% worst and oldest
age = zeros(PS1, 1);
prob_wls = 0.1;
%%
while stop_con==0
    
    iter=iter+1;
    cy=cy+1; % to control CS
    Probs = [1 0]; %% comment this line if the algorithm uses IMODE and CMA-ES
    %% ======================Applying SMSDE ============================
    if (current_eval<Par.Max_FES)
        if rand<=Probs(1)
            
            %% =============================== Linear Reduction of PS1 ===================================================
            UpdPopSize = round((((Par.MinPopSize - InitPop) / Par.Max_FES) * current_eval) + InitPop);
            if PS1 > UpdPopSize
                reduction_ind_num = PS1 - UpdPopSize;
                if PS1 - reduction_ind_num <  Par.MinPopSize
                    reduction_ind_num = PS1 - Par.MinPopSize;
                end
                %% remove the worst ind.
                for r = 1 : reduction_ind_num
                    vv=PS1;
                    EA_1(vv,:)=[];
                    EA_1old(vv,:)=[];
                    EA_obj1(vv)=[];
                    age(vv)=[];
                    PS1 = PS1 - 1; 
                end
                archive.NP = round(arch_rate * PS1);
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
            
            %% apply SMSDE
            [EA_1, EA_1old, EA_obj1,probDE1,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,archive_freq, current_eval,res_det,F,cr,op3_pos,op_score,age,nums] = ...
                SMSDE( EA_1,EA_1old, EA_obj1,probDE1,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,....
                archive_freq, Par.xmin, Par.xmax,  Par.n,  PS1,  current_eval, I_fno,res_det,Par.Printing,Par.Max_FES, Par.Gmax, iter,F,cr,op3_pos,op_score, age, nums);
            Ieval=current_eval;

            if (current_eval / Par.Max_FES) >= (xxx / 10)
                if  xxx > 0
                    curdiv = 0.0;
                    for x1 =1 : Par.n
                        midpoint(x1) = median(x(:,x1));
                    end
                    distobest = 1 : Par.PopSize;
                    for x1 = 1: Par.PopSize
                        distobest (x1)= 0;
                        for y = 1 : Par.n
                            distobest(x1) = distobest(x1) + abs((x(x1,y) - midpoint(y))/(Par.xmax(y) - Par.xmin(y)));
                        end
                        distobest (x1) = distobest (x1) / Par.n;
                        curdiv = curdiv + distobest (x1);
                    end
                    curdiv = curdiv / Par.PopSize;
                    fprintf(fp,'%d %e %e %e\r\n', current_eval, curdiv, mean(EA_obj1-Par.f_optimal), min(EA_obj1-Par.f_optimal));
                    %disp(['nfes:' ,num2str(current_eval),'  curdiv:',num2str( curdiv), '  mean:',num2str(mean(EA_obj1-Par.f_optimal)), '  min:',num2str(min(EA_obj1-Par.f_optimal))] );
                end
                xxx = xxx + 1;
            end
        end
    end
    %% ====================== CMA-ES ======================
    if (current_eval<Par.Max_FES)
        if   rand<Probs(2)
            disp("Use CMA-ES");
            UpdPopSize = round((((Par.MinPopSize1 - InitPop1) / Par.Max_FES) * current_eval) + InitPop1);
            if PS2 > UpdPopSize
                reduction_ind_num = PS2 - UpdPopSize;
                if PS2 - reduction_ind_num <  Par.MinPopSize
                    reduction_ind_num = PS2 - Par.MinPopSize;
                end
                %% remove the worst ind.
                for r = 1 : reduction_ind_num
                    vv=PS2;
                    EA_2(vv,:)=[];
                    %                     EA_1old(vv,:)=[];
                    EA_obj2(vv)=[];
                    PS2 = PS2 - 1;
                end
                %                 archive.NP = round(arch_rate * PS1);
                %                 if size(archive.pop, 1) > archive.NP
                %                     rndpos = randperm(size(archive.pop, 1));
                %                     rndpos = rndpos(1 : archive.NP);
                %                     archive.pop = archive.pop(rndpos, :);
                %                 end
            end
            [ EA_2, EA_obj2, setting,bestold,bestx,bnd,fitness,current_eval,res_det] = ...
                Scout( EA_2, EA_obj2, probSC, setting, iter,bestold,bestx,fitness,bnd,...
                Par.xmin,Par.xmax,Par.n,PS2,current_eval,I_fno,res_det,Par.Printing,Par.Max_FES);
            %             disp([min(EA_obj2) Probs(1) Probs(2)]);
        end
    end
    % ============================ LS2 ====================================
    if current_eval>0.85*Par.Max_FES && current_eval<Par.Max_FES
        if rand<Par.prob_ls
            old_fit_eva=current_eval;
            [bestx,bestold,current_eval,succ] = LS2 (bestx,bestold,Par,current_eval,I_fno,Par.Max_FES,Par.xmin,Par.xmax);
            if succ==1 %% if LS2 was successful
                EA_1(PS1,:)=bestx';
                EA_obj1(PS1)=bestold;
                age(PS1) = 0;
                [EA_obj1, sort_indx]=sort(EA_obj1);
                EA_1= EA_1(sort_indx,:);
                age = age(sort_indx);
                EA_2=repmat(EA_1(1,:), PS2, 1);
                [setting]= init_cma_par(setting,EA_2, Par.n, PS2);
                setting.sigma=1e-05;
                EA_obj2(1:PS2)= EA_obj1(1);
                Par.prob_ls=0.1;
            else
                Par.prob_ls=0.01; %% set p_LS to a small value it  LS was not successful
            end
            %% record best fitness -- set Par.Printing==0 if not
            if Par.Printing==1
                res_det= [res_det repmat(bestold,1,(current_eval-old_fit_eva))];
            end
        end
    end
%%
% 局部搜索，将又差又老的解局部搜索
 if current_eval<Par.Max_FES
    if rand < prob_wls
        PS11 = floor(PS1 * 0.7); % 遍历30%的x
        maxage = 0;
        id = PS1;
        for i = PS1:-1:PS11
            aggee = age(i);
            if(maxage < aggee)
                id = i;
                maxage = aggee;
            end
        end
        wo_x = EA_1(id, :);
        wo_xfit = EA_obj1(id);
%         局部搜索
        [worstx, worstold, current_eval, succ] = LS2(wo_x, wo_xfit, Par, current_eval, I_fno, Par.Max_FES, Par.xmin, Par.xmax); 
        if succ == 1 % 如果局部搜索成功
%             更新种群x和适应度fitx
            EA_1(id,:)=worstx';
            EA_obj1(id)=worstold;
            age(id) = 0;
            [EA_obj1, sort_indx]=sort(EA_obj1);
            EA_1= EA_1(sort_indx,:);
            age = age(sort_indx);
            prob_wls = 0.1;
        else 
            prob_wls = 0.01;
        end
    end
end

%%    
    if (current_eval>=Par.Max_FES && xxx==10)
        for x1 =1 : Par.n
            midpoint(x1) = median(x(:,x1));
        end
        distobest = 1 : Par.PopSize;
        for x1 = 1: Par.PopSize
            distobest (x1)= 0;
            for y = 1 : Par.n
                distobest(x1) = distobest(x1) + abs((x(x1,y) - midpoint(y))/(Par.xmax(y) - Par.xmin(y)));
            end
            distobest (x1) = distobest (x1) / Par.n;
            curdiv = curdiv + distobest (x1);
        end
        curdiv = curdiv / Par.PopSize;
        fprintf(fp,'%d %e %e %e\r\n', Ieval, curdiv, mean(EA_obj1-Par.f_optimal), min(EA_obj1-Par.f_optimal));
        %disp(['nfes:' ,num2str(current_eval),'  curdiv:',num2str( curdiv), '  mean:',num2str(mean(EA_obj1-Par.f_optimal)), '  min:',num2str(min(EA_obj1-Par.f_optimal))] );
        xxx=11;
    end
    
                
    
    %% ====================== stopping criterion check ====================
    if (current_eval>=Par.Max_FES-4*UpdPopSize)
%         stop_con=1;
        avgFE=current_eval;
    end
    if (current_eval>=Par.Max_FES && xxx==11)
        disp(current_eval);
        stop_con=1;
        avgFE=current_eval;
    end
    if ( (abs (Par.f_optimal - bestold)<= thrshold))
%         stop_con=1;
        bestold=Par.f_optimal;
        avgFE=current_eval;
    end
    
    %% =============================== Print ==============================
    %          fprintf('current_eval\t %d fitness\t %d \n', current_eval, abs(Par.f_optimal-bestold));
    if stop_con
        com_time= toc;%cputime-start_time;
        fprintf('run\t %d, fitness\t %d, avg.FFE\t %d\t %d\n', run, abs(Par.f_optimal-bestold),avgFE,indx(1));
        outcome= abs(Par.f_optimal-bestold);
        if (min (bestx))< -100 || (max(bestx))>100 %% make sure  that the best solution is feasible
            fprintf('in problem: %d, there is  a violation',I_fno);
        end
        SR= (outcome==0);
    end
end
end
