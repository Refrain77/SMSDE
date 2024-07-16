                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    %% ============ Two Primary Mutation Strategies and a Group of Secondary Ones Differential Evolution Algorithm (SMSDE) ============
% =========================================================================
clc
format short e;
%%  introductory Definitions
num_prbs = 10;                      %% number of test problems
max_runs=30;                        %% number of runs
outcome=zeros(max_runs,1);          %% to save the solutions of each run
com_time=zeros(max_runs,1);         %% Computational time
SR=zeros(max_runs,1);               %% How many times the optimal solution is obtained
Avg_FES=zeros(max_runs,1);          %% average fitness evaluations
Final_results=zeros(num_prbs,8);    %% to save the final results


%% run on more than one processor
%myCluster = parcluster('local');]

%myCluster.NumWorkers = 10;  % define how many processors to use

%% ========================= main loop ====================================
for dv=[15] % problem dimensions
    for I_fno= [1:10]
        Par= Introd_Par(I_fno,dv); %% set of parameters
        sol=zeros(30*1,Par.n); %% the best solution vector of each run
        vv=[];
        filename2 = strcat(strcat('results_\CEC20_F',num2str(I_fno)),'_SMSDE_',num2str(Par.n),'.txt');
       %parfor run=1:max_runs
       for run=1:max_runs
            nums = zeros(1, 4);
            [outcome(run),com_time(run),SR(run), Avg_FES(run),res, sol(run,:),nums]=SMSDE_main(run,I_fno,dv,nums);
%% 打印频率
%             fp2 = fopen(filename2,'a+');
%             for i = 1:length(nums)
%                 fprintf(fp2, 'op%d: %d ', i, nums(i)); % 写入位置和值
%             end
%             fprintf(fp2,'\n');
           %% to print the convergence of ech run % set 0 if not
            if Par.Printing==1
                res= res- repmat(Par.f_optimal,1,size(res,2));
                res(res<=1e-08)=0; ss=size(res,2);
                endv=res(ss);
                if size(res,2)<Par.Max_FES
                    res(size(res,2):Par.Max_FES)=endv;
                end
                vv(run,:)= res(1:Par.Max_FES);
            end

        end
        filename2 = strcat(strcat('results_\CEC20_results','_SMSDE_',num2str(Par.n),'.txt'));
        fp2 = fopen(filename2,'a+');
        filename = strcat(strcat('Fx_\CEC20_F',num2str(I_fno)),'_SMSDE_',num2str(Par.n),'.txt');
        fp = fopen(filename,'a+');
        fprintf(fp,'%.2e (%.2e): \r\n', mean(outcome), std(outcome)); 
        fprintf(fp2,'F%d: %.2e  (%.2e): ',I_fno , mean(outcome), std(outcome));
        for x = 1 : max_runs
            fprintf(fp,'%e ',outcome(x));   
            fprintf(fp2,'%e ',outcome(x));
        end
        fprintf(fp2,'\n');
        Final_results(I_fno,:)= [min(outcome),max(outcome),median(outcome), mean(outcome),std(outcome),mean(com_time),mean(SR),mean(Avg_FES)];

        disp(Final_results);


        %% fitness values at different levels of the optimization process
        %%% required by the competition
        if Par.Printing==1
            for k=1:16
                lim(k)=Par.n^(((k-1)/5)-3).*Par.Max_FES;
            end
            lim= ceil(lim);
    %         lim= [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].*Par.Max_FES;
            res_to_print= vv(:,lim);
            name1 = 'Results_Record\IMODE';
            name2 = num2str(I_fno);
            name3 = '_';
            name4 = num2str(Par.n);
            name5 = '.dat';
            f_name=strcat(name1,name2,name3,name4,name5);
            res_to_print=res_to_print';
            save(f_name, 'res_to_print', '-ascii');
            name1 = 'Results_Record\seeds_';
            f_name=strcat(name1,name2,name3,name4,name5);
            %% save the seeds used - if needed
    %         myMatrix2= double(seed_run);
    %         save(f_name, 'myMatrix2', '-ascii');

        end  
     end
end
