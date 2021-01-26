% for executing all the subfunctions and producing final plots

%  main workflow for the code:
%     I. Sampling & Preprocessing
%         1. read the generated samples from predefined markov chain dynamics
%         2. estimate the transition kernels from the generated samples
%     II. FW step
%         return the predictors given a prescriptor value
%     III. Plot
clear all
close all
rng(1);
dbstop if error
warning('off');

% parameter setting

k=10; % how many customer segements: i.e., how many different markov chain dynamics
d=12; % how many brands
T=30; % length of each xi^(i)/ sample size
n_exper = 10; % number of independent experiments

P =  .15*rand(k,d); %pricing information for each brand and price sensitivity of each respective customer group
B = 15*rand(k,1);
w=rand(1,k); % weight of each customer segment
w=w./sum(w);

xrange=dec2bin(0:1:2^d-1)-'0'; % decision space
x_feasible=[];
% keeping only the feasible decisions
for row=1:length(xrange(:,1)) % iterate over all possibilities
    x=xrange(row,:)'; %fix one decision
    if P*x<=B
        x_feasible=[x_feasible;x']; % stack feasible solutions on top of each other
    end
end
% read data from mat file

matObj = matfile('xi_new.mat','Writable',true);
xi=zeros(k,T,n_exper); % for each experiment, we have a distinct xi trajectory, each row of xi(:,:,i) represents the trajectory of a customer group
for n=1:n_exper
   xi(:,:,n)=(reshape(matObj.xi(n,:),[T,k])).';
end

epsilon=0.01; % error tolerance for FW gap
iter=10; % maximum iteration for FW alg

cost_fin=zeros(1,n_exper); 
cost_out=zeros(1,n_exper); 
% 
cost_fin_iid=zeros(1,n_exper);
cost_out_iid=zeros(1,n_exper);


r_range=logspace(-4,0,5); %range r
[~,N_r]=size(r_range);
markov_perf=zeros(1,N_r);
markov_perf_lower=zeros(1,N_r);
markov_perf_upper=zeros(1,N_r);

reliability=zeros(1,N_r);
reliability_iid=zeros(1,N_r);

iid_perf=zeros(1,N_r);
iid_perf_lower=zeros(1,N_r);
iid_perf_upper=zeros(1,N_r);

% assuming perfect information:
obj_tran_kernel = matfile('transition_dynamics.mat','Writable',true);
p_real=obj_tran_kernel.P;
[true_x,alpha_real,cost_real] = test_real(d,k,p_real,P,B,w);

for i=1:N_r % for each prescribed radius
    r=r_range(i);
    fprintf('order r %d ',i);
%     disp(r);
    n_disappt=0;
    for n=1:n_exper % run n_exper independent experiments
        fprintf('exp order %d ',n);
        prog=(n_exper*(i-1)+n)/(n_exper*length(r_range));
        fprintf('progress %0.2f\n', prog)
        % estimate stationary distribution
        [alpha0,q]=est_alpha_from_xi(k,d,T,xi(:,:,n));
        % iterate over strategy set
        cost_fin(n)=10^6; % default value for the predicted cost (negative profit)
        tic
        for row=1:length(x_feasible(:,1))
            x_cur=x_feasible(row,:)'; %fix one decision
            
            % apply FW algorithm to get the corresponding prediction
            cost_fin1 = w*FW_main(k,x_cur,epsilon,r,iter,q,alpha0); %return the best minimax prediction given a decision x_cur
            if cost_fin1<cost_fin(n) %compare if it is the best decision so far
                cost_fin(n)=cost_fin1;
                x=x_cur;
            end
        end
        elapsed_time=toc
        % return and store optimal decision and optimal value
        cost_out(n) = -x'*alpha_real*w'; % negative profit in the real situation
        if cost_out(n)>cost_fin(n)
            n_disappt=n_disappt+1;
        end
    end

    reliability(i)=1-n_disappt/n_exper; % empirical reliabitliy

    markov_perf(i)=mean(cost_out);
    markov_perf_lower(i)=markov_perf(i)-2*std(cost_out);
    markov_perf_upper(i)=markov_perf(i)+2*std(cost_out); % 90% confidence bound
end
