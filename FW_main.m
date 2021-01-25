function cost=FW_main(x_cur,epsilon,r,iter,q_all,alpha_all)
    %returns the cost vector k*1, each row represents 
    %the unnormalized cost value for the respective customer group

    %parameters
    k=length(x_cur); % number of customer groups
    d=length(q); % number of available brands

    % initialize cost
    cost=zeros(k,1); %negative profit

    for i=1:k
        %specify the dynamics of ith markov chain
        q=q_all(:,:,i); %transition kernel
        alpha0=alpha_all(:,i); % stationary distribution

        nu=sample(alpha0,r,q,m);
        nu_best=zeros(m); % best worst-case transition kernel

        each


end