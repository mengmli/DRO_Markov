function cost=FW_main(k,x_cur,epsilon,r,iter,q_all,alpha_all)
    %returns the cost vector k*1, each row represents 
    %the unnormalized cost value for the respective customer group

    %parameters
    d=length(q_all); % number of available brands

    % initialize cost
    cost=zeros(k,1); %negative profit

    for i=1:k
        %specify the dynamics of ith markov chain
        q=q_all(:,:,i); %transition kernel
        alpha0=alpha_all(:,i); % stationary distribution

        nu=sample(alpha0,r,q,d);
        nu_best=zeros(d); % best worst-case transition kernel

        for t=1:iter % iterate within the specified iteration number
            c = grad_psi(x_cur,nu); %define cost vector: notice that c is the negative gradient of the psi function  
            stt=linear_sub(alpha0,r,c,d,q); %subproblem, determining argmax for descent direction
        if isnan(stt)
            break
        end
        dir=stt-nu;
        g(t)=trace(dir*c);
        
        if g(t)<epsilon
            nu_best=nu; %disp(nu_best);
            break
        end
        
        buf_lin_search=-10^9;
        gammat=-1;
        for gammax = 0:0.001:1
            nu_mat=nu+gammax*dir;
            if sum(isinf(nu_mat(:)))>=1 || sum(isnan(nu_mat(:)))>=1
            udx=udx+1; break
            end
            nu_mat=nu_mat./sum(nu_mat,2);
            
            f_lin_search = Psi(doterm,m, nu_mat);
            if f_lin_search > buf_lin_search
                gammat = gammax;
                buf_lin_search = f_lin_search;
            end
        end
        if gammat==-2
            break
        end
        nu=nu+gammat*dir;
        if sum(isinf(nu(:)))>=1 || sum(isnan(nu(:)))>=1
                break
        end
    end
    if nu_best==zeros(d)
    nu_best=nu;
    end

    nu_best_mat=nu_best./sum(nu_best,2);


    mc=dtmc(nu_best_mat);
    pi_approx=asymptotics(mc);
    approx_cost=dot(x_cur,pi_approx);
    cost(k)=approx_cost;
    % disp(approx_cost);
end