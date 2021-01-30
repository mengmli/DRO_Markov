function [x_true,alpha_real, val] = test_real(d,k,p_real,P,B,w,x_feasible)

    %calculating stationary distr.
    alpha_real=[];
    for i=1:k
        mc=dtmc(p_real(:,:,i));
        pi_approx = asymptotics(mc);
        alpha_real=[alpha_real pi_approx']; %stack stationary distribution to the right side of all stat. distributions
    end
    
    xrange=x_feasible;
    val=10^6;
    for row=1:length(xrange(:,1))
        x=xrange(row,:)';
        if P*x<=B
            val1 = -x'*alpha_real*w';
        else 
            val1 = 10^6;
        end
        if val1<val
            val=val1;
            x_true=x;
        end
    end
    
    disp('true x');
    disp(x_true);
    end
    