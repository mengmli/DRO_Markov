function [x_true,alpha_real, val] = test_real(d,k,p_real,P,B,w)

    %calculating stationary distr.
    alpha_real=zeros(m,d);
    for i=1:d
    I = eye(m);
    
    Y = null(P_real(:,:,i)'-I);
    pi_approx = Y./(sum(Y));
    alpha_real(:,i)=pi_approx;
    end
    
    % opts = optimoptions('fmincon','Display','off');
    % [x,val]=fmincon(@(x)-(prime'.*x)'*alpha_real*w',ones(m,1),P,B,[],[],zeros(m,1),[]);
    % val0=-sumco*x+1/2*x'*(eye(d).*(1/c_exp(1)))*x; linprog(f,A,b,Aeq,beq,lb,ub)
    xrange=dec2bin(0:1:2^m-1)-'0';   
    
    val=10^6;
    for row=1:length(xrange(:,1))
        x=xrange(row,:)';
        if P*x<=B
            val1 = -(prime'.*x)'*alpha_real*w';
        else 
            val1 = 10^6;
        end
        if val1<val
            val=val1;
            x_true=x;
        end
    end
    % opts = optimoptions('intlinprog','Display','off');
    % [x,val]=intlinprog(@(x)-(prime'.*x)'*alpha_real*w',3,P,B,opts);
    % val0=-sumco*x+1/2*x'*(eye(d).*(1/c_exp(1)))*x; linprog(f,A,b,Aeq,beq,lb,ub)
    disp('true x');
    disp(x_true);
    end
    