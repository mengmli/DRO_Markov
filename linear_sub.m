function s = linear_sub(alpha0,r,c,m,q)
    %outputs arg min <S, c> where c is the negative gradient
    cbar=max(c,[],2); 
    if r<0.1
    eta_0=cbar+10^(-8);
    else
        eta_0=cbar+10^(-1);
    end
    su = @(eta_dual) alpha0'*(sum(q.*(log(repmat( eta_dual, [1,m] )-c)-log(repmat( alpha0, [1,m] ))),2));

    lambda =@(eta_dual) exp(su(eta_dual)-r); %compute lambda
        
    f = @(eta_dual) sum(eta_dual) -lambda(eta_dual);
    y=f(eta_0);
    if isnan(y)
        disp('nan f');
        disp(cbar);
        s=nan;
        return
    end
    % disp('initial value');
    % disp(y);
    % disp('initial value for lambda');    
    % disp(lambda(eta_0));
    % disp('initial value for sum inside lambda');    
    % disp(su(eta_0));
    % disp('initial value for f');    
    % disp(f(eta_0));
    opts = optimoptions('fmincon','Display','off');
    options = optimoptions(opts,'MaxFunctionEvaluations', 5000);
    options = optimoptions(options,'MaxIterations', 5000);
    options = optimoptions(options,'OptimalityTolerance', 0.00001);
    options = optimoptions(options,'FunctionTolerance', 0.00001);
    options = optimoptions(options,'StepTolerance', 0.00001);
    eta_star = fmincon(f,eta_0,[],[],[],[],cbar,[],[],opts);
    
    % disp(eta_star-cbar);
    if (eta_star-cbar>=ones(size(eta_0)))
        s=nan;
        return 
    end
    sum2 = alpha0'*(sum(q.*(log(repmat( eta_star, [1,m] )-c)-log(repmat( alpha0, [1,m] ))),2));
    lambda_star = exp(sum2-r);
    % disp(lambda_star);disp(eta_star);
    
    s=zeros(m);
    
    
    for i=1:m
    %     disp('pi for q');
    %     disp(alpha0(i));
        for j=1:m
            if alpha0(i)~=0
            s(i,j)=lambda_star*alpha0(i)*q(i,j)/(eta_star(i)-c(i,j));
            if isnan(s(i,j))
                s(i,j)=2*10^(-1);
            % error('optimizer nan');
            end
            end
        end
    end
    
    end
    
    
    