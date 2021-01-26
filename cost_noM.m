function cost_val = cost_noM(d,x,prime,alpha_all,r,m)
cost_val=zeros(d,1);
cost=-prime'.*x;
for i=1:d
alpha0=alpha_all(i,:);
% disp(alpha0);
cbar=max(cost); % 
beta0 = 10^7;%(cbar-exp(-r)*dot(cost,alpha0))/(1-exp(-r));
% disp(cbar);
% disp(beta0);
su = @(eta_dual) 0;
for j=1:(m)
    su =@(eta_dual) su(eta_dual) + alpha0(j) * log(eta_dual-cost(j));
end
lambda =@(eta_dual) exp(su(eta_dual)-r); %compute lambda
    
f = @(eta_dual) eta_dual -lambda(eta_dual);
% opts = optimset('fminbnd','Display','off');
[eta_star,cost_val(i)] = fminbnd(f,cbar,beta0);
end
end

