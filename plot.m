
for i=1:N_r % for each prescribed radius
    r=r_range(i);
    fprintf('order r %d ',i);
%     disp(r);
    n_disappt=0;
    for n=1:n_exper % run n_exper independent experiments
    alpha0=naive_est_alpha(d,k,T,xi(:,:,nsample));
    cost_fin_iid(nsample)=10^6;
    for row=1:length(xrange(:,1))
        xrow=xrange(row,:)';
        if P*xrow<=B
            cost_fin1 = w*cost_noM(d,xrow,prime,alpha0,r,m);
        else 
            cost_fin1 = 10^6;
        end
        if cost_fin1<cost_fin_iid(nsample)
            cost_fin_iid(nsample) =cost_fin1;
            x=xrow;
        end
    end
    cost_out_iid(nsample) = -(prime'.*x)'*alpha_real*w';
    if cost_out_iid(nsample)>cost_fin_iid(nsample)
        n_disappt=n_disappt+1;
    end
end
reliability_iid(i)=1-n_disappt/N_sample;

iid_perf(i)=mean(cost_out_iid); 
iid_perf_lower(i)=w0(i)-2*std(cost_out_iid);
iid_perf_upper(i)=w0(i)+2*std(cost_out_iid);
end


figure(1)
hold on;
hmeaniid=plot(xr, reliability_iid, 'LineWidth',2);
hmeaniid.Color='r';
hmean_m=plot(xr, reliability,'LineWidth',2);
hmean_m.Color='b';
set(gca,'XScale','log')

hold off;
%{
figure(2)
hold on;
%markov
x3 = [xr, fliplr(xr)];
inBetween = [z1, fliplr(z2)];
h2=fill(x3, inBetween, 'b','Edgecolor', 'none');
set(h2,'FaceAlpha',0.2)
hmeanout=plot(xr,z, 'b', 'LineWidth', 2);

%iid
x4 = [xr, fliplr(xr)];
inBetween2 = [w1, fliplr(w2)];
h3=fill(x4, inBetween2, 'r','Edgecolor', 'none');
set(h3,'FaceAlpha',0.2)
hmeanout_iid=plot(xr,w0, 'r', 'LineWidth', 2);

real=plot(xr,cost_real+0*xr, 'g', 'LineWidth', 2);
xlabel('r')
ylabel('cost')
% legend([hline],{'True cost'})
set(gca,'XScale','log')
hold off;