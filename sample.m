function val = sample(alpha0,r,q,d)

    N_samples = 10^5;
    for n=1:N_samples
        fprintf('sample order %d',n);
        % randomly sample p if r>=10^(-4)
        p=max(q+0.2*r.*rand(d)-0.1*r,10^(-5));
        p = p./sum(p,2);
        Dc = 0;
        for i=1:d
            sum0=0;
            for j=1:d
                if (p(i,j)~=0) 
                    if (q(i,j)~=0)
                     sum0 = sum0 + q(i,j)*(log(q(i,j))-log(p(i,j)));
                    end
                end
            end
            Dc = Dc + alpha0(i)*sum0;
        end
        fprintf('Dc %d\n',Dc);
        if Dc <= r
            fprintf('Dc %d',Dc);
            val = p; 
            break
        end
        if isnan(Dc)
    %         dbstop if naninf 
            disp('p');
            disp(p);
    %         return
        end
    end
    
        if isnan(Dc)
    %         dbstop if naninf 
            disp('p');
            disp(p);
            error('no right sample');
        end
    end
    