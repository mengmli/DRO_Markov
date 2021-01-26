function val = sample(alpha0,r,q,d)
    % alpha0=alpha0(:,:,1);
    % d=4;
    % q=q(:,:,1);
    N_samples = 10^5;
    for n=1:N_samples
        fprintf('sample order %d',n);
        % randomly sample pif r>=10^(-4)
        p=max(q+0.1*r*rand(d),10^(-7));
        p = p./sum(p,2);
        Dc = 0;
        for i=1:d
            sum0=0;
            for j=1:d
                if (q(i,j)~=0)&&(p(i,j)~=0)
                     sum0 = sum0 + q(i,j)*log(q(i,j)/p(i,j));
                end
            end
            Dc = Dc + alpha0(i)*sum0;
        end
    
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
    