
function val0 = Grad_Psi(prime,d,p)
    %p input matrix
    A = transpose(p);
    A(d,:) = ones(1,d);
    for i=1:(d-1)
        sumj=0;
        for j=1:(d)
            if i~=j
                sumj=sumj+p(i,j);
            end
        end
        A(i,i)= sumj;
    end
    val0=zeros(d);
    
    for i=1:d
        for j=1:d
            for k=1:d %compute partial derivative of  pi_k wrt pm^i_j
                Mdk=minor(A,d,k);
                Mso=so_minor(A,d,i,min(k,i),max(k,i)); %disp(k);disp(i);disp(j);
                if i==j
                    val0(i,i)=val0(i,i);
                elseif i==d
                    base=(-1)^(d+k+1)*det(Mdk)*((-1)^(i+j)*det(minor(A,j,i)))/(det(A)^2);
                    if i==k
                        val0(i,j)=val0(i,j)+ prime(k)*base;
                    else
                        val0(i,j)=val0(i,j)+ prime(k)*(base+(-1)^(i+j+d+k-1)*det(so_minor(A,d,j,min(k,i),max(k,i)))/det(A));
                    end
                elseif j==d
                    base=(-1)^(d+k+1)*det(Mdk)*(-det(minor(A,i,i)))/(det(A)^2);
                    if i==k
                        val0(i,j)=val0(i,j)+ prime(k)*base;
                    elseif i<k
                        val0(i,j)=val0(i,j)+ prime(k)*(base+(-1)^(d+k+1)*det(Mso)/det(A));
                    else
                        val0(i,j)=val0(i,j)+ prime(k)*(base+(-1)^(d+k)*det(Mso)/det(A));
                    end
                else
                    base=(-1)^(d+k)*det(Mdk)*(det(minor(A,i,i))-(-1)^(i+j)*det(minor(A,j,i)))/(det(A)^2);
                    if i==k
                        val0(i,j)=val0(i,j)+ prime(k)*base;
                    elseif i<k
                       val0(i,j)=val0(i,j)+prime(k)*(base+((-1)^(d+k+1)*det(Mso)+(-1)^(d+k+i+j)*det(so_minor(A,d,j,min(k,i),max(k,i))))/det(A));
                    else
                       val0(i,j)=val0(i,j)+prime(k)*(base+((-1)^(d+k)*det(Mso)+(-1)^(d+k+i+j+1)*det(so_minor(A,d,j,min(k,i),max(k,i))))/det(A)); 
                    end 
                
                end
    %             disp(val0(i,j));
            end
        end
    end
        
    end