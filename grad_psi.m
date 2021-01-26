
function val0 = Grad_Psi(L,P)
    % P input TRANSITION matrix
    % L: vector of constant that appears in the loss function: i.e., Psi= sum(L(x,i)*pi_i,i)

    % this is done using chain rule and formula provided in "sensitivity analysis of discrete markov chains" by Hal Caswell

    % first we compute d pi/ d P using the formula
    % to do that we use the fundamental matrix of ergodic MC: Z
    mc=dtmc(P);
    pi_0=asymptotic
    d = length(P); %size of the matrix
    Z=(eye(d)-P+)
    d_pi = ()

        
end