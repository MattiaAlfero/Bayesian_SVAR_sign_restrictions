function [satisfied, Gamma]  = sign_restrictions(S_post, restrictions)

    % Draw Gamma such that: Gamma*Gamma' = S_post

    % Input:
    % S_post: variance-covariance matrix;
    % restrictions: matrix of 1,0,-1 with the same dimension of S_post,
    %   used to identify the shocks.
    %   -1: negative impact, 1: positive impact, 0: not a specific effect;

    % Output:
    % Gamma: matrix drawn with QR decomposition such that Gamma*Gamma' =
    %   S_post;
    % satisfied: dummy variable specifying if the Gamma satisfy the sign
    %   restriction;

    n = size(S_post, 1);
    Lambda = chol(S_post, 'lower'); 
    H = randn(n);
    [Q, R] = qr(H);
    % Normalize R to be positive
    D = diag(sign(diag(R)));
    Q = Q*D; 
    R = D*R;
    % Construct Gamma
    Gamma = Lambda*Q;

    violations = 0;
    n = size(Gamma, 1);
    for i1=1:n
        for i2=1:n
            if restrictions(i1, i2)==-1
                violations = violations + (Gamma(i1,i2) > 0);
            elseif restrictions(i1, i2)==1
                violations = violations + (Gamma(i1,i2) < 0);
            else 
                violations = violations + 0;
            end
        end
    end
    
    if violations == 0
        satisfied = 1;
    else
        satisfied = 0;
    end

end