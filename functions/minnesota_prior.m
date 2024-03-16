function [b, omega] = minnesota_prior(y, x, n, k, T, p, lambda)

    % The function assumes that the x has the intercept

    % Input:
    % n: number of variables;
    % k: number of regressors per equation;
    % T: observations per variable;
    % p: number of lags;

    % Output:
    % b: beta of minnesota prior;
    % omega: omega of minnesota prior;

    k_no_intercept = k-1; %number of betas per equation (n*p) without intercept

    %sigma of AR(1)
    sigma2_ar1 = size(1, n);
    for i=1:n
        y_ar1 = y(:, i);
        x_ar1 = [ones(T-p, 1) x(:, i+1)]; %we don't take in consideration the first variable (intercept)
        theta = inv(x_ar1'*x_ar1)*x_ar1'*y_ar1;
        sigma2_ar1(1, i) = (1/(T-p-1))*((y_ar1 - x_ar1*theta)'*(y_ar1 - x_ar1*theta));
    end
    
    %modify the array to make the construction of minneosta prior easier
    sigma2_ar1 = repmat(sigma2_ar1, 1, p);
    
    b_minnesota = zeros(n, k_no_intercept);
    omega_minnesota = zeros(1, k_no_intercept);
        
    for i=1:n
        s = 1;
        for j=1:k_no_intercept
            if j==i
                b_minnesota(i, j) = 1;
                omega_minnesota(1, j) = (lambda^2)*(1/1)*(1/sigma2_ar1(j));
            else
                if rem(j,p)==0
                    s = s + 1;
                end
                omega_minnesota(1, j) = (lambda^2)*(1/s^2)*(1/sigma2_ar1(j));
            end
        end
    end
    
    % We add the prior on the intercept with high diffuse prior variance
    b_minnesota = [zeros(n,1) b_minnesota]; %Add the prior corresponding to the intercept
    omega_minnesota = [10^6 omega_minnesota]; %Add the variance corresponding to the intercept

    b = b_minnesota';
    omega_minnesota = omega_minnesota';
    omega = diag(omega_minnesota);

end