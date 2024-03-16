%% Bayesian SVAR identified with sign restrictions
% The main sections of the code are the following:
% 1) Posterior distribution of parameters under flat prior
% 2) Identification with sign restrictions under flat prior
% 3) Identification with sign restriction under Minnesota prior
% 4) Reconstruct the remaining part of the sample using only the selected shock

%% Housekeeping
clear;
clc;
rng("default")
addpath("functions")
                                                                                    
%% Import the data and set main parameters

%Import and arrange the data
raw_data = readtable("data.xlsx");
set.dates = table2array(raw_data(:, 1));
set.data = table2array(raw_data(:, 2:end));

% Manual parameters
set.irf_horizon = 20; %IRFs horizon
set.T = 124; %lenght used to estimate the sample
set.growth_rate = 4; %Since the data is quarterly, we are considering annual change
set.p = 4; %number of lags
set.select_shock = 3; 
set.restrictions = [1 1 1; 1 -1 0; 1 1 -1]; % -1: negative, 1: positive, 0: neutral, notice that the matrix has to be set.n x set.n
set.variable_names = {"C", "CoreP", "EnergyP"}; %it has to be the same as the order on x
set.shock_names = {"Demand shock", "Supply shock", "Energy shock"}; %it has to be the same as the order on x

% Default parameters
set.lambda = 0.2;
set.n = size(set.data, 2); %number of variables
set.k = set.n*set.p + 1; %number of regressors per variable + constant
set.mc_sim = 50000;

clear raw_data

%% Rearranging the data

figure('Position',[300 100 900 600]);
for i=1:set.n
    subplot(set.n, 1, i);
    plot(set.dates, set.data(:, i));
    title(set.variable_names(i))
end

set.x = [];
% Construct the X variable with p lags
for i=1:set.p
    j = set.p - i + 1;
    set.x = [set.data(i:end-j, :) set.x];
end
set.x = [ones(size(set.x,1), 1) set.x];
set.y = set.data(set.p+1:end,:);

% Now we select the obs up to set.T, since we don't use the entire sample
% for the estimation, and we construct the matrices used for the estimation
x = set.x(1:set.T-set.p, :);
y = set.y(1:set.T-set.p, :);
X = kron(eye(set.n), x);
Y = y(:);

clear i j

%% OLS estimates - VAR(p)

ols.B = inv(x'*x)*x'*y;
ols.b = inv(X'*X)*X'*Y;
ols.S = (y - x*ols.B)'*(y - x*ols.B);

%% 1): Posterior distribution of parameters under flat prior

flat.b_post = zeros(set.mc_sim, set.k*set.n, 1);
flat.S_post = zeros(set.mc_sim, set.n, set.n);

for i=1:set.mc_sim
    flat.S_post(i, :, :) = iwishrnd(ols.S, set.T-set.p-set.n-set.k-1);
    flat.b_post(i, :, :) = mvnrnd(ols.b, kron(squeeze(flat.S_post(i, :, :)), inv(x'*x)));
end

%Check 
flat.b_post_avg = mean(flat.b_post, 1)';
flat.S_post_avg = mean(flat.S_post); 
flat.S_post_avg = squeeze(flat.S_post_avg(1,:,:)).*(set.T-set.p-2*set.n-set.k-2);

flat.S_post_median = median(flat.S_post); 
flat.S_post_median = squeeze(flat.S_post_median(1,:,:));
disp(flat.S_post_median)

clear i 

%% 2): Identification with sign restrictions

flat.b_post = zeros(set.mc_sim, set.k*set.n, 1); % Posterior draws of beta
flat.S_post = zeros(set.mc_sim, set.n, set.n); % Posterior draws of sigma
flat.Gamma = zeros(set.mc_sim, set.n, set.n); %Matrix with identified shocks
flat.index_sign = []; %Index of matrices satisfying the sign restriction

for i=1:set.mc_sim
    flat.S_post(i, :, :) = iwishrnd(ols.S, set.T-set.p-set.n-set.k-1);
    flat.b_post(i, :, :) = mvnrnd(ols.b, kron(squeeze(flat.S_post(i, :, :)), inv(x'*x)));
    % Check if the matrix satisfy the sign restriction, in case, store it.
    [satisfied, flat.Gamma(i, :, :)] = sign_restrictions(squeeze(flat.S_post(i, :, :)), set.restrictions);
    if satisfied == 1
        flat.index_sign = [flat.index_sign i];
    end
end
flat.b_post_sign = flat.b_post(flat.index_sign, :, :);
flat.Gamma_sign = flat.Gamma(flat.index_sign, :, :);

clear satisfied i

% Put VAR(4) in companion form and compute IRFs
flat.IRF = zeros(size(flat.index_sign, 2), set.n*set.n, set.irf_horizon+1);
for i=1:size(flat.index_sign, 2)
    % Reshape the parameter at every iteration
    B_post_draw = reshape(flat.b_post_sign(i,:,:), [], set.n);
    Gamma_draw = squeeze(flat.Gamma_sign(i, :, :));
    % Put B in companion form, select the betas without the intercept
    B_post_draw = [B_post_draw(2:end, :)'; eye(set.n*(set.p-1)) zeros(set.n*(set.p-1), set.n)];
    Gamma_draw = [Gamma_draw zeros(set.n, set.n*(set.p-1)); zeros(set.n*(set.p-1), set.n) zeros(set.n*(set.p-1), set.n*(set.p-1))];
    % Let's generate the IRFs for the selected horizon
    for h=0:set.irf_horizon
        shocks = (B_post_draw^h)*Gamma_draw; %Use the companion form to generate the shocks
        shocks = shocks(1:set.n, 1:set.n);
        flat.IRF(i, :, h+1) = shocks(:); % 1st shock: 1:set.n, 2nd shock: set.n+1:2*set.n, ...
    end
end

%Plot IRFs
figure('Position',[300 100 900 600]);
j=1;
for i=1:(set.n*set.n)
    subplot(set.n, set.n, i);
    plot(0:set.irf_horizon, prctile(squeeze(flat.IRF(:, i, :)), 50)', ...
        0:set.irf_horizon, prctile(squeeze(flat.IRF(:, i, :)), 5)', "--r", ...
        0:set.irf_horizon, prctile(squeeze(flat.IRF(:, i, :)), 95)', "--r");
    yline(0);
    title(sprintf("%s to %s", set.shock_names{j}, set.variable_names{i - set.n*(j-1)}))
    if rem(i, set.n) == 0; j = j+1; end
end
leg = legend("Median", "90th credible intervals", 'Location','southoutside','orientation','horizontal');
leg.Position(1) = 0.4;
leg.Position(2) = 0.01;
sgtitle("IRFs with flat prior and sign restrictions")
saveas(gcf,'images/IRFs_flat_prior.jpg', 'jpg')

clear B_post_draw Gamma_draw i shocks h j


%% 3): Identification with sign restriction under Minnesota prior

minnesota.b_post = zeros(set.mc_sim, set.k*set.n, 1); % Posterior draws of beta
minnesota.S_post = zeros(set.mc_sim, set.n, set.n); % Posterior draws of sigma
minnesota.Gamma = zeros(set.mc_sim, set.n, set.n); %Matrix with identified shocks
minnesota.index_sign = []; %Index of matrices satisfying the sign restriction

%Compute the Minnesota prior
[minnesota.b_prior, minnesota.omega_prior] = minnesota_prior(y, x, set.n, set.k, set.T, set.p, set.lambda);

% Posterior locations of the Normal(location, scale) and IW(location, scale)
b_location = inv(x'*x + inv(minnesota.omega_prior))*(x'*y + inv(minnesota.omega_prior)*minnesota.b_prior);
S_location = (y - x*b_location)'*(y - x*b_location) + ...
    (b_location - minnesota.b_prior)'*inv(minnesota.omega_prior)*(b_location - minnesota.b_prior);

for i=1:set.mc_sim
    minnesota.S_post(i, :, :) = iwishrnd(S_location, set.T-set.p-set.n-1);
    minnesota.b_post(i, :, :) = mvnrnd(b_location(:), kron(squeeze(minnesota.S_post(i, :, :)), inv(x'*x + inv(minnesota.omega_prior))));
    % Check if the matrix satisfy the sign restriction, in case, store it.
    [satisfied, minnesota.Gamma(i, :, :)] = sign_restrictions(squeeze(minnesota.S_post(i, :, :)), set.restrictions);
    if satisfied == 1
        minnesota.index_sign = [minnesota.index_sign i];
    end
end
minnesota.b_post_sign = minnesota.b_post(minnesota.index_sign, :, :);
minnesota.Gamma_sign = minnesota.Gamma(minnesota.index_sign, :, :);

clear Lambda H Q R satisfied i

% Put VAR(4) in companion form and compute IRFs
minnesota.IRF = zeros(size(minnesota.index_sign, 2), set.n*set.n, set.irf_horizon+1);
for i=1:size(minnesota.index_sign, 2)
    % Reshape the parameter at every iteration
    B_post_draw = reshape(minnesota.b_post_sign(i,:,:), [], set.n);
    Gamma_draw = squeeze(minnesota.Gamma_sign(i, :, :));
    % Put B in companion form, select the betas without the intercept
    B_post_draw = [B_post_draw(2:end, :)'; eye(set.n*(set.p-1)) zeros(set.n*(set.p-1), set.n)];
    Gamma_draw = [Gamma_draw zeros(set.n, set.n*(set.p-1)); zeros(set.n*(set.p-1), set.n) zeros(set.n*(set.p-1), set.n*(set.p-1))];
    % Let's generate the IRFs for the selected horizon
    for h=0:set.irf_horizon
        shocks = (B_post_draw^h)*Gamma_draw; %Use the companion form to generate the shocks
        shocks = shocks(1:set.n, 1:set.n);
        minnesota.IRF(i, :, h+1) = shocks(:); % 1st shock: 1:set.n, 2nd shock: set.n+1:2*set.n, ...
    end
end

%Plot IRFs
figure('Position',[300 100 900 600]);
j=1;
for i=1:(set.n*set.n)
    subplot(set.n, set.n, i);
    plot(0:set.irf_horizon, prctile(squeeze(minnesota.IRF(:, i, :)), 50)', ...
        0:set.irf_horizon, prctile(squeeze(minnesota.IRF(:, i, :)), 5)', "--r", ...
        0:set.irf_horizon, prctile(squeeze(minnesota.IRF(:, i, :)), 95)', "--r");
    yline(0);
    title(sprintf("%s to %s", set.shock_names{j}, set.variable_names{i - set.n*(j-1)}))
    if rem(i, set.n) == 0; j = j+1; end
end
leg = legend("Median", "90th credible intervals", 'Location','southoutside','orientation','horizontal');
leg.Position(1) = 0.4;
leg.Position(2) = 0.01;
sgtitle("IRFs with Minnesota prior and sign restrictions")
saveas(gcf,'images/IRFs_minnesota_prior.jpg', 'jpg')

clear b_location S_location B_post_draw Gamma_draw i shocks h

%% 4): Reconstruct the remaining part of the sample using only the selected shock

% Select the correct sample
x = set.x(set.T+1-set.p-set.growth_rate:end, :);
y = set.y(set.T+1-set.p-set.growth_rate:end, :);
X = kron(eye(set.n), set.x);
Y = set.y(:);

%Pre-allocate variables
counterfactual.eps = zeros(size(minnesota.index_sign,2), size(y,1), set.n);
counterfactual.eps_fake = zeros(size(minnesota.index_sign,2), size(y,1), set.n);
counterfactual.y_fake = zeros(size(minnesota.index_sign,2), size(y,1), set.n);
counterfactual.y_fake_perc = zeros(size(minnesota.index_sign,2), size(y,1)-set.growth_rate, set.n);
counterfactual.y_expected = zeros(size(minnesota.index_sign,2), size(y,1), set.n);
counterfactual.y_expected_perc = zeros(size(minnesota.index_sign,2), size(y,1)-set.growth_rate, set.n);

for i=1:size(minnesota.index_sign,2)
    b = reshape(squeeze(minnesota.b_post_sign(i,:,:)), [], set.n);
    Gamma = squeeze(minnesota.Gamma_sign(i,:,:));
    counterfactual.eps(i, :, :) = (inv(Gamma)*(y - x*b)')';
    counterfactual.eps_fake(i, :, set.select_shock) = counterfactual.eps(i, :, set.select_shock); %Create the fake eps with 0 everywhere else
    counterfactual.y_fake(i, :, :) = x*b + (Gamma*squeeze(counterfactual.eps_fake(i, :, :))')';
    counterfactual.y_fake_perc(i, :, :) = squeeze(counterfactual.y_fake(i, set.growth_rate+1:end, :)) - squeeze(counterfactual.y_fake(i, 1:end-set.growth_rate, :));
    counterfactual.y_expected(i, :, :) = x*b;
    counterfactual.y_expected_perc(i, :, :) = squeeze(counterfactual.y_expected(i, set.growth_rate+1:end, :)) - squeeze(counterfactual.y_expected(i, 1:end-set.growth_rate, :));
end

%Compute the annual change
y_obs_perc = y(set.growth_rate+1:end,:) - y(1:end-set.growth_rate,:);
y_fake_perc = squeeze(median(counterfactual.y_fake_perc(:, :, :)));
y_expected_perc = squeeze(median(counterfactual.y_expected_perc(:, :, :)));

figure('Position',[300 100 900 600]);
for i=1:set.n
    subplot(set.n, 1, i);
    plot(set.dates(set.T+1:end, :), y_obs_perc(:, i), "-o", ...
        set.dates(set.T+1:end, :), y_expected_perc(:, i), "--", ...
        set.dates(set.T+1:end, :), y_fake_perc(:, i), "--o")
    title(sprintf("Counterfactual path of %s with %s", set.variable_names{i}, set.shock_names{set.select_shock}))
end
leg = legend("Observed", "Expected", "Counterfactual with shock", 'Location','southoutside','orientation','horizontal');
leg.Position(1) = 0.3;
leg.Position(2) = 0.01;
saveas(gcf, sprintf("images/Counterfactual_path_%s", set.shock_names{set.select_shock}), 'jpg')

figure;
plot(set.dates(set.T+1-set.growth_rate:end, :), squeeze(median(counterfactual.eps)))
legend(set.shock_names)
saveas(gcf,'images/shocks.jpg', 'jpg')