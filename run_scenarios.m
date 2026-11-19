% RUN_SCENARIOS  Reproduces Fig. 1 and Table II of
%   B. Zambotti, V. Andrieu, D. Astolfi, L. Bako, L. Zaccarian,
%   "Decentralized Stubborn Observer Redesign with Exact Outlier Cancellation".
%
% Compares the Centralized Stubborn Observer (CSO) and the Decentralized
% Stubborn Observer (DSO) on a platoon of N = 15 followers, with outliers
% injected on the distance sensor of Vehicle 9 (output y_17, Block 4).
%   Scenario 1: 1 <= T_dwell <= n_4      (Fig. 1, left column)
%   Scenario 2: T_dwell = n_4 + 1         (Fig. 1, right column)

clearvars; clc; close all;
addpath('src');
rng(1);                 % fixed seed: reproducible disturbance sequences

%% 1. System Setup
p = struct();
p.dt = 0.05;            % Sampling time [s]
p.N = 15;               % Followers
p.n_veh = 3;            % States per vehicle
p.n = p.N * p.n_veh;

% Physics parameters
p.tau = 0.3; p.h = 0.7; p.r = 5; p.kp = 0.2; p.kd = 3.5; p.ka = 0.1;

% Colors for plots
colors = [
    1 0.1 0.1;              % red    (DSO / attacked block)
    0.3, 0.3, 0.75;         % blue   (CSO)
    0.4940, 0.1840, 0.5560; % purple
    0.9290, 0.7940, 0.250   % yellow (DSO healthy blocks)
    ];

% State Space Model
[A, B, C] = model_platoon(p);
sys.A = A; sys.B = B; sys.C = C;
p.A = A;
p.B = B;
p.C = C;

% Observer Design (Decentralized LQR)
L_sparse = transpose(dlqr(A', C', 500,1));


% Numerical cleanup
tol = 1e-3;
L_sparse(abs(L_sparse)<tol) = 0;
fprintf('Spectral radius of A - L*C: %.4f\n', max(abs(eig(A-L_sparse*C))));


% Structural Decomposition
[Pd, n_part, m_part, At, Bt, Ct, Lt, T_global] = findAndTransformObserver(A, B, C, L_sparse);


% Build Distributed Stubborn Observer (DSO)
num_blocks = length(n_part);
r_params = zeros(1, num_blocks);

% Empirical tuning range for r_i (Section V of the paper)
r_min = 0.1;
r_max = 0.3;

n_max = max(n_part);
n_min = min(n_part);

for i = 1:num_blocks
    n_i = n_part(i);
    
    if n_max == n_min
        % All blocks have the same size
        r_params(i) = r_max; 
    else
        % Linear interpolation: r_i scales with n_i in [r_min, r_max]
        r_params(i) = r_min + (n_i - n_min) * (r_max - r_min) / (n_max - n_min);
    end
    
    fprintf('Block %d (size n=%d): r tuned to %.3f\n', i, n_i, r_params(i));
end

[S, G, R] = create_stubborn_matrices_robust(n_part, m_part, r_params);

obs_dist.A = At; obs_dist.B = Bt; obs_dist.C = Ct; 
obs_dist.L = Lt; obs_dist.S = S;  obs_dist.G = G; obs_dist.R = R;
obs_dist.T_global = T_global;

% Build Centralized Stubborn Observer (CSO)
tau_cent = p.n;               
S_cent = diag(ones(tau_cent-1, 1), -1); 
G_cent = zeros(tau_cent, 1); G_cent(1) = 1; 
Theta_cent = 0.2;
R_cent = Theta_cent * ones(1, tau_cent);

obs_cent.A = A; obs_cent.B = B; obs_cent.C = C; 
obs_cent.L = L_sparse; obs_cent.S = S_cent; 
obs_cent.G = G_cent; obs_cent.R = R_cent;
obs_cent.D = zeros(size(C,1), size(B,2));
obs_dist.D = obs_cent.D;
obs_cent.Theta = Theta_cent;
obs_cent.tau = tau_cent;

%% Identify the local block of the attacked vehicle
% Vehicle i sensor indices: 2i-1 (distance), 2i (velocity)
target_veh = 9;
target_sensor_idx = (target_veh-1)*2 +1; 
% Find block containing this sensor
local_blk_idx = 0;
for i = 1:length(Pd)
    if ismember(target_sensor_idx, Pd{i})
        local_blk_idx = i;
        break;
    end
end
if local_blk_idx == 0
    error('Could not find the block of Vehicle %d distance sensor.', target_veh);
end
n_local = n_part(local_blk_idx);
n_global = tau_cent;

fprintf('System Setup Complete.\n');
fprintf('Vehicle %d local memory (n_local): %d\n', target_veh, n_local);
fprintf('Global Memory (n_global): %d\n', n_global);

%% SCENARIO 1: Spatial Isolation 
fprintf('\n--- Running Scenario 1 ---\n');
% Goal: a fault on the attacked vehicle does not propagate to other blocks in the DSO.
% Disturbance: sporadic outliers, |d| in [5, 15] m, dwell time <= n_local.

T_sim_1 = 5.1;
num_steps_1 = floor(T_sim_1 / p.dt);
params_1.C = C;
params_1.target_sensor_idx = target_sensor_idx;
params_1.magnitude = 10;   % |d| = magnitude * U[0.5, 1.5] -> |d| in [5, 15] m
params_1.min_dwell = 2;

if n_local == 0
    params_1.max_dwell = 2;
else
    params_1.max_dwell = max(1, n_local-2); 
end
dist_hist_1 = generate_scenario_disturbance(params_1, num_steps_1, 'spatial_isolation');

% Simulation
[X_1, X_hat_d_1, X_hat_c_1, Delta1, Delta_centr1] = run_simulation(p, A, B, C, obs_dist, obs_cent, num_steps_1, dist_hist_1, 'cruise');

err_norm_d_1 = vecnorm(X_hat_d_1 - X_1);
err_norm_c_1 = vecnorm(X_hat_c_1 - X_1);

% Mean error norm
mean_d1 = mean(err_norm_d_1);
mean_c1 = mean(err_norm_c_1);

% RMSE 
rmse_d1 = sqrt(mean(err_norm_d_1.^2));
rmse_c1 = sqrt(mean(err_norm_c_1.^2));

fprintf('RMSE (DSO): %.4f\n', rmse_d1);
fprintf('RMSE (CSO): %.4f\n', rmse_c1);
%% SCENARIO 2: The Frequency Test
fprintf('\n--- Running Scenario 2 ---\n');
% Goal: DSO recovers, CSO accumulates.
% Disturbance: Period T_attack: n_local < T_attack < n_global.

T_attack = n_local + 1; 
fprintf('Attack Period: %d steps (n_local=%d < T < n_global=%d)\n', T_attack, n_local, n_global);
T_sim_2 = 5.1;
num_steps_2 = floor(T_sim_2 / p.dt);
params_2.C = C;
params_2.target_sensor_idx = target_sensor_idx;
params_2.magnitude = 10;   % |d| in [5, 15] m
params_2.T_attack = T_attack;
dist_hist_2 = generate_scenario_disturbance(params_2, num_steps_2, 'temporal_relaxation');

% Simulation
[X_2, X_hat_d_2, X_hat_c_2, Delta_2, Delta_centr_2] = run_simulation(p, A, B, C, obs_dist, obs_cent, num_steps_2, dist_hist_2, 'cruise');

% Metrics: global estimation error norm
err_norm_d_2 =(vecnorm(X_hat_d_2 - X_2));
err_norm_c_2 = (vecnorm(X_hat_c_2 - X_2));

mean_d2 = mean(err_norm_d_2);
mean_c2 = mean(err_norm_c_2);

% RMSE 
rmse_d2 = sqrt(mean(err_norm_d_2.^2));
rmse_c2 = sqrt(mean(err_norm_c_2.^2));


fprintf('RMSE (DSO): %.4f\n', rmse_d2);
fprintf('RMSE (CSO): %.4f\n', rmse_c2);

fprintf('Mean Error Norm (DSO): %.4f\n', mean(err_norm_d_2));
fprintf('Mean Error Norm (CSO): %.4f\n', mean  (err_norm_c_2));
save('scenario2_results.mat', 'X_2', 'X_hat_d_2', 'X_hat_c_2', 'err_norm_d_2', 'err_norm_c_2', 'dist_hist_2');

% --- PLOT ---
time_1 = (0:num_steps_1-1);
time_2 = (0:num_steps_2-1);

figure(1); 
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
eps_val = 1e-10;

% Plot 1: error Scenario 1
nexttile;
hold on; grid on; box on;
plot(time_1, max(err_norm_c_1, eps_val), 'Marker', 'o','MarkerSize',4,'Color', colors(2,:), 'LineWidth', 2); 
plot(time_1, max(err_norm_d_1, eps_val), 'Marker', 'o','MarkerSize',4,'Color', colors(1,:), 'LineWidth', 2);
ylabel('$|e|$', 'Interpreter', 'latex', 'FontSize', 16);
legend('CSO', 'DSO', 'Location', 'northwest', 'Interpreter', 'latex',  'FontSize', 10);
set(gca, 'FontSize', 14);
xlim([0 25])
ylim([0 15])

% Plot 2: error Scenario 2
nexttile;
hold on; grid on;box on;
plot(time_2, err_norm_c_2, 'Color',colors(2,:),'Marker', 'o','MarkerSize',4, 'LineWidth', 2); 
plot(time_2, err_norm_d_2,'Color', colors(1,:),'Marker', 'o','MarkerSize',4, 'LineWidth', 2);
legend('CSO', 'DSO', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 14);
xlim([0 25])
ylim([0 15])

% Plot 3: delta Scenario 1
nexttile;
hold on; grid on; box on; 
h_attacked = [];
h_healthy = [];

attacked_block_indices = Pd{local_blk_idx};
for i = 1:size(R,1)
    data_plot = max(Delta1(i,:), eps_val); 
      if ismember(i, attacked_block_indices)
            h = plot(time_1, data_plot,'Marker', 'o','MarkerSize',4, 'LineWidth', 2, 'Color', colors(1,:));
        if isempty(h_attacked), h_attacked = h; end 
    else
        h = plot(time_1, data_plot, 'Marker', 'o','MarkerSize',4,'LineWidth', 1.5, 'Color', colors(4,:));
        if isempty(h_healthy), h_healthy = h; end
    end
end
h_cen = plot(time_1, max(Delta_centr1, eps_val), '-','Marker', 'o','MarkerSize',4, 'Color', colors(2,:), 'LineWidth', 2); 
xlabel('$k$','Interpreter', 'latex', 'FontSize', 16); 
ylabel('$\Delta$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 14);
xlim([0 25])
ylim([0 10])
legend([h_cen, h_attacked, h_healthy], ...
       'CSO', ...
       'DSO Attacked', ...
       'DSO Healthy ', ...
       'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

% Plot 4: Delta Scenario 2
nexttile;
hold on; grid on; box on; 

h_attacked = [];
h_healthy = [];

attacked_block_indices = Pd{local_blk_idx};
for i = 1:size(R,1)
   if ismember(i, attacked_block_indices)
        % Attacked block
        h = plot(time_2, max(Delta_2(i,:), eps_val), 'LineWidth', 2, 'Marker', 'o','MarkerSize',4,'Color', colors(1,:));
        if isempty(h_attacked), h_attacked = h; end 
    else
        % Healthy blocks
        h = plot(time_2, max(Delta_2(i,:), eps_val), 'LineWidth', 1.5,'Marker', 'o','MarkerSize',4, 'Color', colors(4,:));
        if isempty(h_healthy), h_healthy = h; end
    end
end
h_cen = plot(time_2, max(Delta_centr_2, eps_val), '-', 'Marker', 'o','MarkerSize',4,'Color', colors(2,:), 'LineWidth', 2); 
xlabel('$k$','Interpreter', 'latex', 'FontSize', 16); 
set(gca, 'FontSize', 14);
xlim([0 25])
ylim([0 10])
legend([h_cen, h_attacked, h_healthy], ...
        'CSO', ...
       'DSO Attacked', ...
       'DSO Healthy ', ...
       'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
