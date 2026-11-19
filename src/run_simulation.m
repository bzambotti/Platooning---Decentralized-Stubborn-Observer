% RUN_SIMULATION  Simulates the platoon together with the Decentralized (DSO)
% and Centralized (CSO) finite-memory stubborn observers.
%   maneuver: 'cruise' (leader accelerates to 20 m/s, then cruises) or
%             'braking' (hard brake for 5 s <= t <= 10 s).
function [X, X_hat_phys, X_hat_centr, Delta, Delta_centr] = run_simulation(p, A, B, C, obs_dist, obs_cent, num_steps, dist_hist, maneuver)

    X = zeros(p.n, num_steps);
    X_hat_z = zeros(p.n, num_steps);
    Sigma = zeros(size(obs_dist.S,1), num_steps);
    Delta = zeros(size(obs_dist.R,1), num_steps); % DSO thresholds
    X_hat_centr = zeros(p.n, num_steps);
    Sigma_centr = zeros(size(obs_cent.S,1), num_steps);
    Delta_centr = zeros(1, num_steps); 
    
    % Initial Conditions
    v_ref = 0;
    for i = 1:p.N
        idx_v = (i-1)*p.n_veh + 2;
        X(idx_v, 1) = v_ref;
        X_hat_z(idx_v, 1) = v_ref; 
        X_hat_centr(idx_v, 1) = v_ref;
    end
    
    X_hat_z(:, 1) = obs_dist.T_global \ X_hat_z(:, 1);
    
    % --- Leader initialisation ---
    v_leader = v_ref; 
    v_target = 20;
    
    for k = 1 : num_steps - 1
        t = (k-1) * p.dt; % Time starts at 0
        
        % 1. Leader acceleration (a_0)
        a_leader = 0;
        if strcmp(maneuver, 'cruise')
        if v_leader < v_target
            a_leader = 1.5; % Accelerate
        else
            a_leader = 0;   % Cruise
        end        
        end
        if strcmp(maneuver, 'braking')
            if t >= 5 && t <= 10
                a_leader = -3; % Hard brake
            end
        end
        
        % 2. Leader velocity update (v_0)
        v_leader_next = v_leader + a_leader * p.dt; 
        
        % 3. Input vector U_k = [v_0; a_0]
        U_k = [v_leader; a_leader];
        
        % 4. Plant evolution
        X(:, k+1) = A * X(:, k) + B * U_k;
        
        % 5. Measurement (with outliers)
        y_meas = C * X(:, k) + dist_hist(:, k);
        
        % 6. Decentralized Stubborn Observer
        [X_hat_z(:, k+1), Sigma(:, k+1), Delta(:, k+1)] = stubborn_finitememory(...
            X_hat_z(:, k), Sigma(:, k), y_meas, U_k, obs_dist);
            
        % 7. Centralized Stubborn Observer
        [X_hat_centr(:, k+1), Sigma_centr(:, k+1), Delta_centr(:, k+1)] = observer_stubborn_finite_memory(...
            X_hat_centr(:, k), Sigma_centr(:, k), y_meas, U_k, obs_cent);
            
        % Leader velocity for the next step
        v_leader = v_leader_next;
    end
    
    % Convert DSO back to physical coordinates
    X_hat_phys = obs_dist.T_global * X_hat_z;
end
