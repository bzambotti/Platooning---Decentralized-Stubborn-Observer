function [X, X_hat_phys, X_hat_centr] = run_simulation(p, A, B, C, obs_dist, obs_cent, num_steps, dist_hist, maneuver)
    X = zeros(p.n, num_steps);
    X_hat_z = zeros(p.n, num_steps);
    Sigma = zeros(size(obs_dist.S,1), num_steps);
    
    X_hat_centr = zeros(p.n, num_steps);
    Sigma_centr = zeros(size(obs_cent.S,1), num_steps);
    
    % Initial Conditions
    % Assuming absolute velocity for simulation realism in Scenario 3.
    v_ref = 20;
    for i = 1:p.N
        idx_v = (i-1)*p.n_veh + 2;
        X(idx_v, 1) = v_ref;
        X_hat_z(idx_v, 1) = v_ref; 
        X_hat_centr(idx_v, 1) = v_ref;
    end
    
    % Transform initial estimate to z (for decentralized)
    % Note: X_hat_z is in transformed coordinates.
    % We assume X_hat_phys initial is same as X.
    % z = T_global \ x
    X_hat_z(:, 1) = obs_dist.T_global \ X_hat_z(:, 1);

    for k = 1 : num_steps - 1
        % Input
        u_scalar = 0;
        if strcmp(maneuver, 'braking')
            t = (k-1) * p.dt; % Time starts at 0
            if t >= 10 && t <= 15
                u_scalar = -4; % Hard brake
            end
        end
        
        % Evolution
        X(:, k+1) = A * X(:, k) + B * u_scalar;
        
        % Measurement
        y_meas = C * X(:, k) + dist_hist(:, k);
        
        % DSO
        [X_hat_z(:, k+1), Sigma(:, k+1)] = stubborn_finitememory(...
            X_hat_z(:, k), Sigma(:, k), y_meas, u_scalar, obs_dist);
            
        % CSO
        [X_hat_centr(:, k+1), Sigma_centr(:, k+1)] = observer_stubborn_finite_memory(...
            X_hat_centr(:, k), Sigma_centr(:, k), y_meas, u_scalar, obs_cent);
    end
    
    % Convert DSO back to physical coordinates
    X_hat_phys = obs_dist.T_global * X_hat_z;
end
