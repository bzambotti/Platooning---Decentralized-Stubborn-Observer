function [X, X_hat_phys, X_hat_centr, Delta, Delta_centr] = run_simulation(p, A, B, C, obs_dist, obs_cent, num_steps, dist_hist, maneuver)

    X = zeros(p.n, num_steps);
    X_hat_z = zeros(p.n, num_steps);
    Sigma = zeros(size(obs_dist.S,1), num_steps);
    Delta = zeros(size(obs_dist.R,1), num_steps); % Memoria Distribuita
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

    % --- INIZIALIZZAZIONE LEADER ---
    v_leader = v_ref;
    v_target = 20;

    for k = 1 : num_steps - 1
        t = (k-1) * p.dt; % Time starts at 0

        % 1. Calcolo accelerazione del leader (a_0)
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

        % 2. Aggiornamento velocità del leader (v_0)
        v_leader_next = v_leader + a_leader * p.dt;

        % 3. Costruzione del vettore di input U_k (2x1)
        U_k = [v_leader; a_leader];

        % 4. Evolution del sistema fisico
        X(:, k+1) = A * X(:, k) + B * U_k;

        % 5. Measurement (con rumore/outlier)
        y_meas = C * X(:, k) + dist_hist(:, k);

        % 6. DSO (Nota: passiamo U_k invece di u_scalar)
        [X_hat_z(:, k+1), Sigma(:, k+1), Delta(:, k+1)] = stubborn_finitememory(...
            X_hat_z(:, k), Sigma(:, k), y_meas, U_k, obs_dist);

        % 7. CSO (Nota: passiamo U_k invece di u_scalar)
        [X_hat_centr(:, k+1), Sigma_centr(:, k+1), Delta_centr(:, k+1)] = observer_stubborn_finite_memory(...
            X_hat_centr(:, k), Sigma_centr(:, k), y_meas, U_k, obs_cent);

        % Prepariamo la velocità del leader per il prossimo step
        v_leader = v_leader_next;
    end

    % Convert DSO back to physical coordinates
    X_hat_phys = obs_dist.T_global * X_hat_z;
end
