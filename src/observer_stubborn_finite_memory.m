function [x_hat_next, Sigma_next, Delta_k] = observer_stubborn_finite_memory(x_hat_k, Sigma_k, y_meas, u_k, params)
    % OBSERVER_STUBBORN_FINITE_MEMORY Updates the centralized finite memory stubborn observer.
    % Computes the next centralized estimation using the current measurement,
    % input, and dynamic saturation threshold computed from its finite memory buffer.
    %
    % Inputs:
    % 1. x_hat_k (Current State Estimation)
    % 2. Sigma_k (Current Error Memory Buffer)
    % 3. y_meas  (Current Measurement)
    % 4. u_k     (Control Input)
    % 5. params  (Observer and System Parameter Struct)

    % Extract System and Observer Matrices
    A = params.A; 
    B = params.B;
    C = params.C;
    L = params.L;
    
    % Parametri Stubborn
    Theta = params.Theta;
    tau = params.tau;
    
    % Output error (tilde_y_k)
    y_hat = C * x_hat_k;
    
    % Aggiunta robustezza: gestiamo D se presente (anche se zero)
    if isfield(params, 'D') && ~isempty(params.D)
         y_hat = y_hat + params.D * u_k;
    end
    
    tilde_y = y_meas - y_hat;
    
    % Input per la memoria (SCALARE per il caso centralizzato)
    abs_y_err = norm(tilde_y, 2); 
    
    % Saturation Threshold (Delta_k)
    % Delta_k = Theta * sum(Sigma_k) equivale a R * Sigma se R = [Theta Theta ...]
    Delta_k = Theta * sum(Sigma_k);
    
    % Error Saturation
    sat_tilde_y = zeros(size(tilde_y));
    for i = 1:length(tilde_y)
        val = tilde_y(i);
        if abs(val) > Delta_k
            sat_tilde_y(i) = sign(val) * Delta_k;
        else
            sat_tilde_y(i) = val;
        end
    end
    
    % Observer Update
    % Ora u_k è uno scalare, quindi B * u_k funziona (n x 1) * (1 x 1)
    x_hat_next = A * x_hat_k + B * u_k + L * sat_tilde_y;
    
    % Threshold Update
    % Extracting S and G matrices directly from params to avoid re-instantiating
    % them at every simulation step, optimizing execution speed.
    S = params.S;
    G = params.G;
    
    Sigma_next = S * Sigma_k + G * abs_y_err;
end