function [x_hat_next, Sigma_next, Delta_k] = observer_stubborn_finite_memory(x_hat_k, Sigma_k, y_meas, u_k, params)
    % Centralized finite-memory stubborn observer (single global threshold).
    % Inputs: estimate x_hat_k, memory buffer Sigma_k, measurement y_meas,
    %         input u_k, params (A, B, C, L, D, Theta, tau).

    % System and observer matrices
    A = params.A; 
    B = params.B;
    C = params.C;
    L = params.L;
    
    % Stubborn parameters
    Theta = params.Theta;
    tau = params.tau;
    
    % Output error (tilde_y_k)
    y_hat = C * x_hat_k;
    
    % Feedthrough term, if present
    if isfield(params, 'D') && ~isempty(params.D)
         y_hat = y_hat + params.D * u_k;
    end
    
    tilde_y = y_meas - y_hat;
    
    % Memory input (scalar in the centralized case)
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
    x_hat_next = A * x_hat_k + B * u_k + L * sat_tilde_y;
    
    % Threshold Update
    % Shift-register matrices S and G
    S = diag(ones(tau-1, 1), -1);
    G = zeros(tau, 1);
    G(1) = 1;
    
    Sigma_next = S * Sigma_k + G * abs_y_err;
end