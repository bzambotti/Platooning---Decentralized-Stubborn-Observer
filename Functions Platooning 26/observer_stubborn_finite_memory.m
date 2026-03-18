function [x_hat_next, Sigma_next, Delta_k] = observer_stubborn_finite_memory(x_hat_k, Sigma_k, y_meas, u_k, params)
    % CORREZIONE: L'ordine degli input ora corrisponde alla chiamata nel main:
    % 1. x_hat_k (Stato)
    % 2. Sigma_k (Memoria)
    % 3. y_meas  (Misura)
    % 4. u_k     (Input scalare)
    % 5. params  (Struttura)

    % Estrai matrici sistema e osservatore
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
    % Ricostruiamo le matrici S e G al volo (o potresti passarle in params)
    S = diag(ones(tau-1, 1), -1);
    G = zeros(tau, 1);
    G(1) = 1;
    
    Sigma_next = S * Sigma_k + G * abs_y_err;
end