function L_dense = generate_observer_gain(Q_weight, R_weight, p)
% Generates 

A = p.A;
C = p.C;
L_dense = zeros(size(A, 1), size(C, 1));

for i = 1:p.N
    % Indici correnti (Ego)
    idx_x_i = (i-1)*p.n_veh + 1 : i*p.n_veh; 
    idx_y_i = (i-1)*2 + 1 : i*2; 
    
    if i == 1
       % --- PRIMO FOLLOWER (i=1) ---
       % Non ha predecessore nello stato X.
        A_11 = A(idx_x_i, idx_x_i);
        C_11 = C(idx_y_i, idx_x_i);
        
        % LQR Locale
        L_loc = dlqr(A_11', C_11', Q_weight*eye(p.n_veh), R_weight*eye(2))';
        
        % Assegna al blocco (1,1)
        L_dense(idx_x_i, idx_y_i) = L_loc;
        
    else 
        % --- ALTRI FOLLOWERS (FOLLOWER) ---
        % Consideriamo il sottosistema "Coppia": [Predecessore; Ego]
        % Indici Predecessore (Prev)
        idx_x_prev = (i-2)*p.n_veh + 1 : (i-1)*p.n_veh;
        idx_y_prev = (i-2)*2 + 1 : (i-1)*2;
        
        % Costruiamo le matrici del sistema AUMENTATO (2 veicoli)
        % A_aug = [A_prev_prev,     0; 
        %          A_ego_prev,  A_ego_ego]
        A_pp = A(idx_x_prev, idx_x_prev);
        A_ee = A(idx_x_i, idx_x_i);
        A_ep = A(idx_x_i, idx_x_prev); % Accoppiamento fisico
        
        A_aug = [A_pp, zeros(p.n_veh); 
                 A_ep, A_ee];
             
        % C_aug è diagonale a blocchi (ogni veicolo ha i suoi sensori)
        C_pp = C(idx_y_prev, idx_x_prev);
        C_ee = C(idx_y_i, idx_x_i);
        C_aug = blkdiag(C_pp, C_ee);
        
        % Q e R aumentati
        Q_aug = blkdiag(Q_weight*eye(p.n_veh), Q_weight*eye(p.n_veh));
        R_aug = blkdiag(R_weight*eye(2), R_weight*eye(2));
        
        % Calcolo LQR distribuito per la coppia
        L_aug = dlqr(A_aug', C_aug', Q_aug, R_aug)';
        
        % --- ESTRAZIONE DEI GUADAGNI ---
        % L_aug ha dimensione (2*n_veh) x (4 output). 
        % A noi interessa solo come aggiornare lo stato di EGO (la metà inferiore delle righe)
        % usando le misure di PREV e EGO.
        
        % Righe relative a Ego nel sistema aumentato: da (n_veh+1) a (2*n_veh)
        rows_ego_in_aug = p.n_veh + 1 : 2*p.n_veh;
        
        % Guadagno L_{i, i-1} (Quanto ascolto il predecessore)
        L_cross = L_aug(rows_ego_in_aug, 1:2); 
        
        % Guadagno L_{i, i} (Quanto ascolto me stesso)
        L_self  = L_aug(rows_ego_in_aug, 3:4);
        
        % Assegna alla matrice globale
        L_dense(idx_x_i, idx_y_prev) = L_cross; % Blocco diagonale inferiore
        L_dense(idx_x_i, idx_y_i)    = L_self;  % Blocco diagonale principale
    end
end

end