
function [P_star, N_star, T_t, At, Ct, Lt, T] = Algorithm1_DecentralizedPartition(A, C, L)
    % Implementazione 1:1 dell'Algorithm 1
    P0 = decentralizeOutputs(A, C)
    n = size(A, 1);
    
    % Line 1: Initialize
    P_star = {};  
    N_star = [];
    
    % Line 2: Initialize intermediate matrices and counters
    At = A; Ct = C; Lt = L; 
    T_t = eye(n); 
    n_pr = 0; 
    ell = 1;
    
    % Line 3: While P_0 is not empty
    while ~isempty(P0)
        found = false; % Line 4
        
        % Line 6: For all candidate blocks L subset of P_0 with |L| = ell
        comb = nchoosek(1:numel(P0), ell);
        for k = 1:size(comb, 1)
            idx_in_P0 = comb(k, :);
            candidate_L = P0(idx_in_P0); 
            
            % Line 8 & 9: Extract unprocessed subsystem outputs
            I_L = [candidate_L{:}]; 
            
            % Line 10: Unprocessed state indices (TOP-LEFT extraction)
            I_u = 1 : (n - n_pr); 
            
            % Line 11: Extract submatrices
            Au = At(I_u, I_u);
            Cu = Ct(I_L, I_u);
            
            % Line 13 & 14: Observable decomposition via Hespanha
            [T_L, n_o] = structured_obs_decomp(Au, Cu);
            
            % Line 15: Set T_c protecting the bottom-right already processed block
            Tc = blkdiag(T_L, eye(n_pr));
            
            % Line 17 & 18: Apply candidate transformation
            Ac = Tc \ At * Tc;
            Cc = Ct * Tc;
            Lc = Tc \ Lt;
            
            % Line 20: Compatibility test setup
            compatible = true;
            L_L = Lc(:, I_L); % Columns of Lc indexed by L
            
            % Line 23: For all blocks B in P_0 \ {L}
            remaining_idx = setdiff(1:numel(P0), idx_in_P0);
            for rem_k = remaining_idx
                I_B = P0{rem_k}; % Row indices corresponding to B
                
                % Line 25 & 26: If C_c * L_L != 0
                if norm(Cc(I_B, :) * L_L, 'fro') > 1e-6
                    compatible = false; 
                    break;
                end
            end
            
            % Line 30: If compatible
            if compatible
                At = Ac; Ct = Cc; Lt = Lc;          % Accept transformation
                P_star{end+1} = I_L;                % Append to P*
                N_star = [N_star, n_o];             % Append to N*
                P0(idx_in_P0) = [];                 % Remove L from P0
                T_t = T_t * Tc;                     % Accumulate transformation
                n_pr = n_pr + n_o;                  % Update processed dimension
                ell = 1;                            % Reset candidate size
                found = true;
                break;
            end
        end
        
        % Line 37: If not found
        if ~found
            ell = ell + 1; % Increase candidate block size
            
            % Line 39: Worst-case, group remaining outputs
            if ell > numel(P0)
                I_L = [P0{:}];
                P_star{end+1} = I_L;
                N_star = [N_star, (n - n_pr)]; % The remaining states form the last block
                break;
            end
        end
    end

    T = buildFinalTransformations(T_t, N_star, P_star);
end

function [T, Pi_s, Pi_y] = buildFinalTransformations(T_t, N_star, P_star)
    % BUILDFINALTRANSFORMATIONS costruisce le matrici di permutazione e la
    % trasformazione globale per forzare la struttura triangolare inferiore.
    
    % Ricaviamo le dimensioni globali direttamente dagli input
    n = size(T_t, 1);
    m = length([P_star{:}]);
    
    % 1. Ribaltiamo le sequenze (da Bottom-Up a Top-Down causale)
    P_reversed = flip(P_star);
    N_reversed = flip(N_star); 
    
    % =====================================================================
    % 2. COSTRUZIONE DI Pi_y (Permutazione delle Uscite)
    % =====================================================================
    % Estraiamo gli indici dei sensori concatenando la partizione ribaltata
    I_y = [P_reversed{:}];
    
    % Pi_y è una matrice identità le cui righe sono state riordinate
    Pi_y = eye(m);
    Pi_y = Pi_y(I_y, :);
    
    % =====================================================================
    % 3. COSTRUZIONE DI Pi_s (Permutazione degli Stati)
    % =====================================================================
    % Raggruppiamo gli indici di stato nei blocchi originali definiti 
    % da N_star, per poi incollarli al contrario.
    
    I_s_blocks = cell(1, length(N_star));
    current_idx = 1;
    
    % Creiamo i blocchi di indici basati sull'ordine originale (Bottom-Up)
    for i = 1:length(N_star)
        I_s_blocks{i} = current_idx : (current_idx + N_star(i) - 1);
        current_idx = current_idx + N_star(i);
    end
    
    % Invertiamo l'ordine dei blocchi di stato e li concateniamo
    I_s_blocks_reversed = flip(I_s_blocks);
    I_s = [I_s_blocks_reversed{:}];
    
    % Pi_s è una matrice identità le cui righe sono state riordinate
    Pi_s = eye(n);
    Pi_s = Pi_s(I_s, :);
    
    % =====================================================================
    % 4. COSTRUZIONE DI T (Trasformazione Globale)
    % =====================================================================

    T = T_t * Pi_s';
end

function [output_partition] = decentralizeOutputs(A, C)
% Helper function to speed the algorithm: looks for outputs with the same
% observable subspace
    [m, ~] = size(C); rem = 1:m; output_partition = {};
    while ~isempty(rem)
        curr = rem(1); span = get_span(A, C(curr,:)); blk = curr; rem(1)=[]; i=1;
        while i<=length(rem)
            next=rem(i); span_c = get_span(A, [C(curr,:); C(next,:)]);
            if size(span_c,2)==size(span,2) && norm(span_c-span,'fro')<1e-9, blk=[blk, next]; rem(i)=[]; else, i=i+1; end
        end
        output_partition{end+1}=blk;
    end
end

function span = get_span(A, C)
    Ob = obsv(A,C); r = rank(Ob, 1e-9); [~,~,V] = svd(Ob, 'econ'); span = V(:, 1:r);
end

function [T, nu] = structured_obs_decomp(A, C)
    Ob = obsv(A, C);
    nu = rank(Ob, 1e-6); 
    n = size(A, 1);
    if nu == n, T = eye(n); elseif nu == 0, T = eye(n);
    else, [U, ~, ~] = svd(Ob'); T = [U(:, 1:nu), U(:, nu+1:end)]; end
end
