% Partition and matrices transformation (Algorithm 1).

function [Pd, n_part, m_part, At, Bt, Ct, Lt, T_global] = findAndTransformObserver(A, B, C, L)
    n = size(A,1);
    P = decentralizeOutputs(A, C);      % P_init   
    Pd = {};                            % Final partition

    At_temp = A; Ct_temp = C; Lt_temp = L;
    n_pr = 0; 
    l = 1; 

    while ~isempty(P)
        found = false;
        if l > numel(P)
             idx_rest = [P{:}];
             Pd{end+1} = idx_rest;
             break; 
        end
        comb = nchoosek(1:numel(P), l);
        for k = 1:size(comb,1)
            idx_in_P = comb(k,:);
            candidate_P = P(idx_in_P); 
            R = [candidate_P{:}];   
            Au = At_temp(n_pr+1:n, n_pr+1:n);
            Cu = Ct_temp(R, n_pr+1:n);

            [T_bar, nu] = structured_obs_decomp(Au, Cu);
            if nu > 0
                Tk = blkdiag(eye(n_pr), T_bar);
                At_test = Tk \ At_temp * Tk;
                Ct_test = Ct_temp * Tk;
                Lt_test = Tk \ Lt_temp;

                % Check Compatibility
                compatible = true;
                Lt_L = Lt_test(:, R); 
                remaining_idx = setdiff(1:numel(P), idx_in_P);
                for rem_k = remaining_idx
                    rows_rem = P{rem_k};
                    if norm(Ct_test(rows_rem, :) * Lt_L, 'fro') > 1e-6
                        compatible = false; break;
                    end
                end 

                if compatible
                    Pd{end+1} = R;
                    P(idx_in_P) = []; 
                    n_pr = n_pr + nu;
                    At_temp = At_test; 
                    Ct_temp = Ct_test; 
                    Lt_temp = Lt_test;
                    found = true;
                    l = 1; 
                    break; 
                end
            end
        end
        if ~found
            l = l + 1; 
        end
    end


    Pd = flip(Pd); 
    At = A; Ct = C; Lt = L; 
    T_global = eye(n);
    n_pr = 0;
    n_part = []; m_part = [];

    % Riordino righe di C e colonne di L in base alla nuova partizione
    idx_outputs_permuted = [Pd{:}];
    Ct = Ct(idx_outputs_permuted, :);
    Lt = Lt(:, idx_outputs_permuted);

    current_row_idx = 1; % per le righe di Ct

    for i = 1:length(Pd)
        block_outputs = Pd{i};
        m_current = length(block_outputs);
        m_part = [m_part, m_current];

        % parte non processata
        Au = At(n_pr+1:n, n_pr+1:n);
        Cu = Ct(current_row_idx : current_row_idx + m_current - 1, n_pr+1:n);

        [T_bar, nu] = structured_obs_decomp(Au, Cu);

        n_current = nu;
        n_part = [n_part, n_current];

        % Trasformazione
        Tk = blkdiag(eye(n_pr), T_bar);
        At = Tk \ At * Tk;
        Ct = Ct * Tk;
        Lt = Tk \ Lt;
        T_global = T_global * Tk;

        n_pr = n_pr + nu;
        current_row_idx = current_row_idx + m_current;
    end

    Bt = T_global \ B;

    % Numerical cleanup
    tol = 1e-6; 
    At(abs(At)<tol) = 0; Bt(abs(Bt)<tol) = 0;
    Ct(abs(Ct)<tol) = 0; Lt(abs(Lt)<tol) = 0;
end


function [T, nu] = structured_obs_decomp(A, C)
    Ob = obsv(A, C);
    nu = rank(Ob, 1e-6); 
    n = size(A, 1);
    if nu == n, T = eye(n); elseif nu == 0, T = eye(n);
    else, [U, ~, ~] = svd(Ob'); T = [U(:, 1:nu), U(:, nu+1:end)]; end
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


% function [Pd, n_part, m_part, At, Bt, Ct, Lt, T_global] = findAndTransformObserver(A, B, C, L)
%     n = size(A,1);
%     P = decentralizeOutputs(A, C);      % P_init   
%     Pd = {};                            % Final partition
% 
%     At_temp = A; Ct_temp = C; Lt_temp = L;
%     n_pr = 0; 
%     l = 1; 
% 
%     % =========================================================
%     % PHASE 1: TOP-DOWN PARTITION DISCOVERY
%     % =========================================================
%     while ~isempty(P)
%         found = false;
%         if l > numel(P)
%              % Fallback: group remaining sensors
%              idx_rest = [P{:}];
%              Pd{end+1} = idx_rest;
%              break; 
%         end
%         comb = nchoosek(1:numel(P), l);
% 
%         for k = 1:size(comb,1)
%             idx_in_P = comb(k,:);
%             candidate_P = P(idx_in_P); 
%             R = [candidate_P{:}];   
% 
%             Au = At_temp(n_pr+1:n, n_pr+1:n);
%             Cu = Ct_temp(R, n_pr+1:n);
% 
%             [T_bar, nu] = structured_obs_decomp(Au, Cu);
% 
%             if nu > 0
%                 Tk = blkdiag(eye(n_pr), T_bar);
%                 % We only need Lt_test to verify causality
%                 Lt_test = Tk \ Lt_temp;
% 
%                 % ---------------------------------------------------
%                 % NEW COMPATIBILITY TEST: Top-Down Causality on L
%                 % ---------------------------------------------------
%                 compatible = true;
% 
%                 % Indices of the newly isolated observable states
%                 new_states_idx = (n_pr + 1) : (n_pr + nu);
% 
%                 remaining_idx = setdiff(1:numel(P), idx_in_P);
%                 if ~isempty(remaining_idx)
%                     remaining_sensors = [P{remaining_idx}];
% 
%                     % Extract the submatrix L_{12} that MUST be zero
%                     L_check = Lt_test(new_states_idx, remaining_sensors);
% 
%                     if norm(L_check, 'inf') > 0.01
%                         compatible = false; 
%                     end
%                 end 
%                 % ---------------------------------------------------
% 
%                 if compatible
%                     Pd{end+1} = R;
%                     P(idx_in_P) = []; 
%                     n_pr = n_pr + nu;
% 
%                     % Update temporary matrices
%                     At_temp = Tk \ At_temp * Tk; 
%                     Ct_temp = Ct_temp * Tk; 
%                     Lt_temp = Lt_test;
% 
%                     found = true;
%                     l = 1; 
%                     break; 
%                 end
%             end
%         end
%         if ~found
%             l = l + 1; 
%         end
%     end
% 
%     % NO FLIP HERE! Top-Down naturally yields the correct lower-triangular order.
% 
%     % =========================================================
%     % PHASE 2: CLEAN TRANSFORMATION RECONSTRUCTION
%     % =========================================================
%     At = A; Ct = C; Lt = L; 
%     T_global = eye(n);
%     n_pr = 0;
%     n_part = []; m_part = [];
% 
%     % Permute C rows and L columns according to the discovered partition
%     idx_outputs_permuted = [Pd{:}];
%     Ct = Ct(idx_outputs_permuted, :);
%     Lt = Lt(:, idx_outputs_permuted);
% 
%     current_row_idx = 1; 
% 
%     for i = 1:length(Pd)
%         block_outputs = Pd{i};
%         m_current = length(block_outputs);
%         m_part = [m_part, m_current];
% 
%         % Unprocessed subspace
%         Au = At(n_pr+1:n, n_pr+1:n);
%         Cu = Ct(current_row_idx : current_row_idx + m_current - 1, n_pr+1:n);
% 
%         [T_bar, nu] = structured_obs_decomp(Au, Cu);
% 
%         n_current = nu;
%         n_part = [n_part, n_current];
% 
%         % Global Transformation
%         Tk = blkdiag(eye(n_pr), T_bar);
%         At = Tk \ At * Tk;
%         Ct = Ct * Tk;
%         Lt = Tk \ Lt;
%         T_global = T_global * Tk;
% 
%         n_pr = n_pr + nu;
%         current_row_idx = current_row_idx + m_current;
%     end
% 
%     Bt = T_global \ B;
% 
%     % Numerical cleanup
%     tol = 1e-6; 
%     At(abs(At)<tol) = 0; Bt(abs(Bt)<tol) = 0;
%     Ct(abs(Ct)<tol) = 0; Lt(abs(Lt)<tol) = 0;
% end

% Helper functions remain identical
% function [T, nu] = structured_obs_decomp(A, C)
%     Ob = obsv(A, C);
%     nu = rank(Ob, 1e-6); 
%     n = size(A, 1);
%     if nu == n, T = eye(n); elseif nu == 0, T = eye(n);
%     else, [U, ~, ~] = svd(Ob'); T = [U(:, 1:nu), U(:, nu+1:end)]; end
% end

% function [output_partition] = decentralizeOutputs(A, C)
%     [m, ~] = size(C); rem = 1:m; output_partition = {};
%     while ~isempty(rem)
%         curr = rem(1); span = get_span(A, C(curr,:)); blk = curr; rem(1)=[]; i=1;
%         while i<=length(rem)
%             next=rem(i); span_c = get_span(A, [C(curr,:); C(next,:)]);
%             if size(span_c,2)==size(span,2) && norm(span_c-span,'fro')<1e-9, blk=[blk, next]; rem(i)=[]; else, i=i+1; end
%         end
%         output_partition{end+1}=blk;
%     end
% end
% 
% function span = get_span(A, C)
%     Ob = obsv(A,C); r = rank(Ob, 1e-9); [~,~,V] = svd(Ob, 'econ'); span = V(:, 1:r);
% end