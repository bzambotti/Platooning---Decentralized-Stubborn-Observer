% Partition and matrices transformation (Algorithm 1)


function [Pd, n_part, m_part, At, Bt, Ct, Lt, T_global] = findAndTransformObserver(A, B, C, L)
    n = size(A,1);
    P = decentralizeOutputs(A, C);      % P_init   
    Pd = {};                            % Final partition (list of sensor indices)
    n_part = [];                        % Vector of state dimensions per block
    m_part = [];                        % Vector of output dimensions per block
    n_pr = 0;                           % States already processed


    At = A; Bt = B; Ct = C; Lt = L;
    T_global = eye(n);

    l = 1; 

    while ~isempty(P)
        found = false;

        if l > numel(P)
             % If we find nothing, we take the rest as a residual block
             idx_rest = [P{:}];
             Pd{end+1} = idx_rest;

             % Calculate dimensions of the residual block
             m_current = length(idx_rest);
             n_current = n - n_pr; % All remaining states

             m_part = [m_part, m_current];
             n_part = [n_part, n_current];
             break; 
        end

        comb = nchoosek(1:numel(P), l);

        for k = 1:size(comb,1)
            idx_in_P = comb(k,:);
            candidate_P = P(idx_in_P); 
            R = [candidate_P{:}];   

            Au = At(n_pr+1:n, n_pr+1:n);
            Cu = Ct(R, n_pr+1:n);

            % Observability decomposition
            [T_bar, nu] = structured_obs_decomp(Au, Cu);

            if nu > 0
                Tk = blkdiag(eye(n_pr), T_bar);

                At_test = Tk \ At * Tk;
                Ct_test = Ct * Tk;
                Lt_test = Tk \ Lt;

                %Check Stubborn Compatibility
                compatible = true;
                Lt_L = Lt_test(:, R); 

                remaining_idx = setdiff(1:numel(P), idx_in_P);
                for rem_k = remaining_idx
                    rows_rem = P{rem_k};
                    if norm(Ct_test(rows_rem, :) * Lt_L, 'fro') > 1e-6
                        compatible = false; break;
                    end
                end 

                % Check Stubborn Compatibility (NUOVA VERSIONE - Zeri Strutturali in L)

                if compatible
                    % --- BLOCK ACCEPTED ---
                    Pd{end+1} = R;
                    P(idx_in_P) = []; 

                    % SAVE DIMENSIONS
                    n_current = nu;         % Number of observable states of the block
                    m_current = length(R);  % Number of outputs of the block

                    n_part = [n_part, n_current];
                    m_part = [m_part, m_current];

                    n_pr = n_pr + nu;

                    % Matrix Update
                    At = At_test; 
                    Ct = Ct_test; 
                    Lt = Lt_test;
                    T_global = T_global * Tk; 

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

    Bt = T_global \ B;

    % Numerical cleanup
    tol = 1e-4;
    At(abs(At)<tol) = 0; Bt(abs(Bt)<tol) = 0;
    Ct(abs(Ct)<tol) = 0; Lt(abs(Lt)<tol) = 0;
    Pd = flip(Pd);
end

function [T, nu] = structured_obs_decomp(A, C)
    Ob = obsv(A, C);
    nu = rank(Ob, 1e-6); 
    n = size(A, 1);
    if nu == n, T = eye(n); elseif nu == 0, T = eye(n);
    else, [U, ~, ~] = svd(Ob'); T = [U(:, 1:nu), U(:, nu+1:end)]; end
end

function [output_partition] = decentralizeOutputs(A, C)
    [m, ~] = size(C); rem = 1:m; output_partition = {};
    while ~isempty(rem)
        curr = rem(1); span = get_span(A, C(curr,:)); blk = curr; rem(1)=[]; i=1;
        while i<=length(rem)
            next=rem(i); span_c = get_span(A, [C(curr,:); C(next,:)]);
            if size(span_c,2)==size(span,2) && norm(span_c-span,'fro')<1e-9, blk=[blk, next]; rem(i)=[]; else, i=i+1; end
        end
        output_partition{end+1}=blk;
    end
    num_blocks = length(output_partition);
end

function span = get_span(A, C)
    Ob = obsv(A,C); r = rank(Ob, 1e-9); [~,~,V] = svd(Ob, 'econ'); span = V(:, 1:r);
end
