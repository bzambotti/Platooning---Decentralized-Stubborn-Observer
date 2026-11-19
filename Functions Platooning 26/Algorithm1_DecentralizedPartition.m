function [A_tilde, C_tilde, L_tilde, P_star, N_star, T, Pi_y, Pi_s] = Algorithm1_DecentralizedPartition(A, C, L, tol)
    % BLOCK_TRIANGULARIZATION: Algorithm 1 
    %
    % Inputs:
    %   A, C, L : System matrices
    %   tol     : Numerical tolerance for compatibility test (default: 1e-8)

    if nargin < 4
        tol = 1e-8;
    end

    n = size(A, 1);
    m = size(C, 1);

    % Initialize
    P0 = decentralizeOutputs(A, C);
    P_star = {};
    N_star = [];

    At = A; Ct = C; Lt = L;
    T_t = eye(n);
    n_pr = 0;

    while ~isempty(P0)
        found = false;
        l = 1; % Start with single outputs

        while l <= length(P0) && ~found
            % Generate all combinations of candidate blocks of size l
            combos = nchoosek(1:length(P0), l);

            for idx = 1:size(combos, 1)
                L_idx = combos(idx, :);
                candidate_L = [P0{L_idx}]; % Output indices for L

                I_u = 1:(n - n_pr);
                Au = At(I_u, I_u);
                Cu = Ct(candidate_L, I_u);

                % Observable decomposition 
                [~, ~, ~, T_mat, k] = obsvf(Au, zeros(length(I_u),length(candidate_L)), Cu, tol);
                no = sum(k);
                T_L = T_mat'; 
                Tc = blkdiag(T_L, eye(n_pr));

                % Apply candidate transformation
                Ac = Tc \ (At * Tc);
                Cc = Ct * Tc;
                Lc = Tc \ Lt;

                % Compatibility test
                compatible = true;
                L_L = Lc(:, candidate_L);

                remaining_idx = setdiff(1:length(P0), L_idx);
                for r_idx = remaining_idx
                    B_indices = P0{r_idx};
                    % Check if innovation injects into remaining blocks
                    if norm(Cc(B_indices, :) * L_L, 'inf') > tol
                        compatible = false;
                        break;
                    end
                end




                if compatible
                    At = Ac; Ct = Cc; Lt = Lc;
                    P_star{end+1} = candidate_L;
                    N_star(end+1) = no;

                    P0(L_idx) = []; % Remove accepted block
                    T_t = T_t * Tc;
                    n_pr = n_pr + no;

                    found = true;
                    break; % Exit combo loop
                end
            end

            if ~found
                l = l + 1;
                % Worst-case scenario: group remaining outputs
                if l > length(P0)
                    candidate_L = [P0{:}];
                    P_star{end+1} = candidate_L;

                    % Extract final observable part to isolate A_star_star
                    I_u = 1:(n - n_pr);
                    Au = At(I_u, I_u);
                    Cu = Ct(candidate_L, I_u);
                    [~, ~,~, T_mat, k] = obsvf(Au, zeros(length(I_u),length(candidate_L)), Cu, tol);
                    no = sum(k);

                    Tc = blkdiag(T_mat', eye(n_pr));
                    At = Tc \ (At * Tc);
                    Ct = Ct * Tc;
                    Lt = Tc \ Lt;
                    T_t = T_t * Tc;

                    N_star(end+1) = no;
                    P0 = {};
                    break;
                end
            end
        end
    end

    % --- Post-Processing: Reversing the order to achieve Property 1 ---

    % Output permutation Pi_y
    idx_y = [];
    for i = length(P_star):-1:1 % Reverse sequence
        idx_y = [idx_y, P_star{i}];
    end
    Pi_y = eye(m);
    Pi_y = Pi_y(idx_y, :);

    % State permutation Pi_s
    n_bar_o = n - sum(N_star); % Fully unobservable part
    state_blocks = [n_bar_o, fliplr(N_star)]; % Current block sizes from top to bottom

    idx_s = [];
    curr_idx = n; % Start from the bottom of At
    for i = 2:length(state_blocks) % Iterate through N_star (ignoring n_bar_o for now)
        bs = state_blocks(i);
        idx_s = [(curr_idx - bs + 1):curr_idx, idx_s];
        curr_idx = curr_idx - bs;
    end
    idx_s = [idx_s, 1:n_bar_o]; % Put unobservable part at the end (bottom right)

    Pi_s = eye(n);
    Pi_s = Pi_s(idx_s, :);

    % Final transformations
    T = T_t * Pi_s';
    A_tilde = T \ (A * T);
    C_tilde = Pi_y * C * T;
    L_tilde = T \ (L * Pi_y');

    % Numerical cleanup (zeros out tiny floating point artifacts)
    A_tilde(abs(A_tilde) < 1e-6) = 0;
    C_tilde(abs(C_tilde) < 1e-6) = 0;
    L_tilde(abs(L_tilde) < 1e-6) = 0;
end


% function [A_tilde, C_tilde, L_tilde, P_star, N_star, T] = ...
%     Algorithm1_DecentralizedPartition(A, C, L, tol)
% 
% if nargin < 4
%     tol = 1e-9;
% end
% 
% n = size(A,1);
% m = size(C,1);
% 
% % --- inizializzazione
% A_bar = A;
% C_bar = C;
% L_bar = L;
% 
% T = eye(n);
% 
% P0 = flip(decentralizeOutputs(A,C));   % prima partizione grossolana
% P_star = {};
% N_star = [];
% 
% n_pr = 0;
% 
% while ~isempty(P0)
% 
%     found = false;
% 
%     for l = 1:length(P0)
% 
%         combos = nchoosek(1:length(P0),l);
% 
%         for k = 1:size(combos,1)
% 
%             idx = combos(k,:);
%             candidate = [P0{idx}];
% 
%             % sottosistema ancora non estratto
%             Iu = 1:(n - n_pr);
%             Au = A_bar(Iu,Iu);
%             Cu = C_bar(candidate,Iu);
% 
%             % decomposizione osservabile 
%             [~,~,~,Tloc,kobs] = obsvf(Au,zeros(size(Au,1),length(candidate)),Cu,tol);
%             no = sum(kobs);
% 
%             Tloc = Tloc';
%             Tc = blkdiag(Tloc, eye(n_pr));
% 
%             % trasformazione candidata
%             Atest = Tc \ (A_bar * Tc);
%             Ctest = C_bar * Tc;
%             Ltest = Tc \ L_bar;
% 
%             % ---- test di compatibilità ----
%             compatible = true;
%             Lcand = Ltest(:,candidate);
% 
%             remaining = setdiff(1:length(P0),idx);
% 
%             for r = remaining
%                 ind = P0{r};
%                 if norm(Ctest(ind,:) * Lcand,'inf') > tol
%                     compatible = false;
%                     break
%                 end
%             end
% 
%             if compatible
%                 disp(candidate)
% disp(no)
% disp(norm(Atest(1:(n-n_pr-no),(n-n_pr-no+1):(n-n_pr))))
%                 A_bar = Atest;
%                 C_bar = Ctest;
%                 L_bar = Ltest;
% 
%                 T = T * Tc;
% 
%                 P_star{end+1} = candidate;
%                 N_star(end+1) = no;
% 
%                 P0(idx) = [];
%                 n_pr = n_pr + no;
% 
%                 found = true;
%                 break
%             end
%         end
% 
%         if found
%             break
%         end
%     end
% 
%     % caso limite: prendo tutto insieme
%     if ~found
% 
%         candidate = [P0{:}];
% 
%         Iu = 1:(n - n_pr);
%         Au = A_bar(Iu,Iu);
%         Cu = C_bar(candidate,Iu);
% 
%         [~,~,~,Tloc,kobs] = obsvf(Au,zeros(size(Au,1),length(candidate)),Cu,tol);
%         no = sum(kobs);
% 
%         Tc = blkdiag(Tloc',eye(n_pr));
% 
%         A_bar = Tc \ (A_bar * Tc);
%         C_bar = C_bar * Tc;
%         L_bar = Tc \ L_bar;
%         T = T * Tc;
% 
%         P_star{end+1} = candidate;
%         N_star(end+1) = no;
% 
%         break
%     end
% end
% 
% A_tilde = A_bar;
% C_tilde = C_bar;
% L_tilde = L_bar;
% 
% % --- Post-Processing: Reversing the order to achieve Property 1 ---
% 
%     % Output permutation Pi_y
%     idx_y = [];
%     for i = length(P_star):-1:1 % Reverse sequence
%         idx_y = [idx_y, P_star{i}];
%     end
%     Pi_y = eye(m);
%     Pi_y = Pi_y(idx_y, :);
% 
%     % State permutation Pi_s
%     n_bar_o = n - sum(N_star); % Fully unobservable part
%     state_blocks = [n_bar_o, fliplr(N_star)]; % Current block sizes from top to bottom
% 
%     idx_s = [];
%     curr_idx = n; % Start from the bottom of At
%     for i = 2:length(state_blocks) % Iterate through N_star (ignoring n_bar_o for now)
%         bs = state_blocks(i);
%         idx_s = [(curr_idx - bs + 1):curr_idx, idx_s];
%         curr_idx = curr_idx - bs;
%     end
%     idx_s = [idx_s, 1:n_bar_o]; % Put unobservable part at the end (bottom right)
% 
%     Pi_s = eye(n);
%     Pi_s = Pi_s(idx_s, :);
% 
%     % Final transformations
%     T = T * Pi_s';
%     A_tilde = T \ (A * T);
%     C_tilde = Pi_y * C * T;
%     L_tilde = T \ (L * Pi_y');
% 
%     % Numerical cleanup (zeros out tiny floating point artifacts)
%     A_tilde(abs(A_tilde) < 1e-6) = 0;
%     C_tilde(abs(C_tilde) < 1e-6) = 0;
%     L_tilde(abs(L_tilde) < 1e-6) = 0;
% 
% end



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
