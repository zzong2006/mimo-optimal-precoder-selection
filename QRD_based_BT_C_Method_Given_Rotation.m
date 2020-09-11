% Made by Woosung 2020.Jan.14th

function [F, index] = QRD_based_BT_C_Method_Given_Rotation(H, K, delta, perms_array, perms_tree)
% < Input >
% K : The number of the bit streams
% H : Channel, the size of nR X nT
% delta : Orthogonality of the transformed channel matrix H_{LR}
% < Output >
% F : Selected Sub-optimal Precoder
% index : And its index

%% Default Value
if nargin < 1
    disp('[Message] QRD_BT-C_Precoder_Selection  : Default set H, K, delta');
    nT = 4; nR = 4; 
    K = 4;
    H = 1/sqrt(K)*sqrt(1/2)*(randn(nT,nR)+1i*randn(nT,nR));
    delta = 0.99;
elseif nargin < 2
    disp('[Message] QRD_BT-C_Precoder_Selection : Default set K, delta');
    K = 4;
    delta = 0.99;
elseif nargin < 3
    disp('[Message] QRD_BT-C_Precoder_Selection : Default set delta only');
    delta = 0.99;
end


%% Create Precoder based on the LTE-A Codebook
j = sqrt(-1);
u=[ 1,  1,  1,  1,  1,              1,             1,             1,              1,  1,  1,  1,  1,  1,  1, 1;
   -1, -j,  1,  j, (-1-j)/sqrt(2), (1-j)/sqrt(2), (1+j)/sqrt(2), (-1+j)/sqrt(2), -1, -j,  1,  j, -1, -1,  1, 1;
   -1,  1, -1,  1, -j,              j,            -j,             j,              1, -1,  1, -1, -1,  1, -1, 1;
   -1,  j,  1, -j, (1-j)/sqrt(2),  (-1-j)/sqrt(2),(-1+j)/sqrt(2),(1+j)/sqrt(2),   1, -j, -1,  j,  1, -1, -1, 1];

W=zeros(4,4,16);

for i=1:length(W)
    a = u(:, i) * u(:, i)';
    b = u(:, i)' * u(:, i);
    W(:, :, i) = eye(4) - (2 * a) / b;
end    

F4_matrix_order = ...
[[1 2 3 4];[1 2 3 4];[3 2 1 4];[3 2 1 4];
 [1 2 3 4];[1 2 3 4];[1 3 2 4];[1 3 2 4];
 [1 2 3 4];[1 2 3 4];[1 3 2 4];[1 3 2 4];
 [1 2 3 4];[1 3 2 4];[3 2 1 4];[1 2 3 4]];

% Precoder Matrix 'F4' For QAM-16
F4=zeros(4,4,16);
for i=1:length(W)
    F4(:, :, i) = W(:, F4_matrix_order(i,:), i) / 2 ;
end

F2_matrix_order = ...
[[1 4];[1 2];[1 2];[1 2];
 [1 4];[1 4];[1 3];[1 3];
 [1 2];[1 4];[1 3];[1 3];
 [1 2];[1 3];[1 3];[1 2]];

% Precoder Matrix 'F2' For QAM-16
F2=zeros(4,2,16);
for i=1:length(W)
    F2(:, :, i) = W(:, F2_matrix_order(i,:), i) / sqrt(2) ;
end

if K == 2
    precoder_idx = 1:1:16;
elseif K == 4
    precoder_idx = [1, 2, 5, 6, 13];
end

% In case of K == 4, note that we only consider 5 precoder 
% because the rest of them are identical : F0, F1, F4, F5, F12

%% Basic Setup For searching the minimum distance

min_dist_per_precoder = zeros(length(precoder_idx), 1);
curr_precoder_idx = 1;
[~, perm_size_total] = size(perms_array);

for i = precoder_idx
    R_matries = zeros( 4, K, perm_size_total);
    Q_matries = zeros( 4, 4, perm_size_total);
    
    if K == 2
        Hp = H * F2(:, :, i);
    elseif K == 4
        Hp = H * F4(:, :, i);
    end
    % Complex_LLL function returns [Q, R, P]
    [~, ~, P] = Complex_LLL(Hp, delta);
    Hp_P = Hp * P;
    
    % perms_matrix = eye(K); perms_matrix = perms_matrix(:, perms_idx);
    % Utilize that -> 'Hp * perms_matrix == Hp(:, perms_idx);'
    
    R_cas_min = realmin;  
    for v=1:perm_size_total
        % count = count + 1;
        perm_idx = perms_tree(1, v);
        neighbor_perm_idx = perms_tree(2, v);
        perms_pattern = perms_array(:, perm_idx);
        perms_pattern_parent = perms_array(:, neighbor_perm_idx);
        % Find perm matrix
        perms_matrix_child = eye(K); perms_matrix_child = perms_matrix_child(:, perms_pattern);
        perms_matrix_parent = eye(K); perms_matrix_parent = perms_matrix_parent(:, perms_pattern_parent);
        
        % First time to start permumation pattern
        if perm_idx == 1
            perm_Hp_P = Hp_P(:, perms_pattern);
            [Q_matries(:, :, perm_idx), R_matries(:, :, perm_idx)] = qr(perm_Hp_P);
            R_cas_min = max(R_cas_min, min(diag(R_matries(:, :, perm_idx)).^2));
        else
            perms_matrix = perms_matrix_parent \ perms_matrix_child;
            w = find(perms_pattern - perms_pattern_parent, 1);
            % Generate Givens Rotation Matrix
            R_rot_matries = R_matries(:, :, neighbor_perm_idx) * perms_matrix;
            alpha = R_rot_matries(w,w); beta = R_rot_matries(w+1,w);
            non_zero_r = sqrt(abs(alpha)^2 + abs(beta)^2);
            given_c = alpha / non_zero_r; given_s = beta / non_zero_r;
            G = eye(4); 
            G(w,w) = conj(given_c); G(w, w+1) = conj(given_s);
            G(w+1,w) = -given_s; G(w+1, w+1) = given_c;
            Q_matries(:, :, perm_idx) = Q_matries(:, :, neighbor_perm_idx) * G';
            R_matries(:, :, perm_idx) = G * R_rot_matries;

            R_cas_min = max(R_cas_min, min(diag(R_matries(:, :, perm_idx)).^2));
        
        end
    end   
    min_dist_per_precoder(curr_precoder_idx) = R_cas_min;
    curr_precoder_idx = curr_precoder_idx + 1;
end

%% Find the precoder having mini-max distance from the candidate error vectors
[~, max_idx] = max(min_dist_per_precoder);
if K == 2
    F = F2(:, :, precoder_idx(max_idx));
elseif K == 4
    F = F4(:, :, precoder_idx(max_idx));
end
index = max_idx;
end