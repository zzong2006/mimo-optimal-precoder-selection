% Made by Woosung 2020.Jan.14th

function [F, index] = QRD_based_BT_E_Method(H, K, K_R, delta)
% < Input >
% K : The number of the bit streams
% K_R : Reduced number of the bit streams (must be smaller than K)
% H : Channel, the size of nR X nT
% delta : Orthogonality of the transformed channel matrix H_{LR}
% < Output >
% F : Selected Sub-optimal Precoder
% index : And its inde

%% Default Value
if nargin < 1
    % rng('default');
    % rng(2254);
    disp('[Message] QRD_BT-E_Precoder_Selection : Default set H, K, K_R, delta');
    nT = 4; nR = 4; 
    K = 4;
    K_R = 2;
    H = 1/sqrt(K)*sqrt(1/2)*(randn(nT,nR)+1i*randn(nT,nR));
    delta = 0.99;
elseif nargin < 2
    disp('[Message] QRD_BT-E_Precoder_Selection : Default set K, K_R, delta');
    K = 4;
    K_R = 2;
    delta = 0.99;
elseif nargin < 3
    disp('[Message] QRD_BT-E_Precoder_Selection : Default set K_R, delta');
    K_R = 2;
    delta = 0.99;
elseif nargin < 4
    disp('[Message] QRD_BT-E_Precoder_Selection : Default set delta only');
    delta = 0.99;
end

assert(K_R < K, "3rd parameter K_R must be smaller than K (K_R < K)");

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

% In case of K == 4, note that we only consider 5 precoder 
% because the rest of them are identical : F0, F1, F4, F5, F12

%% Basic Setup For searching the minimum distance

if K == 2
    precoder_idx = 1:1:16;
elseif K == 4
    precoder_idx = [1, 2, 5, 6, 13];
end

min_dist_per_precoder = zeros(length(precoder_idx), 1);
curr_precoder_idx = 1;

for i = precoder_idx
    perms_array = perms(1:1:K_R).';
    reduced_perms_array = cat(1,perms_array,repmat((K_R + 1:1:K)', 1, size(perms_array,2)));
    %   curr_min = realmax;
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
    
    R_eff_cas_min = realmin;
    for perms_idx = reduced_perms_array      
        perm_Hp_P = Hp_P(:, perms_idx);
        % As permutation matrices are orthogonal matrices, 
        % the inverse matrix exists and can be written as (perms_matrix * perms_matrix' = I).
        % Perform the QR Decomposition
        
        [Q, R] = qr(perm_Hp_P);
        R_eff_cas_min = max(R_eff_cas_min, min(diag(R).^2));
    end   
    min_dist_per_precoder(curr_precoder_idx) = R_eff_cas_min;
    curr_precoder_idx = curr_precoder_idx + 1;
end

%% Find the precoder having mini-max distance from the candidate error vectors
[~, max_idx] = max(min_dist_per_precoder);
if K == 2
    F = F2(:, :, precoder_idx(max_idx));
elseif K == 4
    F = F4(:, :, precoder_idx(max_idx));
end
index = precoder_idx(max_idx);
end