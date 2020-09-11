% Made by Woosung 2020.Jan.14th

function [F, index] = Conventional_Precoder_Selection(H, K)
% < Input >
% K : The number of the bit streams
% H : Channel, the size of nR X nT
% < Output >
% F : Selected Sub-optimal Precoder
% index : And its index

%% Default Value
if nargin < 1
    disp('[Message] Conventional_Selection : Default set H, K');
    nT = 4; nR = 4; 
    K = 4;
    H = 1/sqrt(K)*sqrt(1/2)*(randn(nT,nR)+1i*randn(nT,nR));
elseif nargin < 2
    disp('[Message] Conventional_Precoder_Selection : Default set only K');
    K = 4;
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

% In case of K == 4, note that we only consider 5 precoder 
% because the rest of them are identical : F0, F1, F4, F5, F12

%% Basic Setup For searching the minimum distance

if K == 2
    precoder_idx = 1:1:16;
    omega_possible = [-1, 1];
    half_omega_possible = [1];
    mapped_x = zeros(1, 4);
    half_mapped_x = zeros(1, 2);
    xq = zeros(K, 4^4);
    xp = zeros(K, 4^3 * 2); 
elseif K == 4
    precoder_idx = [1, 2, 5, 6, 13];
    omega_possible = [-3, -1, 1, 3];
    half_omega_possible = [1, 3];
    mapped_x = zeros(1, 16);
    half_mapped_x = zeros(1, 8);
    xq = zeros(K, 16^4);
    xp = zeros(K, 16^3 * 8); 
end
min_dist_per_precoder = zeros(length(precoder_idx), 1);
curr_precoder_idx = 1;

curr_idx = 1;
% mapped_x
for a = omega_possible
    for b = omega_possible
        mapped_x(:, curr_idx) = a + b*j;
        curr_idx = curr_idx + 1;
    end
end

curr_idx = 1;
% half_mapped_x
for a = half_omega_possible
    for b = omega_possible
        half_mapped_x(:, curr_idx) = a + b*j;
        curr_idx = curr_idx + 1;
    end
end

% make xq
curr_idx = 1;
for a = mapped_x
    for b = mapped_x
        for c = mapped_x
            for d = mapped_x 
                xq(:, curr_idx) = [a, b, c, d];
                curr_idx = curr_idx + 1;
            end
        end
    end
end

% make xp
curr_idx = 1;
for a = half_mapped_x
    for b = mapped_x
        for c = mapped_x
            for d = mapped_x 
                xp(:, curr_idx) = [a, b, c, d];
                curr_idx = curr_idx + 1;
            end
        end
    end
end

for i = precoder_idx
    d_squared_min = realmax;
    if K == 2
        Hp = H * F2(:, :, i);
    elseif K == 4
        Hp = H * F4(:, :, i);
    end
    tic
    for xp_ins = xp
        for xq_ins = xq
            if xp_ins ~= xq_ins
                d_squared_min = min(d_squared_min, (norm(Hp * (xp_ins - xq_ins)))^2);
            end
        end
    end
    toc
    min_dist_per_precoder(curr_precoder_idx) = d_squared_min;
    curr_precoder_idx = curr_precoder_idx + 1;
end
[~, max_idx] = max(min_dist_per_precoder);
if K == 2
    F = F2(:, :, precoder_idx(max_idx));
elseif K == 4
    F = F4(:, :, precoder_idx(max_idx));
end
index = max_idx;
end