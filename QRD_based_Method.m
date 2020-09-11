% Made by Woosung 2020.Jan.14th

function [F, index] = QRD_based_Method(H, K)
% < Input >
% K : The number of the bit streams
% H : Channel, the size of nR X nT
% < Output >
% F : Selected Sub-optimal Precoder
% index : And its index

%% Default Value
if nargin < 1
    disp('[Message] QRD_Precoder_Selection : Default set H, K');
    nT = 2; nR = 4; 
    K = 2;
    H = 1/sqrt(K)*sqrt(1/2)*(randn(nR,nR)+1i*randn(nR,nR));
elseif nargin < 2
    disp('[Message] QRD_Precoder_Selection : Default set only K');
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

% Precoder Matrix 'F2' For diag(R)QAM-16
F2=zeros(4,2,16);
for i=1:length(W)
    F2(:, :, i) = W(:, F2_matrix_order(i,:), i) / sqrt(2) ;
end

% In case of K == 4, note that we only consider 5 precoder 
% because the rest of them are identical : F0, F1, F4, F5, F12

%% Basic Setup For searching the minimum distance

if K == 2
    precoder_idx = [1:1:16];
elseif K == 4
    precoder_idx = [1, 2, 5, 6, 13];
end
min_dist_per_precoder = zeros(length(precoder_idx), 1);
curr_precoder_idx = 1;


for i = precoder_idx
%   curr_min = realmax;
    if K == 2
        Hp = H * F2(:, :, i);
    elseif K == 4
        Hp = H * F4(:, :, i);
    end
    
    % Perform the QR Decomposition
    [~, R] = qr(Hp);
    min_dist_per_precoder(curr_precoder_idx) = min(diag(R).^2);
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