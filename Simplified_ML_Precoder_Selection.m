% Made by Woosung 2020.Jan.

function [F, index] = Simplified_ML_Precoder_Selection(H, K)
% < Input >
% K : The number of the bit streams
% H : Channel, the size of nR X nT
% < Output >
% F : Selected Sub-optimal Precoder
% index : And its index

%% Default Value
if nargin < 1
    disp('[Message] Simplified_ML_Precoder_Selection : Default set H, K');
    nT = 4; nR = 4; 
    K = 4;
    H = 1/sqrt(K)*sqrt(1/2)*(randn(nT,nR)+1i*randn(nT,nR));
elseif nargin < 2
    disp('[Message] Simplified_ML_Precoder_Selection : Default set only K');
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
    curr_min = realmax;
%% 1. Error vectors where one entry is non-zero.
    if K == 2
        G = H * F2(:, :, i);
    elseif K == 4
        G = H * F4(:, :, i);
    end
% 1-1. For the error vectors in the case of (1) with the norm of 1 (m_1)
% The search space : K
    for w = 1:K
        m1_val = G(:,w)' * G(:,w);
        curr_min = min(curr_min, m1_val);
    end

% 1-2. For the error vectors in the case of (2) with the norm of 2 (m_1')
% (It isn't mentioned in the paper, so I don't implement for fair game.)

%% 2. Error vectors where two entries are non-zero.

% 2-1. For the error vectors in the case of (1,1) with the norm of 2 (m_2)
% The search space : 2K(K-1) (including (2,2) case)
% But below implementation only covers K(K-1). Then where is the rest (2) ?
    for w = 1:K
        for z = 1:K
            if w == z 
                continue;
            end
            upGamma_1 = max(abs(real(G(:,w)' * G(:, z))), abs(imag(G(:,w)' * G(:,z))));
            m2_val = (G(:,w)' * G(:,w)) + (G(:,z)' * G(:,z)) - 2 * (upGamma_1);
            curr_min = min(curr_min, m2_val);
        end
    end

    
% 2-2. For the error vectors in the case of (1,2) & (2,1) with the norm of
% 3 (m_3)
% The search space : 4K(K-1)
    for w = 1:K
        for z = 1:K
            if w == z 
                continue;
            end
            upThetaMin = min((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upThetaMax = max((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upGamma_2 = max(abs(real(G(:,w)' * G(:, z)) + imag(G(:,w)' * G(:, z))) ...
                        , abs(real(G(:,w)' * G(:, z)) - imag(G(:,w)' * G(:, z))));
            m3_val = 2 * upThetaMin + upThetaMax - 2 * (upGamma_2);
            curr_min = min(curr_min, m3_val);
        end
    end

% From now, we assume that we are using Higher Order Modulation. ( > 4-QAM)
% 2-3. For the error vectors in the case of (2,2) with the norm of 4 
% (This case is just twice that of the distance in the (1,1) case.)

% 2-4. For the error vectors in the case of (4,1) & (1,4) with the norm of
% 5 (m_4)
    for w = 1:K
        for z = 1:K
            if w == z 
                continue;
            end
            upThetaMin = min((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upThetaMax = max((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upGamma_1 = max(abs(real(G(:,w)' * G(:, z))), abs(imag(G(:,w)' * G(:,z))));
            m4_val = 4 * upThetaMin + upThetaMax - 4 * (upGamma_1);
            curr_min = min(curr_min, m4_val);
        end
    end

% 2-5. For the error vectors in the case of (5,1) & (1,5) with the norm of
% 6 (m_5)
    for w = 1:K
        for z = 1:K
            if w == z 
                continue;
            end
            upThetaMin = min((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upThetaMax = max((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upGamma_3 = max(abs(2 * real(G(:,w)' * G(:, z)) + imag(G(:,w)' * G(:,z))), ...
                           abs(2 * real(G(:,w)' * G(:, z)) - imag(G(:,w)' * G(:,z))));
            upGamma_3 = max(upGamma_3, abs(real(G(:,w)' * G(:, z)) + 2 * imag(G(:,w)' * G(:,z))));
            upGamma_3 = max(upGamma_3, abs(real(G(:,w)' * G(:, z)) - 2 * imag(G(:,w)' * G(:,z))));
            m5_val = 5 * upThetaMin + upThetaMax - 2 * (upGamma_3);
            curr_min = min(curr_min, m5_val);
        end
    end
% 2-6. For the error vectors in the case of (5,2) & (2,5) with the norm of
% 7 (m_6) 
    for w = 1:K
        for z = 1:K
            if w == z 
                continue;
            end
            upThetaMin = min((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upThetaMax = max((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upGamma_4 = max(abs(3 * real(G(:,w)' * G(:, z)) + imag(G(:,w)' * G(:,z))), ...
                           abs(3 * real(G(:,w)' * G(:, z)) - imag(G(:,w)' * G(:,z))));
            upGamma_4 = max(upGamma_4, abs(real(G(:,w)' * G(:, z)) + 3 * imag(G(:,w)' * G(:,z))));
            upGamma_4 = max(upGamma_4, abs(real(G(:,w)' * G(:, z)) - 3 * imag(G(:,w)' * G(:,z))));
            m6_val = 5 * upThetaMin + 2 * upThetaMax - 2 * (upGamma_4);
            curr_min = min(curr_min, m6_val);
        end
    end
% 2-7. For the error vectors in the case of (5,4) & (4,5) with the norm of
% 9 (m_7) 
    for w = 1:K
        for z = 1:K
            if w == z 
                continue;
            end
            upThetaMin = min((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upThetaMax = max((G(:,w)' * G(:,w)), (G(:,z)' * G(:,z)));
            upGamma_3 = max(abs(2 * real(G(:,w)' * G(:, z)) + imag(G(:,w)' * G(:,z))), ...
                           abs(2 * real(G(:,w)' * G(:, z)) - imag(G(:,w)' * G(:,z))));
            upGamma_3 = max(upGamma_3, abs(real(G(:,w)' * G(:, z)) + 2 * imag(G(:,w)' * G(:,z))));
            upGamma_3 = max(upGamma_3, abs(real(G(:,w)' * G(:, z)) - 2 * imag(G(:,w)' * G(:,z))));
            m7_val = 5 * upThetaMin + 4 * upThetaMax - 4 * (upGamma_3);
            curr_min = min(curr_min, m7_val);
        end
    end
    
    min_dist_per_precoder(curr_precoder_idx) = curr_min;
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
