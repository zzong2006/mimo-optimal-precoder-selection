j = sqrt(-1);
seed = 200;
comparison = false;
nS = 4;
nR = 4;
nT = 4;
dataSize = 1;

%LTE-A codebook
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

F4=zeros(4,4,16);
for i=1:length(W)
    F4(:, :, i) = W(:, F4_matrix_order(i,:), i) / 2;
end

F2_matrix_order = ...
[[1 4];[1 2];[1 2];[1 2];
 [1 4];[1 4];[1 3];[1 3];
 [1 2];[1 4];[1 3];[1 3];
 [1 2];[1 3];[1 3];[1 2]];
F2=zeros(4,2,16);
for i=1:length(W)
    F2(:, :, i) = W(:, F2_matrix_order(i,:), i) / sqrt(2);
end

x = [[0,0,0,0]; [0,0,0,1]; [0,0,1,0]; [0,0,1,1]; 
    [0,1,0,0]; [0,1,0,1]; [0,1,1,0]; [0,1,1,1]; 
    [1,0,0,0]; [1,0,0,1]; [1,0,1,0]; [1,0,1,1];
    [1,1,0,0]; [1,1,0,1]; [1,1,1,0]; [1,1,1,1]];
mapped_x = zeros(1, 16);

% QAM16 Mapping
parfor idx=1:16
    cv = QAM_mapper(x(idx, :), 4);
    mapped_x(idx) = cv;
end

possible_x = zeros(nS, 16^(nS));
x_idx = 0;

% 가능한 모든 x vector(possible_x)를 구함
if nS == 4
    for a = mapped_x
        for b = mapped_x
            for c = mapped_x
                for d = mapped_x 
                    x_idx = x_idx + 1;
                    possible_x(:, x_idx) = [a, b, c, d];
                end
            end
        end
    end
% in case of F2
else    
    for a = mapped_x
        for b = mapped_x
            x_idx = x_idx + 1;
            possible_x(:, x_idx) = [a, b];
        end
    end
end



if comparison 
    rng('default');
    rng(129);
    dataSize = 0;
    Ex_list_Accorded_H = zeros(1, dataSize);
    Sd_list_Accorded_H = zeros(1, dataSize); 
    BTC_list_Accorded_H = zeros(1, dataSize); 

    % [1000] Acc (BTC) : 81.80 %
    % [1000] Acc (Complex) : 99.30 %
    % [1000] Acc (Real) : 99.40 %
    diff_H = zeros(nR, nT, 100);
    curr_idx = 1;

    for idx = 1:dataSize 
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        Ex_list_Accorded_H(idx) = get_max_F_idx(H, F, nS, possible_x);
        % Sd_list_Accorded_H(idx) = precoder_select_ML_kim_ver1(H, nS);
        % Sd_list_Accorded_H(idx) = precoder_select_ML_ver6(H, nS);
        Sd_list_Accorded_H(idx) = precoder_select_ML_ver6_complex_ver1(H, nS);
        % BTC_list_Accorded_H(idx) = QRD_based_BT_C_Method(H, nS, 0.99);
        if Ex_list_Accorded_H(idx) ~= Sd_list_Accorded_H(idx)
            diff_H(:, :, curr_idx) = H;
            curr_idx = curr_idx + 1;
        end
        if mod(idx, 10) == 0
            fprintf('[%d] Acc : %.2f %%\n', idx, sum(Ex_list_Accorded_H == Sd_list_Accorded_H)/dataSize * 100.0);
            % fprintf('[%d] Acc : %.2f %%\n', idx, sum(Ex_list_Accorded_H == BTC_list_Accorded_H)/dataSize * 100.0);
        end
    end
end



nS=4;   %number of streams
nT=4;   %number of transmit antenna
nR=4;   %number of receive antenna

k=4;    %2,4,6

frame_size=nS*k;      % frame_size = 8;
MAX_frame_NUM=50000000;
EbNo_Start=7;%[dB]
EbNo_End=12;
EbNo_Step=1;
EbNo=EbNo_Start:EbNo_Step:EbNo_End;       %[dB]
BER_target=5 * 10^(-6);
bit_error_target=1000;
threshold =  bit_error_target / (frame_size * BER_target) ; disp(threshold);
seed = 6807;

% ---- exhasutive search
% EbNo_idx 7 is over ) 0.0157588
% EbNo_idx 8 is over ) 0.00827875
% EbNo_idx 9 is over ) 0.0040462

% ---- Proposed
% EbNo_idx 7 is over ) 0.0155822
% EbNo_idx 8 is over ) 0.00821084
% EbNo_idx 9 is over ) 0.00407932

% ---- QRD-BTE
% EbNo_idx 7 is over ) 0.0144037
% EbNo_idx 8 is over ) 0.00799265
% EbNo_idx 9 is over ) 0.00387936

% nS = 4 or nS = 2
if nS == 4
    F = F4;
else
    F = F2;
end

if nS == 4
    precoder_index=[1,2,5,6,13];
elseif nS == 2
    precoder_index = 1:1:16;
end

BER=zeros(1,length(EbNo));
i = 0;
for EbNo_idx=EbNo % needs to modify
    i = i + 1;
    bit_error=0;
    rng('default');
    rng(seed);
    %noise variance calculation
    n_var=10^(-EbNo(i)/10)/k;
    tic
    for frame_idx=1:MAX_frame_NUM
        disp(frame_idx)
        %random frame generation
        bits=randi([0 1],1,frame_size);
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column

        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR, nR)+j*randn(nR, nR));
        
        F = get_max_F_idx(H, F, nS, possible_x, precoder_index);
        % F = precoder_select_ML_kim_ver1(H, nS);
        % F = QRD_based_BT_E_Method(H, nS, 1, 0.99);
        % Use selected Random Precoder
        % n : Gaussian noise
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y =H*F*x+n;
        
        %receiver
        %Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        
        s_hat=Sphere_Decoder_Complex(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);

        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if mod(bit_error, 10) == 0 && bit_error ~= 0
            fprintf('[%d %%] EbNo_idx %d ) %g \n', bit_error/bit_error_target * 100,EbNo(EbNo_idx) ,bit_error/(frame_idx*frame_size));
        end
        if bit_error > bit_error_target
            BER(i)=bit_error/(frame_idx*frame_size);
            fprintf('EbNo_idx %d is over ) %g\n',EbNo(i), bit_error/(frame_idx*frame_size));
            break;
        end

        if bit_error_target /(frame_idx*frame_size) < BER_target
            BER(i)=bit_error/(frame_idx*frame_size);
            fprintf('EbNo_idx %d is over ) %g\n',EbNo(i), bit_error/(frame_idx*frame_size));
            break;
        end
    end
    toc
    if BER(i) < BER_target
        break; 
    end
end


function F_return = get_max_F_idx(H, F, nS, possible_x, precoder_index)
    if nS == 2 
        F_dist_list = zeros(1, 16);
        parfor F_idx = 1:16
            min_val = intmax;
            for x_1 = 1:length(possible_x)
                temp_possible_x = possible_x;
                picked_possible_x = possible_x(:, x_1);
                temp_possible_x(:, x_1) = [];
                diff_x = temp_possible_x - picked_possible_x;
                dist = vecnorm(H * F(:, :, F_idx) * (diff_x),2);
                [min_dist, min_idx] = min(dist);
                if min_val > min_dist
                    min_val = min_dist;
                    % fprintf('%d th precoder Min value changed : %f\n', F_idx, min_val);
                    % min_x_list (:, F_idx) = diff_x(:, min_idx);
                end
            end
            F_dist_list(F_idx) = min_val;
            % fprintf('%d th precoder is over. %f\n', F_idx, min_val);
        end
    else
        F_dist_list = zeros(1, length(precoder_index));
        parfor idx = 1:length(precoder_index)
            F_idx = precoder_index(idx)
            min_val = intmax;
            for x_1 = 1:length(possible_x)
                temp_possible_x = possible_x;
                picked_possible_x = possible_x(:, x_1);
                temp_possible_x(:, x_1) = [];
                diff_x = temp_possible_x - picked_possible_x;
                dist = vecnorm(H * F(:, :, F_idx) * (diff_x),2);
                [min_dist, ~] = min(dist);
                if min_val > min_dist
                    min_val = min_dist;
                    % fprintf('%d th precoder Min value changed : %f\n', F_idx, min_val);
                    % min_x_list (:, F_idx) = diff_x(:, min_idx);
                end
            end
            F_dist_list(idx) = min_val;
            % fprintf('%d th precoder is over. %f\n', F_idx, min_val);
        end
    end
    [~, max_idx] = max(F_dist_list);
    F_return = F(:, :, max_idx);
end


