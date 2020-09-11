%filename:precoded_MIMO_communication_ML_ver2.m

j=sqrt(-1);
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
    F2(:, :, i) = W(:, F2_matrix_order(i,:), i) /sqrt(2);
end

nS=2;   %number of streams
nT=4;   %number of transmit antenna
nR=2;   %number of receive antenna

k=4;    %2,4,6

frame_size=nS*k;                    %bits
% frame_size = 8;
MAX_frame_NUM=50000000;
EbNo_Start=1;%[dB]
EbNo_End=30;
EbNo_Step=1;
EbNo=[EbNo_Start:EbNo_Step:EbNo_End];       %[dB]
BER_target=5 * 10^(-6);
bit_error_target=1000;
threshold =  bit_error_target / (frame_size * BER_target) ; disp(threshold);
seed = randi(10000);

% Single (No) precoder (with Sphere Decoding)
tic
Random_precoder_SD_BER=zeros(1,length(EbNo));
for EbNo_idx=1:0
    bit_error=0;

    rng('default');
    rng(seed);

    for frame_idx=1:MAX_frame_NUM
        %random frame generation
        bits=randi([0 1],1,frame_size);

        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column

        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;

        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

        % Use selected Random Precoder
        % n : Gaussian noise
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;

        %receiver
        %Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);

        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if mod(bit_error, 100) == 0
            fprintf('[%d %%] EbNo_idx %d ) %g \n', bit_error/bit_target * 100,EbNo(EbNo_idx) ,bit_error/(frame_idx*frame_size));
        end
        if bit_error > bit_error_target
            Random_precoder_SD_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('EbNo_idx %d is over ) %g\n',EbNo(EbNo_idx), bit_error/(frame_idx*frame_size));
            break;
        end

        if bit_error_target /(frame_idx*frame_size) < BER_target
            Random_precoder_SD_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('EbNo_idx %d is over ) %g\n',EbNo(EbNo_idx), bit_error/(frame_idx*frame_size));
            break;
        end
    end
    if Random_precoder_SD_BER(EbNo_idx) < BER_target
        break; 
    end
end
% ML Precoder
disp('--------------- Single (No) Precoder Selection performance measurement is over ---------------');
toc

SMLPS_BER=zeros(1,length(EbNo));
for EbNo_idx=1:0
    bit_error=0;
    rng('default');
    rng(seed);
    for frame_idx=1:MAX_frame_NUM
        if mod(frame_idx/MAX_frame_NUM * 100, 10.0) == 0.0
            fprintf('[Message] SMLPS Detection, EbNo_idx : %d, %g \%, Bit Error : %d\n' ...
                , EbNo_idx, frame_idx/MAX_frame_NUM * 100, bit_error); 
        end
        %random frame generation
        bits=randn(1,frame_size);
        neg_idx=find(bits<0);
        bits=ones(1,frame_size);
        bits(neg_idx)=0;
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column
        
        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;
        
        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        
        % Select Precoder
        [F, ~] = Simplified_ML_Precoder_Selection(H, nS);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;
        
        %receiver
        %Rayleigh fading channel
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if bit_error > bit_error_target
            SMLPS_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[SMLPS_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
        
        if bit_error_target /(frame_idx*frame_size) < BER_target
            SMLPS_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[SMLPS_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
    end
    
    if SMLPS_BER(EbNo_idx) < BER_target
        break; 
    end
end
disp('--------------- Simplified ML Selection performance measurement is over ---------------');
toc



tic
ML_BER=zeros(1,length(EbNo));
for EbNo_idx=1:length(EbNo)
    bit_error=0;
    rng('default');
    rng(6807);
    % sprintf('EbNo_idx : %d \n', EbNo_idx); 
    for frame_idx=1:MAX_frame_NUM
        %random frame generation
        bits=randi([0 1],1,frame_size);
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column
        
        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;
        
        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        
        % Select Precoder
        F = precoder_select_ML_ver6_complex_ver1(H,nS);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;
        
        %receiver
        %Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if bit_error > bit_error_target
            ML_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[ML_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
        
        if bit_error_target /(frame_idx*frame_size) < BER_target
            ML_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[ML_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
    end
    
    if ML_BER(EbNo_idx) < BER_target
        break; 
    end
end

disp('--------------- ML Selection performance measurement is over ---------------');
toc


tic
QRD_BER=zeros(1,length(EbNo));
for EbNo_idx=1:length(EbNo)
    bit_error=0;
    rng('default');
    rng(seed);
    
    for frame_idx=1:MAX_frame_NUM
        %random frame generation
        bits=randn(1,frame_size);
        neg_idx=find(bits<0);
        bits=ones(1,frame_size);
        bits(neg_idx)=0;
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column
        
        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;
        
        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        
        % Select Precoder
        [F, ~] = QRD_based_Method(H,nS);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;
        
        %receiver
        %Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if bit_error > bit_error_target
            QRD_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
        
        if bit_error_target /(frame_idx*frame_size) < BER_target
            QRD_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
    end
    
    if QRD_BER(EbNo_idx) < BER_target
        break; 
    end
end
disp('--------------- QRD based Selection performance measurement is over ---------------');
toc


perms_array = perms(1:1:K).';
[perms_size_case, perm_size_total] = size(perms_array);
perms_idx_refer = zeros(1, 10^(perms_size_case-1) * (perms_size_case+1));

for i = 1:perm_size_total
    if perms_size_case == 2
        perms_idx = perms_array(1,i) * 10  +perms_array(2,i); 
    elseif perms_size_case == 3
        perms_idx = perms_array(1,i) * 100 + perms_array(2,i) * 10 +perms_array(3,i);
    else         
        % Considered only perms_size_case == 4
        perms_idx = perms_array(1,i) * 1000 + perms_array(2,i) * 100 +perms_array(3,i) * 10 + perms_array(4,i);
    end
    perms_idx_refer(perms_idx) = i;
end

% Generate Permutation Tree
perms_tree = zeros(2, perm_size_total);
perms_queue = zeros(size(perms_array));
tree_idx = 1; queue_curr_idx = 1;
perms_queue(:, queue_curr_idx) = perms_array(:, 1);
tree_checker = zeros(1, perm_size_total);

while queue_curr_idx ~= 0
    perms_pattern = perms_queue(: , queue_curr_idx);
    queue_curr_idx = queue_curr_idx - 1;
    perms_idx = perms_pattern(1) * 1000 + perms_pattern(2) * 100 +perms_pattern(3) * 10 + perms_pattern(4);
    parent_checker_idx = perms_idx_refer(perms_idx);
    if parent_checker_idx == 1
        tree_checker(1) = true;
        perms_tree(1, tree_idx) = 1;    % child
        perms_tree(2, tree_idx) = 0;    % parent
        tree_idx = tree_idx + 1;
    end
    for w=1:perms_size_case-1
        swapped_pattern = perms_pattern;
        % Swap
        temp_element = swapped_pattern(w+1); 
        swapped_pattern(w+1) = swapped_pattern(w); 
        swapped_pattern(w) = temp_element;
        swapped_perms_idx = swapped_pattern(1) * 1000 + swapped_pattern(2) * 100 +swapped_pattern(3) * 10 + swapped_pattern(4);
        child_checker_idx = perms_idx_refer(swapped_perms_idx);
        if ~tree_checker(child_checker_idx)
            tree_checker(child_checker_idx) = true;
            perms_tree(1, tree_idx) = child_checker_idx;    % child
            perms_tree(2, tree_idx) = parent_checker_idx;    % parent
            tree_idx = tree_idx + 1;
            queue_curr_idx = queue_curr_idx + 1;
            perms_queue(:, queue_curr_idx) = perms_array(:, child_checker_idx);
        end
    end
end

tic
QRD_BTP_BER=zeros(1,length(EbNo));
for EbNo_idx=1:length(EbNo)
    bit_error=0;
    rng('default');
    rng(seed);
    
    for frame_idx=1:MAX_frame_NUM

        %random frame generation
        bits=randn(1,frame_size);
        neg_idx=find(bits<0);
        bits=ones(1,frame_size);
        bits(neg_idx)=0;
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column
        
        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;
        
        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        
        % Select Precoder
        [F, ~] = QRD_based_BT_P_Method(H,nS);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;
        
        %receiver
        %Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if bit_error > bit_error_target
            QRD_BTP_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BTP_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
        
        if bit_error_target /(frame_idx*frame_size) < BER_target
            QRD_BTP_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BTP_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
    end
    
    if QRD_BTP_BER(EbNo_idx) < BER_target
        break; 
    end
end
disp('--------------- QRD_BT-P based Selection performance measurement is over ---------------');
toc


tic
QRD_BTC_BER=zeros(1,length(EbNo));
for EbNo_idx=1:length(EbNo)
    bit_error=0;
    rng('default');
    rng(seed);
    
    for frame_idx=1:MAX_frame_NUM
        %random frame generation
        bits=randn(1,frame_size);
        neg_idx=find(bits<0);
        bits=ones(1,frame_size);
        bits(neg_idx)=0;
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column
        
        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;
        
        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        
        % Select Precoder
        [F, ~] = QRD_based_BT_C_Method(H, nS, 0.99);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;
        
        %receiver
        %Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if bit_error > bit_error_target
            QRD_BTC_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BTC_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
        
        if bit_error_target /(frame_idx*frame_size) < BER_target
            QRD_BTC_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BTC_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
    end
    
    if QRD_BTC_BER(EbNo_idx) < BER_target
        break; 
    end
end
disp('--------------- QRD_BT-C based Selection performance measurement is over ---------------');
toc


tic
QRD_BTE_BER=zeros(1,length(EbNo));
for EbNo_idx=1:length(EbNo)
    bit_error=0;
    rng('default');
    rng(seed);
    
    for frame_idx=1:MAX_frame_NUM
        %random frame generation
        bits=randn(1,frame_size);
        neg_idx=find(bits<0);
        bits=ones(1,frame_size);
        bits(neg_idx)=0;
        
        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column
        
        %noise variance calculation
        n_var=10^(-EbNo(EbNo_idx)/10)/k;
        
        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
        
        % Select the proper Precoder
        [F, ~] = QRD_based_BT_E_Method(H, nS, 2, 0.99);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        y=H*F*x+n;
        
        % Receiver / Rayleigh fading channel
        %If F is selected by ML algorithm, then use sphere decoder to get
        % original stream(signal).
        HF = H*F;
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error=bit_error+sum(xor(bits,bits_hat));
        if bit_error > bit_error_target
            QRD_BTE_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BTE_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
        
        if bit_error_target /(frame_idx*frame_size) < BER_target
            QRD_BTE_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
            fprintf('[QRD_BTE_BER] EbNo_idx %d is over ) %g\n',EbNo_idx, bit_error/(frame_idx*frame_size));
            break;
        end
    end
    
    if QRD_BTE_BER(EbNo_idx) < BER_target
        break; 
    end
end
disp('--------------- QRD_BT-E based Selection performance measurement is over ---------------');
toc

disp('All over');
figure;
semilogy(EbNo,Random_precoder_SD_BER, 'rx-');  hold on;
semilogy(EbNo, ML_BER, 'bo-');         hold on;
semilogy(EbNo, SMLPS_BER, 'm*--');      hold on;
semilogy(EbNo, QRD_BER, 'k+-');        hold on;
semilogy(EbNo, QRD_BTP_BER, 'k^-');    hold on;
semilogy(EbNo, QRD_BTC_BER, 'ks-');    hold on;
semilogy(EbNo, QRD_BTE_BER, 'kd-.');    hold on;
title('Precoder Selection Processing Time');
xlabel('SNR[dB]');
ylabel('BER');
legend('No precoder','ML precoder', 'SML precoder', 'QRD precoder', ...
'QRD BT-P precoder', 'QRD BT-C precoder', 'QRD BT-E precoder');

saveas(gcf,'MIMO_4x4_QAM16_nS4');