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

%% HyperParameter setting

% 4x4 MIMO, 4 Streams
nS=4;   %number of streams
nT=4;   %number of transmit antenna
nR=4;   %number of receive antenna

% 4x2 MIMO, 2 Streams
nS=2;   %number of streams
nT=2;   %number of transmit antenna
nR=4;   %number of receive antenna

k=4;    %2,4,6

frame_size=nS*k;                    %bits
if nS == 4
    MAX_frame_NUM=79000000;
else 
    MAX_frame_NUM=165250000;
end

EbNo_Start=1;%[dB]
EbNo_End=20;
EbNo_Step=1;
EbNo=EbNo_Start:EbNo_Step:EbNo_End;       %[dB]
BER_target=8 * 10^(-7);
bit_error_target=1000;
threshold =  bit_error_target / (frame_size * BER_target) ; 
seed = 6807;   
file_name = 'MIMO_4x2_QAM16_nS2';
disp(threshold);
disp(seed);

NO_BER = zeros(1,length(EbNo));
SML_BER = zeros(1,length(EbNo));
ML_BER = zeros(1,length(EbNo));
QRD_BER = zeros(1,length(EbNo));
QRD_BTP_BER =zeros(1,length(EbNo));
QRD_BTC_BER =zeros(1,length(EbNo));
QRD_BTE_BER =zeros(1,length(EbNo));
SVD_BER = zeros(1, length(EbNo));

%% Read Figure which is already exists (SEED: 6807)
fig_name = 'MIMO_4x2_QAM16_nS2';
try
    fig = openfig(fig_name);
    h = findobj(gca,'Type','line');
    y_data_cell = get(h,'Ydata');
    y_data_mat = cell2mat(y_data_cell);
    Loaded_NO_BER = y_data_mat(8, :);
    Loaded_ML_BER = y_data_mat(7, :);
    Loaded_SML_BER = y_data_mat(6, :);
    Loaded_QRD_BER = y_data_mat(5, :);
    Loaded_QRD_BTP_BER = y_data_mat(4, :);
    Loaded_QRD_BTC_BER = y_data_mat(3, :);
    Loaded_QRD_BTE_BER = y_data_mat(2, :);
    Loaded_SVD_BER =y_data_mat(1, :);
catch
    disp('Cannot find the figure : set all BER value to zero (list)');
    Loaded_NO_BER = zeros(1,length(EbNo));
    Loaded_SML_BER = zeros(1,length(EbNo));
    Loaded_ML_BER = zeros(1,length(EbNo));
    Loaded_QRD_BER = zeros(1,length(EbNo));
    Loaded_QRD_BTP_BER =zeros(1,length(EbNo));
    Loaded_QRD_BTC_BER =zeros(1,length(EbNo));
    Loaded_QRD_BTE_BER =zeros(1,length(EbNo));
    Loaded_SVD_BER =zeros(1,length(EbNo));
end

precoder_selection_name = ["No_Precoder", "SPS", "Proposed", "QRD", "QRD_BT_P", "QRD_BT_C", "QRD_BT_E", "SVD"];
BER_list = zeros(length(precoder_selection_name),length(EbNo));

%%
parfor function_idx = 0:0
    tic
    BER_Bucket = zeros(1, length(EbNo));
    if function_idx == 1     % No Precoder
        BER_Bucket(1:length(Loaded_NO_BER)) = Loaded_NO_BER;
    elseif function_idx == 2 % SML
        BER_Bucket(1:length(Loaded_SML_BER)) = Loaded_SML_BER;
    elseif function_idx == 3 % ML
        BER_Bucket(1:length(Loaded_ML_BER)) = Loaded_ML_BER;
    elseif function_idx == 4 % QRD
        BER_Bucket(1:length(Loaded_QRD_BER)) = Loaded_QRD_BER;
    elseif function_idx == 5 % QRD_BT_P
        BER_Bucket(1:length(Loaded_QRD_BTP_BER)) = Loaded_QRD_BTP_BER;
    elseif function_idx == 6 % QRD_BT_C
        BER_Bucket(1:length(Loaded_QRD_BTC_BER)) = Loaded_QRD_BTC_BER;
    elseif function_idx == 7 % QRD_BT_E
        BER_Bucket(1:length(Loaded_QRD_BTE_BER)) = Loaded_QRD_BTE_BER;
    elseif function_idx == 8 % SVD
        BER_Bucket(1:length(Loaded_SVD_BER)) = Loaded_SVD_BER;
    end
    
    for EbNo_idx=1:length(EbNo)
        bit_error=0;
        rng('default');
        rng(seed);
        
        % For saving execute time, if BER value already exists, then skip.
        
        if BER_Bucket(EbNo_idx) ~= 0 
           fprintf(strcat(precoder_selection_name(function_idx), ' : EbNo_idx %d) is already exist -> %g\n'), EbNo_idx, BER_Bucket(EbNo_idx));
           continue;
        end
        for frame_idx=1:MAX_frame_NUM
            %random frame generation
            bits=randi([0 1],1,frame_size);
%             bits=randn(1,frame_size);
%             neg_idx=find(bits<0);
%             bits=ones(1,frame_size);
%             bits(neg_idx)=0;

            %QAM mapping
            x=QAM_mapper(bits,k);
            x=transpose(x);%row->column

            %noise variance calculation
            n_var=10^(-EbNo(EbNo_idx)/10)/k;

            %wireless transmission
            %Rayleigh fading channel
            % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
            H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
            
            if function_idx == 1        % No Precoder
                if nS == 2
                    F = F2(:, :, 1);             
                else                % nS == 4
                    F = F4(:, :, 1);
                end
            elseif function_idx == 2 % SML
                [F, ~] = Simplified_ML_Precoder_Selection(H, nS);
            elseif function_idx == 3    % ML
                 F = precoder_select_ML_kim_ver1(H, nS);
            elseif function_idx == 4  % QRD
                [F, ~] = QRD_based_Method(H, nS);
            elseif function_idx == 5   % QRD_BT_P
                [F, ~] = QRD_based_BT_P_Method(H, nS);
            elseif function_idx == 6 % QRD_BT_C
                [F, ~] = QRD_based_BT_C_Method(H, nS, 0.99);
            elseif function_idx == 7 % QRD_BT_E
                if nS ~= 2
                    [F, ~] = QRD_based_BT_E_Method(H, nS, 2, 0.99);
                else
                    [F, ~] = QRD_based_BT_E_Method(H, nS, 1, 0.99);
                end
            elseif function_idx == 8 % SVD
                [F, ~] = SVD_based_Method(H, nS);
            end
            % Use selected Random Precoder
            % n : Gaussian noise
            n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
            y=H*F*x+n;

            %receiver
            %Rayleigh fading channel
            %If F is selected by ML algorithm, then use sphere decoder to get
            % original stream(signal).

            HF = H*F;
            s_hat=Sphere_Decoder_Complex(y,HF,nS,k,n_var);
            bits_hat=QAM_demapper(s_hat,k);

            %BER calculation
            bit_error=bit_error+sum(xor(bits,bits_hat));
            if bit_error > bit_error_target
                BER_Bucket(EbNo_idx)=bit_error/(frame_idx*frame_size);
                fprintf(strcat(precoder_selection_name(function_idx), ' : EbNo_idx %d) is over -> %g\n'),EbNo_idx, bit_error/(frame_idx*frame_size));
                break;
            end

            if bit_error_target /(frame_idx*frame_size) < BER_target
                BER_Bucket(EbNo_idx)=bit_error/(frame_idx*frame_size);
                fprintf(strcat(precoder_selection_name(function_idx), ' : EbNo_idx %d) is over -> %g\n'),EbNo_idx, bit_error/(frame_idx*frame_size));
                break;
            end
        end
        if BER_Bucket(EbNo_idx) < BER_target
            break; 
        end
    end
    
    disp(strcat('---------------', precoder_selection_name(function_idx), ...
    ' Precoder Selection performance measurement is over ---------------'));
    if function_idx == 1        % No Precoder
       NO_BER = NO_BER + BER_Bucket;
    elseif function_idx == 2 % SML
       SML_BER =  SML_BER + BER_Bucket;
    elseif function_idx == 3    % ML
        ML_BER = ML_BER + BER_Bucket;
    elseif function_idx == 4  % QRD
        QRD_BER = QRD_BER + BER_Bucket;
    elseif function_idx == 5   % QRD_BT_P
        QRD_BTP_BER = QRD_BTP_BER + BER_Bucket;
    elseif function_idx == 6 % QRD_BT_C
        QRD_BTC_BER = QRD_BTC_BER + BER_Bucket;
    elseif function_idx == 7 % QRD_BT_E
        QRD_BTE_BER = QRD_BTE_BER + BER_Bucket;
    elseif function_idx == 8 % QRD_BT_E
        SVD_BER = SVD_BER + BER_Bucket;
    end
    toc
end


disp('All over');
figure;
semilogy(EbNo,NO_BER, 'rx-');  hold on;
semilogy(EbNo, ML_BER, 'bo-');         hold on;
semilogy(EbNo, SML_BER, 'm*--');      hold on;
semilogy(EbNo, QRD_BER, 'k+-');        hold on;
semilogy(EbNo, QRD_BTP_BER, 'k^-');    hold on;
semilogy(EbNo, QRD_BTC_BER, 'ks-');    hold on;
semilogy(EbNo, QRD_BTE_BER, 'kd-.');    hold on;
semilogy(EbNo, SVD_BER, 'k>--');    hold on;
title('BER Performance comparison for precoder selection');
xlabel('SNR[dB]');
ylabel('BER');
lgd = legend('No precoder','ML precoder', 'SML precoder', 'QRD precoder', ...
'QRD BT-P precoder', 'QRD BT-C precoder', 'QRD BT-E precoder', 'SVD precoder','Location','southwest');
lgd.NumColumns = 2;
saveas(gcf,file_name);


figure;
semilogy(EbNo,NO_BER, 'rx-');  hold on;
semilogy(EbNo, ML_BER, 'bo-');         hold on;
semilogy(EbNo, SML_BER, 'm*--');      hold on;
semilogy(EbNo, QRD_BER, 'k+-');        hold on;
% semilogy(EbNo, QRD_BTP_BER, 'k^-');    hold on;
% semilogy(EbNo, QRD_BTC_BER, 'ks-');    hold on;
semilogy(EbNo, SVD_BER, 'k>--');    hold on;
semilogy(EbNo, QRD_BTE_BER, 'kd-.');    hold on;

xlabel('SNR[dB]');
ylabel('BER');
lgd = legend('No precoder','ML precoder', 'Simplified', 'QRD'...
,'SVD','LR','Location','southwest');
lgd.NumColumns = 2;
saveas(gcf,file_name);