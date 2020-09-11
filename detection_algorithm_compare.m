%filename:
clear all; %close all;

j=sqrt(-1);
%LTE-A codebook
u=[ 1,  1,  1,  1,  1,              1,             1,             1,              1,  1,  1,  1,  1,  1,  1, 1;
    -1, -j,  1,  j, (-1-j)/sqrt(2), (1-j)/sqrt(2), (1+j)/sqrt(2), (-1+j)/sqrt(2), -1, -j,  1,  j, -1, -1,  1, 1;
    -1,  1, -1,  1, -j,              j,            -j,             j,              1, -1,  1, -1, -1,  1, -1, 1;
    -1,  j,  1, -j, (1-j)/sqrt(2),  (-1-j)/sqrt(2),(-1+j)/sqrt(2),(1+j)/sqrt(2),   1, -j, -1,  j,  1, -1, -1, 1];

W=zeros(4,4,16);
W(:,:,1)=eye(4)-(2*u(:,1)*(u(:,1)'))/((u(:,1)')*u(:,1));
W(:,:,2)=eye(4)-(2*u(:,2)*(u(:,2)'))/((u(:,2)')*u(:,2));
W(:,:,3)=eye(4)-(2*u(:,3)*(u(:,3)'))/((u(:,3)')*u(:,3));
W(:,:,4)=eye(4)-(2*u(:,4)*(u(:,4)'))/((u(:,4)')*u(:,4));
W(:,:,5)=eye(4)-(2*u(:,5)*(u(:,5)'))/((u(:,5)')*u(:,5));
W(:,:,6)=eye(4)-(2*u(:,6)*(u(:,6)'))/((u(:,6)')*u(:,6));
W(:,:,7)=eye(4)-(2*u(:,7)*(u(:,7)'))/((u(:,7)')*u(:,7));
W(:,:,8)=eye(4)-(2*u(:,8)*(u(:,8)'))/((u(:,8)')*u(:,8));
W(:,:,9)=eye(4)-(2*u(:,9)*(u(:,9)'))/((u(:,9)')*u(:,9));
W(:,:,10)=eye(4)-(2*u(:,10)*(u(:,10)'))/((u(:,10)')*u(:,10));
W(:,:,11)=eye(4)-(2*u(:,11)*(u(:,11)'))/((u(:,11)')*u(:,11));
W(:,:,12)=eye(4)-(2*u(:,12)*(u(:,12)'))/((u(:,12)')*u(:,12));
W(:,:,13)=eye(4)-(2*u(:,13)*(u(:,13)'))/((u(:,13)')*u(:,13));
W(:,:,14)=eye(4)-(2*u(:,14)*(u(:,14)'))/((u(:,14)')*u(:,14));
W(:,:,15)=eye(4)-(2*u(:,15)*(u(:,15)'))/((u(:,15)')*u(:,15));
W(:,:,16)=eye(4)-(2*u(:,16)*(u(:,16)'))/((u(:,16)')*u(:,16));

F4=zeros(4,4,16);
F4(:,:,1)=W(:,[1 2 3 4],1);
F4(:,:,2)=W(:,[1 2 3 4],2);
F4(:,:,3)=W(:,[3 2 1 4],3);
F4(:,:,4)=W(:,[3 2 1 4],4);
F4(:,:,5)=W(:,[1 2 3 4],5);
F4(:,:,6)=W(:,[1 2 3 4],6);
F4(:,:,7)=W(:,[1 3 2 4],7);
F4(:,:,8)=W(:,[1 3 2 4],8);
F4(:,:,9)=W(:,[1 2 3 4],9);
F4(:,:,10)=W(:,[1 2 3 4],10);
F4(:,:,11)=W(:,[1 3 2 4],11);
F4(:,:,12)=W(:,[1 3 2 4],12);
F4(:,:,13)=W(:,[1 2 3 4],13);
F4(:,:,14)=W(:,[1 3 2 4],14);
F4(:,:,15)=W(:,[3 2 1 4],15);
F4(:,:,16)=W(:,[1 2 3 4],16);

F2=zeros(4,2,16);
F2(:,:,1)=W(:,[1 4],1);
F2(:,:,2)=W(:,[1 2],2);
F2(:,:,3)=W(:,[1 2],3);
F2(:,:,4)=W(:,[1 2],4);
F2(:,:,5)=W(:,[1 4],5);
F2(:,:,6)=W(:,[1 4],6);
F2(:,:,7)=W(:,[1 3],7);
F2(:,:,8)=W(:,[1 3],8);
F2(:,:,9)=W(:,[1 2],9);
F2(:,:,10)=W(:,[1 4],10);
F2(:,:,11)=W(:,[1 3],11);
F2(:,:,12)=W(:,[1 3],12);
F2(:,:,13)=W(:,[1 2],13);
F2(:,:,14)=W(:,[1 3],14);
F2(:,:,15)=W(:,[1 3],15);
F2(:,:,16)=W(:,[1 2],16);


nS=2;%number of streams
nT=4;%number of transmit antenna
nR=4;%number of receive antenna

k=4;%2,4,6

frame_size=nS*k;%bits
MAX_frame_NUM=100000;
EbNo_Start=0;%[dB]
EbNo_End=15;
EbNo_Step=1;
EbNo=[EbNo_Start:EbNo_Step:EbNo_End];       %[dB]
BER_target=10^(-6);
bit_error_target=0;

seed = randi(10000);

% ZF Precoder 

rng('default');
rng(seed);
ZF_BER=zeros(1,length(EbNo));
H_list = zeros(length(EbNo), MAX_frame_NUM, nR, nT);

tic
parfor EbNo_idx=1:length(EbNo)
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
        ZF_index = precoder_select_ver4(H, nS, 1);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        error_list = do_all_precoders_ZF(x, H, F2, n,k, bits);
        
        %BER calculation
        bit_error= bit_error + abs(min(error_list) - error_list(ZF_index));
%         if bit_error > bit_error_target
%             disp(frame_idx)
%             ZF_BER(EbNo_idx)=bit_error/(frame_idx*frame_size);
%             break;
%         end
    end
    ZF_BER(EbNo_idx) = bit_error/(MAX_frame_NUM*frame_size)
end
disp('ZF Precoder Finished');
toc




rng('default');
rng(seed);
ML_BER=zeros(1,length(EbNo));

tic
parfor EbNo_idx=1:length(EbNo)
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
        ML_index = precoder_select_ML_ver4(H, k, nS, 1);
        n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
        error_list = do_all_precoders_ML(x, H, F2, n,k, nS, n_var, bits);
        
        %BER calculation
        bit_error= bit_error + abs(min(error_list) - error_list(ML_index));
    end
    ML_BER(EbNo_idx) = bit_error/(MAX_frame_NUM*frame_size);
end
disp('ML Precoder Finished');
toc


semilogy(EbNo,ZF_BER ,EbNo, ML_BER);
title('MIMO ZF and ML Comparison');
xlabel('Eb/No[dB]');
ylabel('Best BER - Selected BER');
legend('Zero Forcing','Maximum Likelihood');

function [error_list]=do_all_precoders_ZF(x, H, F, n, k, bits)
    bit_error = zeros(1, 16);
    
    parfor idx= 1:16
        y=H*F(:, :, idx)*x+n;        
        HF = H*F(:, :, idx);
        pinvHF=(HF'*HF)\HF';
        s_hat=QAM_slicer(transpose(pinvHF*y),k);
        bits_hat=QAM_demapper(s_hat,k);
        bit_error(idx) = sum(xor(bits, bits_hat));
    end
    
    error_list = bit_error;
end

function [error_list] = do_all_precoders_ML(x, H, F, n, k, nS, n_var, bits)
    bit_error = zeros(1, 16);
    
    parfor idx= 1:16
        y=H*F(:, :, idx)*x+n;
        HF = H*F(:, :, idx);
        s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
        bits_hat=QAM_demapper(s_hat,k);
        
        %BER calculation
        bit_error(idx) = sum(xor(bits, bits_hat));
    end
    
    error_list = bit_error;
end