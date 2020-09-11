j=sqrt(-1);

nS=4;   %number of streams
nT=4;   %number of transmit antenna
nR=4;   %number of receive antenna

k=4;    %2,4,6

frame_size=nS*k;                    %bits
seed = randi(10000);
symbol_vector_num = 1000;
EbNo_lvl = 25;

prs_time_ML = zeros(1,symbol_vector_num);
rng('default');
rng(seed);
for i=1:symbol_vector_num
    bits=randn(1,frame_size);
    neg_idx=find(bits<0);
    bits=ones(1,frame_size);
    bits(neg_idx)=0;

    %QAM mapping
    x=QAM_mapper(bits,k);
    x=transpose(x);%row->column

    %noise variance calculation
    n_var=10^(-EbNo_lvl/10)/k;

    %wireless transmission
    %Rayleigh fading channel
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    tStart = tic;
    F = precoder_select_ML_ver5(H,k,nS,0);
    n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
    y=H*F*x+n;
    
    %receiver
    %Rayleigh fading channel
    %If F is selected by ML algorithm, then use sphere decoder to get
    % original stream(signal).
    HF = H*F;
    s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
    bits_hat=QAM_demapper(s_hat,k);
    prs_time_ML(i) = toc(tStart);
end


prs_time_SML = zeros(1,symbol_vector_num);
rng('default');
rng(seed);
for i=1:symbol_vector_num
    bits=randn(1,frame_size);
    neg_idx=find(bits<0);
    bits=ones(1,frame_size);
    bits(neg_idx)=0;

    %QAM mapping
    x=QAM_mapper(bits,k);
    x=transpose(x);%row->column

    %noise variance calculation
    n_var=10^(-EbNo_lvl/10)/k;

    %wireless transmission
    %Rayleigh fading channel
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    tStart = tic;
    [F, ~] = Simplified_ML_Precoder_Selection(H, nS);
    n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
    y=H*F*x+n;
    
    %receiver
    %Rayleigh fading channel
    %If F is selected by ML algorithm, then use sphere decoder to get
    % original stream(signal).
    HF = H*F;
    s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
    bits_hat=QAM_demapper(s_hat,k);
    prs_time_SML(i) = toc(tStart);
end


prs_time_QRD = zeros(1,symbol_vector_num);
rng('default');
rng(seed);
for i=1:symbol_vector_num
    bits=randn(1,frame_size);
    neg_idx=find(bits<0);
    bits=ones(1,frame_size);
    bits(neg_idx)=0;

    %QAM mapping
    x=QAM_mapper(bits,k);
    x=transpose(x);%row->column

    %noise variance calculation
    n_var=10^(-EbNo_lvl/10)/k;

    %wireless transmission
    %Rayleigh fading channel
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    tStart = tic;
    [F, ~] = QRD_based_Method(H, nS);
    n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
    y=H*F*x+n;
    
    %receiver
    %Rayleigh fading channel
    %If F is selected by ML algorithm, then use sphere decoder to get
    % original stream(signal).
    HF = H*F;
    s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
    bits_hat=QAM_demapper(s_hat,k);
    prs_time_QRD(i) = toc(tStart);
end


prs_time_QRD_BTP = zeros(1,symbol_vector_num);
rng('default');
rng(seed);
for i=1:symbol_vector_num
    bits=randn(1,frame_size);
    neg_idx=find(bits<0);
    bits=ones(1,frame_size);
    bits(neg_idx)=0;

    %QAM mapping
    x=QAM_mapper(bits,k);
    x=transpose(x);%row->column

    %noise variance calculation
    n_var=10^(-EbNo_lvl/10)/k;

    %wireless transmission
    %Rayleigh fading channel
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    tStart = tic;
    [F, ~] = QRD_based_BT_P_Method(H, nS);
    n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
    y=H*F*x+n;
    
    %receiver
    %Rayleigh fading channel
    %If F is selected by ML algorithm, then use sphere decoder to get
    % original stream(signal).
    HF = H*F;
    s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
    bits_hat=QAM_demapper(s_hat,k);
    prs_time_QRD_BTP(i) = toc(tStart);
end

prs_time_QRD_BTC = zeros(1,symbol_vector_num);
rng('default');
rng(seed);
for i=1:symbol_vector_num
    bits=randn(1,frame_size);
    neg_idx=find(bits<0);
    bits=ones(1,frame_size);
    bits(neg_idx)=0;

    %QAM mapping
    x=QAM_mapper(bits,k);
    x=transpose(x);%row->column

    %noise variance calculation
    n_var=10^(-EbNo_lvl/10)/k;

    %wireless transmission
    %Rayleigh fading channel
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    tStart = tic;
    [F, ~] = QRD_based_BT_C_Method(H, nS,0.99);
    toc(tStart)
    n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
    y=H*F*x+n;
    
    %receiver
    %Rayleigh fading channel
    %If F is selected by ML algorithm, then use sphere decoder to get
    % original stream(signal).
    HF = H*F;
    s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
    bits_hat=QAM_demapper(s_hat,k);
    prs_time_QRD_BTC(i) = toc(tStart);
end

prs_time_QRD_BTE = zeros(1,symbol_vector_num);
rng('default');
rng(seed);
for i=1:symbol_vector_num
    bits=randn(1,frame_size);
    neg_idx=find(bits<0);
    bits=ones(1,frame_size);
    bits(neg_idx)=0;

    %QAM mapping
    x=QAM_mapper(bits,k);
    x=transpose(x);%row->column

    %noise variance calculation
    n_var=10^(-EbNo_lvl/10)/k;

    %wireless transmission
    %Rayleigh fading channel
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    tStart = tic;
    [F, ~] = QRD_based_BT_E_Method(H, nS, 3, 0.99);
    n=sqrt(n_var/2)*(randn(nR,1)+j*randn(nR,1));
    y=H*F*x+n;
    
    %receiver
    %Rayleigh fading channel
    %If F is selected by ML algorithm, then use sphere decoder to get
    % original stream(signal).
    HF = H*F;
    s_hat=Sphere_Decoder(y,HF,nS,k,n_var);
    bits_hat=QAM_demapper(s_hat,k);
    prs_time_QRD_BTE(i) = toc(tStart);
end


disp('All over');

figure;
x_range = 1:1:1000; window_size = 10;
semilogy(x_range,smoothdata(prs_time_SML,'gaussian',window_size), '-','LineWidth',2);    hold on;
semilogy(x_range,smoothdata(prs_time_ML,'gaussian',window_size), '-','LineWidth',2);  hold on;
semilogy(x_range,smoothdata(prs_time_QRD,'gaussian',window_size), '-','LineWidth',2);  hold on;
semilogy(x_range,smoothdata(prs_time_QRD_BTP,'gaussian',window_size), '-','LineWidth',2);  hold on;
semilogy(x_range,smoothdata(prs_time_QRD_BTC,'gaussian',window_size), '-','LineWidth',2);  hold on;
semilogy(x_range,smoothdata(prs_time_QRD_BTE,'gaussian',window_size), '-','LineWidth',2);  hold on;

% semilogy(1:1:1000,prs_time_SML, '-','LineWidth',2);    hold on;
% semilogy(1:1:1000,prs_time_ML, '-','LineWidth',2);  hold on;
% semilogy(1:1:1000,prs_time_QRD, '-','LineWidth',2);  hold on;

title('Precoder Selection Processing Time');
xlabel('Symbol vector index');
ylabel('Processing time [sec]');
legend('SML precoder','ML precoder','QRD precoder','QRD_BTP precoder','QRD_BTC precoder','QRD_BTE precoder');

