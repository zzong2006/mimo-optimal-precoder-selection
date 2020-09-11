% Measure the processing time with timeit function.
% It seems to be more accurate to me (Personal Opinion)

j=sqrt(-1);

% 4x4 MIMO, 4 Streams
nS=4;   %number of streams
nT=4;   %number of transmit antenna
nR=4;   %number of receive antenna

% 4x2 MIMO, 2 Streams
% nS=2;   %number of streams
% nT=2;   %number of transmit antenna
% nR=4;   %number of receive antenna

k=4;    %2,4,6

frame_size=nS*k;                    %bits
seed = 6807;
symbol_vector_num = 1000;
EbNo_lvl = 28;  % originally 25

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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));

    % Select Precoder
    % mlf = @() precoder_select_ML_kim_ver1(H, nS);
    % mlf = @() precoder_select_ML_kim_ver2(H, nS);
    mlf = @() precoder_select_ML_kim_ver3(H, nS);
    % mlf = @() precoder_select_ML_ver6(H,nS);
    % mlf = @() precoder_select_ML_ver6_rev1(H,nS);
    % mlf = @() precoder_select_ML_ver7p2(H,nS);
    % mlf = @() precoder_select_ML_ver6_complex(H,nS);
    % mlf = @() precoder_select_ML_ver6_complex_ver1(H, nS);
    prs_time_ML(i) = timeit(mlf);
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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
    % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    smlf = @() Simplified_ML_Precoder_Selection(H, nS);
    prs_time_SML(i) = timeit(smlf);
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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
    % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    qrd_f = @() QRD_based_Method(H, nS);
    prs_time_QRD(i) = timeit(qrd_f);
end


prs_time_SVD = zeros(1,symbol_vector_num);
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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
    % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    qrd_f = @() SVD_based_Method(H, nS);
    prs_time_SVD(i) = timeit(qrd_f);
end

perms_array = perms(1:1:nS).';
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

if nS ~= 2
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
            perms_tree(2, tree_idx) = 1;    % parent
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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
    % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    qrd_btp_f = @() QRD_based_BT_P_Method(H, nS);
    prs_time_QRD_BTP(i) = timeit(qrd_btp_f);
end


prs_time_previous = zeros(1,symbol_vector_num);
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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
    % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));

    % Select Precoder
    previousf = @() precoder_select_ML_conference(H, nS);
    prs_time_previous(i) = timeit(previousf);
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
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
    % H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
    
    % Select Precoder
    if nS ~= 2
        qrd_bte_f = @() QRD_based_BT_E_Method(H, nS, 3, 0.99);
    else
        qrd_bte_f = @() QRD_based_BT_E_Method(H, nS, 1, 0.99);
    end
    prs_time_QRD_BTE(i) = timeit(qrd_bte_f);
end


disp('All over');

c = categorical({'Proposed', 'Proposed(Intermediate)',  'SVD-based', 'QRD-based','LR-based' ,'Limited Search Space'});
prices = [mean(prs_time_ML), mean(prs_time_previous), mean(prs_time_SVD), mean(prs_time_QRD), ...
    mean(prs_time_QRD_BTE), mean(prs_time_SML)];
std_prices = [std(prs_time_ML), std(prs_time_previous), std(prs_time_SVD), std(prs_time_QRD), ...
     std(prs_time_QRD_BTE), std(prs_time_SML)]/sqrt(symbol_vector_num) ;

[sorted_val, sort_idx] = sort(prices);
sorted_precoder_selection = c(sort_idx);
sorted_std = std_prices(sort_idx);

figure();

barh(sorted_val);
set(gca, 'YTickLabel', sorted_precoder_selection);
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.4f'));

xlim([0 2.5*10^(-3)]);
ylabel('Precoder Selection Technique');
xlabel('Processing time [sec]');

% hold on
% er = errorbar(sorted_precoder_selection, sorted_val, sorted_std);
% er.Color=[0 0 0];
% er.LineStyle = 'non';
% hold off
% title('Precoder Selection Processing Time');

% x_range = 1:1:1000; window_size = 1;
% semilogy(x_range,smoothdata(prs_time_SML,'gaussian',window_size), '-','LineWidth',2);    hold on;
% semilogy(x_range,smoothdata(prs_time_ML,'gaussian',window_size), '-','LineWidth',2);  hold on;
% semilogy(x_range,smoothdata(prs_time_QRD,'gaussian',window_size), '-','LineWidth',2);  hold on;
% semilogy(x_range,smoothdata(prs_time_QRD_BTP,'gaussian',window_size), '-','LineWidth',2);  hold on;
% semilogy(x_range,smoothdata(prs_time_QRD_BTC,'gaussian',window_size), '-','LineWidth',2);  hold on;
% semilogy(x_range,smoothdata(prs_time_QRD_BTE,'gaussian',window_size), '-','LineWidth',2);  hold on;

% legend('SML precoder','ML precoder','QRD precoder','QRD-BTP precoder','QRD-BTC precoder','QRD-BTE precoder');

saveas(gcf,'MIMO_test_fig_bar');
