% filename:make_label_ML_detection.m
% 라벨링 말고 스코어 개념을 도입한 버전
nS=4;
k=4;
sample_num = 1000000;
test_sample_num = 1000000;
small_label = zeros(sample_num, 5);     % Specific precoders (1,2,5,6,13)
large_label = zeros(sample_num ,16);    % All precoders
H_list = zeros(sample_num, 32);
H_original_list = zeros(sample_num, 16);

if k==2
    delta_x_set=[-1:1:1];
elseif k==4
    delta_x_set=[-3:1:3];
elseif k==6
    delta_x_set=[-7:1:7];
end

tic;
% ticBytes(gcp);
parfor ii=1:sample_num
    H = 1/sqrt(nS)*sqrt(1/2)*(randn(4,4)+1i*randn(4,4));
    % H = reshape(arr_H(ii, :), 4, 4).';
    
    [F, small_label(ii, :)] = precoder_select_ML_ver5(H, nS, k, 1);

    % disp(reshape(large_label(ii, :),4,4).');
    H_input=vertcat(real(H),imag(H));
    training_example = reshape(transpose(H_input),1, []);
    
    H_list(ii, : ) = training_example;
    H_original_list(ii, : ) = reshape(H, 1, []);
    
    % disp(H);
end
% small_label = large_label(:, [1, 2, 5, 6, 13]);
% tocBytes(gcp)
toc

% label_filename = string(strjoin({'4x4MIMO_ML_label_', num2str(sample_num,'%d'), '.mat'}));
% save(['MIMO_ML_train_all_label_k_' num2str(k, '%d') '_' num2str(sample_num,'%d') '.mat'], 'large_label')
save(['MIMO_ML_train_some_label_k_' num2str(k, '%d') '_' num2str(sample_num,'%d') '.mat'], 'small_label')
save(['MIMO_ML_train_data_k_' num2str(k, '%d') '_' num2str(sample_num,'%d') '.mat'], 'H_list')
save(['MIMO_ML_train_original_data_k_' num2str(k, '%d') '_' num2str(sample_num,'%d') '.mat'], 'H_original_list')

% test_data

test_small_label = zeros(sample_num, 5);     % Specific precoders (1,2,5,6,13)
test_large_label = zeros(test_sample_num ,16);    % All precoders
test_H_list = zeros(test_sample_num, 32);
test_H_original_list = zeros(test_sample_num, 16);

tic;
parfor ii=1:test_sample_num
    H = 1/sqrt(nS)*sqrt(1/2)*(randn(4,4)+1i*randn(4,4));
    % H = reshape(arr_H(ii, :), 4, 4).';
    
    [F, test_small_label(ii, :)] = precoder_select_ML_ver5(H, nS, k, 1);

    H_input=vertcat(real(H),imag(H));
    training_example = reshape(transpose(H_input),1, []);
    
    test_H_list(ii, : ) = training_example;
    test_H_original_list(ii, : ) = reshape(H, 1, []);
    % disp(H);
end
% test_small_label = test_large_label(:, [1, 2, 5, 6, 13]);
toc

% label_filename = string(strjoin({'4x4MIMO_ML_label_', num2str(sample_num,'%d'), '.mat'}));
% save(['MIMO_ML_test_all_label_k_' num2str(k, '%d') '_' num2str(test_sample_num,'%d') '.mat'], 'test_large_label')
save(['MIMO_ML_test_some_label_k_' num2str(k, '%d') '_' num2str(test_sample_num,'%d') '.mat'], 'test_small_label')
save(['MIMO_ML_test_data_k_' num2str(k, '%d') '_' num2str(test_sample_num,'%d') '.mat'], 'test_H_list')
save(['MIMO_ML_test_original_data_k_' num2str(k, '%d') '_' num2str(test_sample_num,'%d') '.mat'], 'test_H_original_list')
