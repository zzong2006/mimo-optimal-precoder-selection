% x = randsample(5, 5, true);
% w = randsample(5, 5, true);
% a = clock;
% filename = "xwfile_"+ a(1)+ a(2)+ a(3)+ "_"  + a(4) + a(5) + ".mat";
% save(filename,'x','w');
% % 

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

nR = 4;
nT = 4;
j = sqrt(-1);
rng('default');
rng(129);
dataSize = 1000;
Ex_list_Accorded_H = zeros(1, dataSize);
Sd_list_Accorded_H = zeros(1, dataSize); 
BTC_list_Accorded_H = zeros(1, dataSize); 

% [1000] Acc (BTC) : 81.80 %
% [1000] Acc (Complex) : 99.30 %
% [1000] Acc (Real) : 99.40 %
diff_H = zeros(nR, nT, 100);
curr_idx = 1;
diff_count = 0;
for idx = 1:dataSize 
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
    [~, proposed_idx] = precoder_select_ML_ver4(H, nS);
    [~, previous_idx] = precoder_select_ML_conference(H, nS);
    fprintf('proposed : %d, previous : %d\n', proposed_idx, previous_idx);
    if proposed_idx ~= previous_idx
        disp('Diff!')
        diff_count = diff_count + 1;
    end
end
disp(diff_count)