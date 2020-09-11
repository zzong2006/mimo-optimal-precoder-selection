clear all;
j = sqrt(-1);
seed = 200;
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
    a = 2 * u(:, i) * u(:, i)';
    b = u(:, i)' * u(:, i);
    W(:, :, i) = eye(4) - (a) / b;
end    

F4_matrix_order = ...
[[1 2 3 4];[1 2 3 4];[3 2 1 4];[3 2 1 4];
 [1 2 3 4];[1 2 3 4];[1 3 2 4];[1 3 2 4];
 [1 2 3 4];[1 2 3 4];[1 3 2 4];[1 3 2 4];
 [1 2 3 4];[1 3 2 4];[3 2 1 4];[1 2 3 4]];

F4=zeros(4,4,16);
for i=1:length(W)
    F4(:, :, i) = W(:, F4_matrix_order(i,:), i);
end

F2_matrix_order = ...
[[1 4];[1 2];[1 2];[1 2];
 [1 4];[1 4];[1 3];[1 3];
 [1 2];[1 4];[1 3];[1 3];
 [1 2];[1 3];[1 3];[1 2]];
F2=zeros(4,2,16);
for i=1:length(W)
    F2(:, :, i) = W(:, F2_matrix_order(i,:), i);
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

possible_x = zeros(4, 16^4);
x_idx = 0;

% 가능한 모든 x vector(possible_x)를 구함
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

rng('default');
rng(seed);
H_list = ones(nR, nT, dataSize);
dist_list = ones(nR, nT, dataSize);
F_list_Accorded_H = ones(1, dataSize);

% nS = 4 or nS = 2
if nS == 4
    F = F4;
else
    F = F2;
end

matrix_step_size = 100;

tic
for idx = 1:dataSize 
    H=1/sqrt(nS)*sqrt(1/2)*(randn(nR,nT)+j*randn(nR,nT));
    H_list(:, :, idx) = H;
    [F_list_Accorded_H(idx), dist_list(:, :, idx)] = get_max_F_idx(idx, H, F, possible_x, matrix_step_size);
end
toc



function [max_F_idx, F_dist_list] = get_max_F_idx(outter_idx, H, F, possible_x, matrix_step_size)
    possible_x_device = gpuArray(possible_x);
    F_device = gpuArray(F);
    H_device = gpuArray(H);

%     possible_x_device = possible_x;
%     F_device = F;
%     H_device = H;
    
    possible_num = length(possible_x);
    F_dist_list = zeros(1, 16);
    tic
    for F_idx = 1:16
        min_val = intmax;
        for idx_x = 1:floor(possible_num / matrix_step_size) + 1
            if idx_x ~= floor(possible_num / matrix_step_size) + 1
                repeated_x = repmat(possible_x_device, 1, 1, matrix_step_size);
                [x_1, x_2, x_3]=size(repeated_x);
                repeated_x(:,1 + (idx_x - 1) * matrix_step_size:x_2+1:end)=[];
                repeated_x = reshape(repeated_x, x_1, [], x_3);
                temp_repelemed = repelem(possible_x_device(:,1 + (idx_x - 1) * matrix_step_size:(idx_x) * matrix_step_size), ...
                                    1,possible_num - 1);
                single_repeated_x = reshape(temp_repelemed,[],possible_num-1,matrix_step_size);
                subbed_repeated_x = repeated_x - single_repeated_x;
            else
                remained_size = mod(possible_num, matrix_step_size);
                repeated_x = repmat(possible_x, 1, 1, remained_size);
                [x_1, x_2, x_3]=size(repeated_x);
                repeated_x(:,1 + (idx_x - 1) * matrix_step_size:x_2+1:end)=[];
                repeated_x = reshape(repeated_x, x_1, [], x_3);
                temp_repelemed = repelem(possible_x_device(:,1 + (idx_x - 1) * matrix_step_size:end), ...
                                    1,possible_num - 1);
                single_repeated_x = reshape(temp_repelemed,[],possible_num-1,remained_size);
                subbed_repeated_x = repeated_x - single_repeated_x;
            end
            Y = H_device * F_device(: , : , F_idx);
            celled_repeated_x = cellfun(@(x) Y * x, num2cell(subbed_repeated_x, [1 2]),'UniformOutput',false);
            celled_repeated_x = cat(3, celled_repeated_x{:});
            dist = vecnorm(celled_repeated_x, 2);
            min_dist = min(dist,[],'all');
            % fprintf('%d %f\n', idx_x, min_dist);
            % temp_possible_x = possible_x;
            % picked_possible_x = possible_x(:, idx_x);
            % temp_possible_x(:, idx_x) = [];
            % dist = vecnorm(H * F(:, :, F_idx) * (temp_possible_x - picked_possible_x),2);
            % min_dist = min(dist);
            if min_val > min_dist
                min_val = min_dist;
            end
            % fprintf('%d %d %f\n',F_idx, x_1, min_dist);
        end
        F_dist_list(F_idx) = gather(min_val);
    end
    [true_max_val, max_idx] = max(F_dist_list);
    F_dist_list = reshape(F_dist_list,4,4);
    fprintf('------ %d ------\n', outter_idx);
    disp(H);
    disp(F_dist_list);
    fprintf('Max F : %d %f\n',true_max_val, max_idx);
    max_F_idx = max_idx;
    toc
end
