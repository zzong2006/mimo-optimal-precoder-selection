j = sqrt(-1);
seed = 200;
nS = 2;
nR = 2;
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

else

    for a = mapped_x
        for b = mapped_x
            x_idx = x_idx + 1;
            possible_x(:, x_idx) = [a, b];
        end
    end
end

% nS = 4 or nS = 2
if nS == 4
    F = F4;
else
    F = F2;
end


get_max_F_idx(diff_H(:, :, 5), F, nS, possible_x)

function max_F_idx = get_max_F_idx(H, F, nS, possible_x)
    F_dist_list = zeros(1, 16);
    min_x_list = zeros(nS, 16);
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
    [true_max_val, max_idx] = max(F_dist_list);
    % disp(F_dist_list);
    % disp(min_x_list);
    % fprintf('Max F : %d %f\n',true_max_val, max_idx);
    % disp(F(:, :, max_idx));
    max_F_idx = max_idx;
end