% Optimal 이 Proposed 랑 다른게 있는데 실제로 다른지 확인해보는 실험
% 아무래도 이러한 버그는 distance의 중복에 의해 발생하는것 같음
% nT = nR = nS = 4에서 진행하도록 설계
clearvars;
j = sqrt(-1);

comparison = false;

nS=4;   %number of streams
nT=4;   %number of transmit antenna
nR=4;   %number of receive antenna

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

time_a = clock;
save_filename = "selectionFile_"+ time_a(1)+ time_a(2)+ time_a(3) + ...
            "_"  + time_a(4) + time_a(5) + ".mat";
load_filename = "selectionFile_202082_2328.mat";


k=4;    %2,4,6
frame_size=nS*k;      % frame_size = 8;
MAX_frame_NUM=1000;
EbNo_Start=28;%[dB]
EbNo_End=28;
EbNo_Step=1;
EbNo=EbNo_Start:EbNo_Step:EbNo_End;       %[dB]
BER_target=5 * 10^(-6);
seed = 6807;

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

load(load_filename);
% optimal_selection=zeros(1,MAX_frame_NUM);
% proposed_selection=zeros(1,MAX_frame_NUM);
% LR_selection=zeros(1,MAX_frame_NUM);

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
        fprintf('[%d] \n', frame_idx);
        
        %random frame generation
        bits=randi([0 1],1,frame_size);

        %QAM mapping
        x=QAM_mapper(bits,k);
        x=transpose(x);%row->column

        %wireless transmission
        %Rayleigh fading channel
        H=1/sqrt(nS)*sqrt(1/2)*(randn(nR, nR)+j*randn(nR, nR));
        if optimal_selection(frame_idx) ~= proposed_selection(frame_idx)
            
            fprintf('Conflict) Optimal : %d, Proposed : %d\n', ...
                    optimal_selection(frame_idx), proposed_selection(frame_idx));
            

            % parfor 로 수행시, precoder 당, 약 500초(8.3분) 정도 수행
            % GPU 로 수행시, precoder 당, 약 128초(2.13분) 정도 수행
            % [opt_F, opt_dist_list] = optimal_precoder_select_with_GPUs(H, F4, nS, possible_x, precoder_index);
            [proposed_F, proposed_idx] = precoder_select_ML_kim_ver3(H, nS);
            new_proposed_idx = find(precoder_index==proposed_idx);
            for w = 1: length(opt_dist_list)
                fprintf('dist(%d) : %.30f\n', w, opt_dist_list(w))
            end

            if optimal_selection(frame_idx) == new_proposed_idx || ...
               opt_dist_list(optimal_selection(frame_idx)) == new_proposed_idx || ...
               opt_dist_list(optimal_selection(frame_idx)) == opt_dist_list(new_proposed_idx)
               if opt_dist_list(optimal_selection(frame_idx)) == opt_dist_list(new_proposed_idx)
                   new_proposed_idx = optimal_selection(frame_idx);
               end
                fprintf('Solved) Optimal : %d, Proposed : %d\n', ...
                    optimal_selection(frame_idx), new_proposed_idx);
                proposed_selection(frame_idx) = new_proposed_idx;
                save(save_filename,'optimal_selection','proposed_selection','LR_selection');
            end
            
            if opt_dist_list(optimal_selection(frame_idx)) == new_proposed_idx
                
            end
        else
            % fprintf('Skip \n');
            continue;
        end
        
        

    end
    toc
    
    % optimal_selection = proposed_selection;
end

plotting_selections(optimal_selection(1:100), proposed_selection(1:100), LR_selection(1:100));
plotting_accuracy(optimal_selection, proposed_selection, LR_selection);

function plotting_accuracy(optimal_selection, proposed_selection, LR_selection)

    proposed_accuracy = sum(optimal_selection == proposed_selection) / length(optimal_selection) * 100;
    LR_accuracy = sum(optimal_selection == LR_selection) / length(optimal_selection) * 100;
    
    figure;
    X = categorical({'Proposed', 'LR-based'});
    bar(X,[proposed_accuracy LR_accuracy],'BarWidth', 0.4);
    ylabel('Accuracy');
    xlabel('Technique');

end

function plotting_selections(optimal, proposed, LR)
    [~, trials] = size(optimal);
    y = 1:trials;
    
    figure;
    plot(y,optimal,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    hold on;
    grid on;
    grid minor;
    plot(y,proposed,'ko', 'MarkerSize', 8);
    hold on;

    ylabel('Selected precoder index');
    xlabel('Frame index');

    xlim([0 trials + 1])
    ylim([0 6])

    xbounds = ylim;
    set(gca,'YTick',xbounds(1):xbounds(2));

    hold off;
    legend('Optimal','Proposed');

    figure;

    plot(y,optimal,'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    hold on;
    grid on;
    grid minor;
    plot(y,LR,'ko', 'MarkerSize', 8);
    hold on;

    ylabel('Selected precoder index');
    xlabel('Frame index');

    xlim([0 trials + 1])
    ylim([0 6])

    xbounds = ylim;
    set(gca,'YTick',xbounds(1):xbounds(2));

    hold off;

    legend('Optimal','LR-based');
end


function [F_return, F_list] = optimal_precoder_select_with_GPUs(H, F, nS, possible_x, precoder_index)
    matrix_step_size = 150;
    if nS == 4 
    %     possible_x_device = possible_x;
    %     F_device = F;
    %     H_device = H;

        possible_num = length(possible_x);
        F_dist_list = zeros(1, length(precoder_index));                    
        
        parfor idx = 1:length(precoder_index)
            possible_x_device = gpuArray(possible_x);
            F_device = gpuArray(F);
            H_device = gpuArray(H);
            
            tic
            F_idx = precoder_index(idx);
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
                
                if min_val > min_dist
                    min_val = min_dist;
                end
                
                % fprintf('%d %d %f\n',F_idx, x_1, min_dist);
            F_dist_list(idx) = gather(min_val);
            end
            toc
        end
        disp(F_dist_list);
        [max_val, ~] = max(F_dist_list);
        max_idx = find(F_dist_list==max_val, 1, 'last');
        F_return = F(:, :, max_idx);
        F_list = F_dist_list;
        
    end
end
