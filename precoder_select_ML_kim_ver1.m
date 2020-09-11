function [F]=precoder_select_ML_kim_ver1(H,nS)
%precoder_select_ML_ver9 
j= sqrt(-1);
%% Default Value
if nargin < 1
    rng('default');
    rng(1512);
    disp('[Message] Proposed_Precoder_Selection : Default set H, nS');
    nT = 4; nR = 4; 
    nS = 4;
    H = 1/sqrt(nS)*sqrt(1/2)*(randn(nR,nR)+j*randn(nR,nR));
elseif nargin < 2
    disp('[Message] Proposed_Precoder_Selection : Default set only nS');
    nS = 4;
end

%%

%LTE-A codebook
u=[ 1,  1,  1,  1,  1,              1,             1,             1,              1,  1,  1,  1,  1,  1,  1, 1;
   -1, -j,  1,  j, (-1-j)/sqrt(2), (1-j)/sqrt(2), (1+j)/sqrt(2), (-1+j)/sqrt(2), -1, -j,  1,  j, -1, -1,  1, 1;
   -1,  1, -1,  1, -j,              j,            -j,             j,              1, -1,  1, -1, -1,  1, -1, 1;
   -1,  j,  1, -j, (1-j)/sqrt(2),  (-1-j)/sqrt(2),(-1+j)/sqrt(2),(1+j)/sqrt(2),   1, -j, -1,  j,  1, -1, -1, 1];

W=zeros(4,4,16);

for ii=1:length(W)
    a = u(:, ii) * u(:, ii)';    b = u(:, ii)' * u(:, ii);    W(:, :, ii) = eye(4) - (2 * a) / b;
end    

F4_matrix_order = ...
[[1 2 3 4];[1 2 3 4];[3 2 1 4];[3 2 1 4];
 [1 2 3 4];[1 2 3 4];[1 3 2 4];[1 3 2 4];
 [1 2 3 4];[1 2 3 4];[1 3 2 4];[1 3 2 4];
 [1 2 3 4];[1 3 2 4];[3 2 1 4];[1 2 3 4]];

F4=zeros(4,4,16);
for i=1:length(W), F4(:, :, i) = W(:, F4_matrix_order(i,:), i) 	; end

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


%%

delta_x_set = zeros(49, 1); 
% plus_cond = false((2 * nS - 1)^2, 1); 
max_of_min_radius = 0;

x_list_init= complex(ones(nS,length(delta_x_set)).*(100));
x_now_init = complex(ones(nS,1).*(100));

% for rx = x_set
%     for ix = x_set
%         delta_x_set(curr_delta_x_idx) = rx + j * ix;
%         if rx >= 0
%             plus_cond(curr_delta_x_idx) = true;
%         end
%         curr_delta_x_idx = curr_delta_x_idx + 1;
%     end
% end

delta_x_set=[  -3.0000 - 3.0000i;
  -3.0000 - 2.0000i;
  -3.0000 - 1.0000i;
  -3.0000 + 0.0000i;
  -3.0000 + 1.0000i;
  -3.0000 + 2.0000i;
  -3.0000 + 3.0000i;
  -2.0000 - 3.0000i;
  -2.0000 - 2.0000i;
  -2.0000 - 1.0000i;
  -2.0000 + 0.0000i;
  -2.0000 + 1.0000i;
  -2.0000 + 2.0000i;
  -2.0000 + 3.0000i;
  -1.0000 - 3.0000i;
  -1.0000 - 2.0000i;
  -1.0000 - 1.0000i;
  -1.0000 + 0.0000i;
  -1.0000 + 1.0000i;
  -1.0000 + 2.0000i;
  -1.0000 + 3.0000i;
   0.0000 - 3.0000i;
   0.0000 - 2.0000i;
   0.0000 - 1.0000i;
   0.0000 + 0.0000i;
   0.0000 + 1.0000i;
   0.0000 + 2.0000i;
   0.0000 + 3.0000i;
   1.0000 - 3.0000i;
   1.0000 - 2.0000i;
   1.0000 - 1.0000i;
   1.0000 + 0.0000i;
   1.0000 + 1.0000i;
   1.0000 + 2.0000i;
   1.0000 + 3.0000i;
   2.0000 - 3.0000i;
   2.0000 - 2.0000i;
   2.0000 - 1.0000i;
   2.0000 + 0.0000i;
   2.0000 + 1.0000i;
   2.0000 + 2.0000i;
   2.0000 + 3.0000i;
   3.0000 - 3.0000i;
   3.0000 - 2.0000i;
   3.0000 - 1.0000i;
   3.0000 + 0.0000i;
   3.0000 + 1.0000i;
   3.0000 + 2.0000i;
   3.0000 + 3.0000i];

% plus_cond([22:49]) = true;


if nS == 4
    precoder_index=[1,2,5,6,13];
elseif nS == 2
    precoder_index = 1:1:16;
end

%%
% abs_delta_x_set = abs(delta_x_set);

%%

for precoder_idx=precoder_index
    candidate_count = zeros(nS, 1); 
    candidate_idx = ones(nS, 1);
    searchOver = false; 
    
    if nS == 4
        HF=H*F4(:,:,precoder_idx);
    else
        HF=H*F2(:,:,precoder_idx);
    end
        
    [~,RR] = qr(HF); 
    size_HF = size(HF);
    if size_HF(1) ~= size_HF(2)
        RR = RR(1:min(size_HF),:);
    end
    RR = complex(RR);
    RR = diag(sign(diag(RR))) * RR;
    min_radius = min(vecnorm(RR)); 
    
    if min_radius<=max_of_min_radius
        continue; 
    end
    
    stage_idx = 1; Is_New_Stage = true; x_now = x_now_init; x_list = x_list_init;
    
    while ~searchOver
        if stage_idx == nS
            x_1_temp=round(-RR(1,2:end)*x_now(2:end)/RR(1,1));
            x_now(1)=min(max(real(x_1_temp),-3),3)+j*min(max(imag(x_1_temp),-3),3);
            
            if ~any(x_now)
                stage_idx=stage_idx-1; 
                Is_New_Stage=false; % if metric_temp is zero vector
            else
                metric_temp=norm(RR*x_now); % calculate the radius
                % if metric_temp is less than max_of_min_radius, stop and start next procoder_idx 
                if metric_temp <= max_of_min_radius
                    min_radius=metric_temp; 
                    break; 
                end
                if metric_temp < min_radius
                    min_radius=metric_temp; 
                    stage_idx=1; 
                    Is_New_Stage=true; 
                else
                    stage_idx=stage_idx-1; 
                    Is_New_Stage=false;
                end
            end
        else
            if Is_New_Stage  
                % Refer "MIMO-OFDM Wireless Communications with MATLAB" 332p
                R_value_stage = RR(nS-stage_idx+1,nS-stage_idx+1);
                RR_mul_x_now = RR(nS-stage_idx+1:end, nS-stage_idx+2:end)*x_now(nS-stage_idx+2:end);
                if stage_idx == 1, middle = min_radius^2; coord = 0;
                else
                    middle = min_radius^2-sum(abs(RR_mul_x_now(2:end)).^2);
                    coord = RR_mul_x_now(1) / R_value_stage;                    
                end
                squared_x_set = abs(delta_x_set + coord).^2 ;
                range = middle / (R_value_stage^2);
                another_delta_x_cond = squared_x_set <= range;
              
                if stage_idx==1    
                    % another_delta_x_cond = another_delta_x_cond & plus_cond;
                    another_delta_x_cond(1:21) = false;
                end
                
                len_B = sum(another_delta_x_cond);
                x_list(nS-stage_idx+1, 1:len_B)=delta_x_set(another_delta_x_cond);
                candidate_count(stage_idx) = len_B;
                candidate_idx(stage_idx) = 1;

            end
            if candidate_count(stage_idx) == 0 
                %if the stage does not contain any candidate symbol
                if stage_idx==1, searchOver=true; 
                else
                    Is_New_Stage=false; stage_idx=stage_idx-1;
                end
            else % candidates exist
                x_now(nS-stage_idx+1)=x_list(nS-stage_idx+1, candidate_idx(stage_idx));
                Is_New_Stage=true; 
                candidate_idx(stage_idx) = candidate_idx(stage_idx) + 1;
                candidate_count(stage_idx) = candidate_count(stage_idx) - 1;
                stage_idx=stage_idx+1;
            end % if candidate_count(stage_idx) == 0 
        end
    end     % While
    if min_radius > max_of_min_radius
       max_of_min_radius = min_radius;
       opt_precoder_idx = precoder_idx;
    end
end
if nS ~= 2
    F = F4(:,:,opt_precoder_idx);
else
    F = F2(:,:,opt_precoder_idx);
end
% disp(time_temp / time_curr);