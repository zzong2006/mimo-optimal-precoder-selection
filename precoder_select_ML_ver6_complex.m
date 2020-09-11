function [F]=precoder_select_ML_ver6_complex(H,nS)
%PRECODER_SELECT_ML_VER6_COMPLEX 
%   자세한 설명 위치

%% Default Value
if nargin < 1
    rng('default');
    rng(2254);
    disp('[Message] Proposed_Precoder_Selection : Default set H, nS');
    nT = 4; nR = 4; 
    nS = 4;
    H = 1/sqrt(nS)*sqrt(1/2)*(randn(nT,nR)+1i*randn(nT,nR));
elseif nargin < 2
    disp('[Message] Proposed_Precoder_Selection : Default set only nS');
    nS = 4;
end

%%
j= sqrt(-1);

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
    F4(:, :, i) = W(:, F4_matrix_order(i,:), i) / 2 ;
end

%%
x_set = -3:1:3;
delta_x_set = zeros( 1, (2 * nS - 1)^2); 
plus_cond = zeros( 1, (2 * nS - 1)^2); 
curr_delta_x_idx = 1;
max_of_min_radius = 0;

x_list_init= ones(nS,length(delta_x_set)).*100;
x_now_init = ones(nS,1).*100;

for rx = x_set
    for ix = x_set
        delta_x_set(curr_delta_x_idx) = rx + j * ix;
        if rx >= 0
            plus_cond(curr_delta_x_idx) = 1;
        end
        curr_delta_x_idx = curr_delta_x_idx + 1;
    end
end

real_delta_x_set=  real(delta_x_set);
imag_delta_x_set=  imag(delta_x_set);

if nS == 4
    precoder_index=[1,2,5,6,13];
elseif nS == 2
    precoder_index = 1:1:16;
end

%%

for precoder_idx=precoder_index
    searchOver = false;
    HF=H*F4(:,:,precoder_idx);
    [~,RR] = qr(HF);
    % U = chol((HF)' * HF);
    min_radius = min(vecnorm(RR)); 
    stage_idx = 1;
    Is_New_Stage = true;
    x_now = x_now_init;
    x_list = x_list_init;
    
    while ~searchOver
        if stage_idx == nS+1
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
                    x_list=x_list_init; 
                    x_now=x_now_init;
                else
                    stage_idx=stage_idx-1; 
                    Is_New_Stage=false;
                end
            end
        else
            if Is_New_Stage  
                % Refer "MIMO-OFDM Wireless Communications with MATLAB" 332p
                R_value_stage = RR(nS-stage_idx+1,nS-stage_idx+1);
                temp_sqrt=sqrt(min_radius^2-norm(RR(nS-stage_idx+2:end, nS-stage_idx+2:end)*...
                    x_now(nS-stage_idx+2:end))^2);
                temp_no_sqrt = (RR(nS-stage_idx+1, nS-stage_idx+2:end)* x_now(nS-stage_idx+2:end))/R_value_stage;

                % Compute the bounds
                 
                 range1 = -temp_sqrt/R_value_stage;
                 range2 = temp_sqrt/R_value_stage;
                 %if stage_idx is 1, delta_x_set->delta_x_set+
                if range1 > range2 
                    lower = range2;
                    upper = range1;
                elseif range1 < range2
                    lower = range1;
                    upper = range2;                    
                else
                    upper = range1;
                    lower = 0;
                end

                real_upper = upper - real(temp_no_sqrt);
                imag_upper = upper - imag(temp_no_sqrt);
                if stage_idx==1               
                    real_lower = 0;
                else
                    real_lower = lower - real(temp_no_sqrt);
                end
                imag_lower = lower - imag(temp_no_sqrt);
                % mod
                real_delta_x_cond = ((real_delta_x_set <= real_upper) & (real_delta_x_set >= real_lower));
                imag_delta_x_cond = ((imag_delta_x_set <= imag_upper) & (imag_delta_x_set >= imag_lower));
                delta_x_cond = real_delta_x_cond & imag_delta_x_cond;
                B = delta_x_set(delta_x_cond); 
                
                x_list(nS-stage_idx+1, 1:length(B))=B;
            end
            if ~any(x_list(nS-stage_idx+1,:) - 100) 
                %if the stage does not contain any candidate symbol
                if stage_idx==1
                    searchOver=true; 
                else
                    Is_New_Stage=false; 
                    stage_idx=stage_idx-1; 
                end
            else % candidates exist
                x_now(nS-stage_idx+1)=x_list(nS-stage_idx+1,1);
                x_list(nS-stage_idx+1,1:end-1) = x_list(nS-stage_idx+1,2:end);
                x_list(nS-stage_idx+1, end) = 100;
                Is_New_Stage=true; 
                stage_idx=stage_idx+1;
            end % if ~any(x_list(2*nS-stage_idx+1,:) - 100) 
        end
    end     % While
    if min_radius > max_of_min_radius
        max_of_min_radius = min_radius;
        opt_precoder_idx = precoder_idx;
    end
end

F = F4(:,:,opt_precoder_idx);