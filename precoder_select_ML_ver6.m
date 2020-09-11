function [F, F_index]=precoder_select_ML_ver6(H,nS)
% Input parameters
%     H: Channel matrix (4 x 4 Complex matrix), 
%    nS: number of Tx antennas (default = 4)

% Output parameter
%    F : Precoder itself(when getIndex==0), or precoder index (when getIndex==1)
%    radius_scores : Calculated minimum radius of all(16) precoder

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

delta_x_set=[-3 -2 -1 0 1 2 3]; plus_cond = logical([0 0 0 1 1 1 1]);
% delta_x_set=[0 1 -1 2 -2 3 -3]; plus_cond = [1 1 0 1 0 1 0];

YES=1; NO=0; 
termination_idx=100; 
x_list_init=ones(2*nS,length(delta_x_set)).*100;
x_now_init =ones(2*nS,1).*100;
max_of_min_radius=0;

if nS == 4
    precoder_index=[1,2,5,6,13];
elseif nS == 2
    precoder_index = 1:1:16;
end


for precoder_idx=precoder_index
    candidate_count = zeros(nS * 2, 1);
    candidate_idx = ones(nS, 1);
    zero_count = 0;
    
    if nS == 4
        HF=H*F4(:,:,precoder_idx);
    elseif nS == 2
        HF=H*F2(:,:,precoder_idx);
    end
    H_R =[real(HF) -(imag(HF)); imag(HF) real(HF)]; % complex system -> real system
    [~,R_R] = qr(H_R);              % QR Decomposition
    min_radius = min(vecnorm(R_R)); %initialize min_radius
    
    if min_radius<=max_of_min_radius 
        continue;
    end
    
    stage_idx=1; 
    Is_New_Stage=YES; 
    x_list=x_list_init; 
    x_now=x_now_init;
    
    while stage_idx~=termination_idx
        if stage_idx==2*nS+1 %if a vector of length 2xnT is found
            if zero_count == 2 * nS
                stage_idx=stage_idx-1;
                zero_count = zero_count - 1;
                Is_New_Stage=NO; % if metric_temp is zero vector
            else
                metric_temp=norm(R_R*x_now); % calculate the radius
                % if metric_temp is less than max_of_min_radius, stop and start next procoder_idx 
                if metric_temp<=max_of_min_radius 
                    min_radius=metric_temp; 
                    break; 
                end
                if metric_temp<min_radius
                    min_radius=metric_temp;
                    stage_idx=1; 
                    Is_New_Stage=YES; 
                    zero_count = 0;
                    x_list=x_list_init; 
                    x_now=x_now_init;
                else
                    if x_now(1) == 0
                        zero_count = zero_count - 1;
                    end
                    stage_idx=stage_idx-1; 
                    Is_New_Stage=NO;
                end
            end
        else %if the length of the vector x_now is shorter than 2nT
            if Is_New_Stage==YES % In case of a new stage, we choose candidate symbols for the stage
                %tic
                RR_mul_x_now = R_R(2 *nS-stage_idx+1:end, 2 *nS-stage_idx+2:end)*...
                    x_now(2 *nS-stage_idx+2:end);
                if stage_idx == 1 
                    middle = min_radius^2;
                else
                    middle = min_radius^2-sum(abs(RR_mul_x_now(2:end)).^2);                  
                end

                R_value_stage = R_R(2*nS-stage_idx+1,2*nS-stage_idx+1);
                if R_value_stage > 0
                    bound_lower=(-sqrt(middle)-RR_mul_x_now(1))/R_value_stage;
                    bound_upper=(sqrt(middle)-RR_mul_x_now(1))/R_value_stage;
                else
                    bound_lower=(sqrt(middle)-RR_mul_x_now(1))/R_value_stage;
                    bound_upper=(-sqrt(middle)-RR_mul_x_now(1))/R_value_stage;
                end
                
                delta_x_cond = (bound_lower <= delta_x_set)&(delta_x_set <= bound_upper);
                if stage_idx==1               
                    delta_x_cond(1:3) = false;
                end %if stage_idx is 1, delta_x_set->delta_x_set+

                len_B = sum(delta_x_cond);
                x_list(end-stage_idx+1, 1:len_B) = delta_x_set(delta_x_cond);
                candidate_count(stage_idx) = len_B;
                candidate_idx(stage_idx) = 1;
                %time_temp = time_temp + toc;
                % time_curr = time_curr  + 1;
            end
            if candidate_count(stage_idx) == 0 
                %if the stage does not contain any candidate symbol
                if stage_idx==1
                    stage_idx=termination_idx; 
                else
                    Is_New_Stage=NO; 
                    stage_idx=stage_idx-1; 
                    if x_now(end - stage_idx + 1) == 0
                        zero_count = zero_count - 1;
                    end

                end
            else % candidates exist
                candidate_entry = x_list(end-stage_idx+1,candidate_idx(stage_idx));
                x_now(end -stage_idx+1)= candidate_entry ;
                if candidate_entry == 0
                    zero_count = zero_count + 1;
                end
                Is_New_Stage=YES; 
                candidate_idx(stage_idx) = candidate_idx(stage_idx) + 1;
                candidate_count(stage_idx) = candidate_count(stage_idx) - 1;
                stage_idx=stage_idx+1;

            end % if ~any(x_list(2*nS-stage_idx+1,:) - 100) 
        end % if stage_idx==2*nS+1
    end % while
    if min_radius > max_of_min_radius
        max_of_min_radius = min_radius;
        opt_precoder_idx = precoder_idx;
    end
    % precoder_score(precoder_idx) = min_radius;
end % for

if nS == 4
    F = F4(:,:,opt_precoder_idx);
    F_index = opt_precoder_idx;
elseif nS == 2
    F = opt_precoder_idx;
end
% disp(time_temp / time_curr);
% radius_scores = precoder_score;