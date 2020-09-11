%filename:precoder_select_ML_ver4.m
function [F, F_idx]=precoder_select_ML_ver4(H,nS)
% Input parameters
%     H: Channel matrix (4 x 4 Complex matrix), 
%     k: delta_x range (default=4) 
%    nT: number of Tx antennas (default = 4)
%    getIndex : Selected precoder output type (1 == index, 0 == matrix)

% Output parameter
%    F : Precoder itself(when getIndex==0), or precoder index (when getIndex==1)
%    radius_scores : Calculated minimum radius of all(16) precoder

j = sqrt(-1);
u=[ 1,  1,  1,  1,  1,              1,             1,             1,              1,  1,  1,  1,  1,  1,  1, 1;
   -1, -j,  1,  j, (-1-j)/sqrt(2), (1-j)/sqrt(2), (1+j)/sqrt(2), (-1+j)/sqrt(2), -1, -j,  1,  j, -1, -1,  1, 1;
   -1,  1, -1,  1, -j,              j,            -j,             j,              1, -1,  1, -1, -1,  1, -1, 1;
   -1,  j,  1, -j, (1-j)/sqrt(2),  (-1-j)/sqrt(2),(-1+j)/sqrt(2),(1+j)/sqrt(2),   1, -j, -1,  j,  1, -1, -1, 1];

W=zeros(4,4,16);
W(:,:,1)=eye(4)-(2*u(:,1)*(u(:,1)'))/((u(:,1)')*u(:,1));
W(:,:,2)=eye(4)-(2*u(:,2)*(u(:,2)'))/((u(:,2)')*u(:,2));
W(:,:,3)=eye(4)-(2*u(:,3)*(u(:,3)'))/((u(:,3)')*u(:,3));
W(:,:,4)=eye(4)-(2*u(:,4)*(u(:,4)'))/((u(:,4)')*u(:,4));
W(:,:,5)=eye(4)-(2*u(:,5)*(u(:,5)'))/((u(:,5)')*u(:,5));
W(:,:,6)=eye(4)-(2*u(:,6)*(u(:,6)'))/((u(:,6)')*u(:,6));
W(:,:,7)=eye(4)-(2*u(:,7)*(u(:,7)'))/((u(:,7)')*u(:,7));
W(:,:,8)=eye(4)-(2*u(:,8)*(u(:,8)'))/((u(:,8)')*u(:,8));
W(:,:,9)=eye(4)-(2*u(:,9)*(u(:,9)'))/((u(:,9)')*u(:,9));
W(:,:,10)=eye(4)-(2*u(:,10)*(u(:,10)'))/((u(:,10)')*u(:,10));
W(:,:,11)=eye(4)-(2*u(:,11)*(u(:,11)'))/((u(:,11)')*u(:,11));
W(:,:,12)=eye(4)-(2*u(:,12)*(u(:,12)'))/((u(:,12)')*u(:,12));
W(:,:,13)=eye(4)-(2*u(:,13)*(u(:,13)'))/((u(:,13)')*u(:,13));
W(:,:,14)=eye(4)-(2*u(:,14)*(u(:,14)'))/((u(:,14)')*u(:,14));
W(:,:,15)=eye(4)-(2*u(:,15)*(u(:,15)'))/((u(:,15)')*u(:,15));
W(:,:,16)=eye(4)-(2*u(:,16)*(u(:,16)'))/((u(:,16)')*u(:,16));


F4=zeros(4,4,16);
F4(:,:,1)=W(:,[1 2 3 4],1) / 2;
F4(:,:,2)=W(:,[1 2 3 4],2) / 2;
F4(:,:,3)=W(:,[3 2 1 4],3) / 2;
F4(:,:,4)=W(:,[3 2 1 4],4) / 2;
F4(:,:,5)=W(:,[1 2 3 4],5) / 2;
F4(:,:,6)=W(:,[1 2 3 4],6) / 2;
F4(:,:,7)=W(:,[1 3 2 4],7) / 2;
F4(:,:,8)=W(:,[1 3 2 4],8) / 2;
F4(:,:,9)=W(:,[1 2 3 4],9) / 2;
F4(:,:,10)=W(:,[1 2 3 4],10) / 2;
F4(:,:,11)=W(:,[1 3 2 4],11) / 2;
F4(:,:,12)=W(:,[1 3 2 4],12) / 2;
F4(:,:,13)=W(:,[1 2 3 4],13) / 2;
F4(:,:,14)=W(:,[1 3 2 4],14) / 2;
F4(:,:,15)=W(:,[3 2 1 4],15) / 2;
F4(:,:,16)=W(:,[1 2 3 4],16) / 2;

% if k==2
%     delta_x_set=[-1:1:1];
% elseif k==4
%     delta_x_set=[-3:1:3];
% elseif k==6
%     delta_x_set=[-7:1:7];
% end

delta_x_set=[-7:1:7];

YES=1; NO=0; 
termination_idx=100; 
x_list=ones(2*nS,length(delta_x_set)).*100;
x_now=ones(2*nS,1).*100;
precoder_score = zeros(5, 1);
max_of_min_radius=0;
precoder_index=[1,2,5,6,13];
i = 0;

for precoder_idx=precoder_index
    i = i + 1;
    if nS == 4 
        HF=H*F4(:,:,precoder_idx);
    end
    
    H_R =[real(HF) -(imag(HF)); imag(HF) real(HF)]; % complex system -> real system
    [~,R_R] = qr(H_R);              % QR Decomposition
    min_radius=min(vecnorm(R_R)); %initialize min_radius
    if min_radius<=max_of_min_radius 
        precoder_score(i) = min_radius;
        continue;
    end
    stage_idx=1; 
    Is_New_Stage=YES; 
    x_list=(x_list.*0)+100; 
    x_now=(x_now.*0)+100;
    
    while stage_idx~=termination_idx
        if stage_idx==2*nS+1 %if a vector of length 2xnT is found
            if prod(x_now==0)==1
                stage_idx=stage_idx-1; 
                Is_New_Stage=NO; % if metric_temp is zero vector
            else
                metric_temp=norm(R_R*x_now); % calculate the radius
                if metric_temp<=max_of_min_radius % if metric_temp is less than max_of_min_radius, stop and start next procoder_idx 
                    min_radius=metric_temp; 
                    break; 
                end
                if metric_temp<min_radius
                    min_radius=metric_temp;
                    stage_idx=1; 
                    Is_New_Stage=YES; 
                    x_list=(x_list.*0)+100; 
                    x_now=(x_now.*0)+100;
                else
                    stage_idx=stage_idx-1; 
                    Is_New_Stage=NO;
                end
            end
        else %if the length of the vector x_now is shorter than 2nT
            if Is_New_Stage==YES % In case of a new stage, we choose candidate symbols for the stage
                temp_sqrt=sqrt(min_radius-norm(R_R(2*nS-stage_idx+2:end,2*nS-stage_idx+2:end)*...
                    x_now(2*nS-stage_idx+2:end))^2);
                temp_no_sqrt=R_R(2*nS-stage_idx+1,2*nS-stage_idx+2:end)*...
                    x_now(2*nS-stage_idx+2:end);
                if R_R(2*nS-stage_idx+1,2*nS-stage_idx+1)>0
                    bound_lower=(-temp_sqrt-temp_no_sqrt)/R_R(2*nS-stage_idx+1,2*nS-stage_idx+1);
                    bound_upper=(temp_sqrt-temp_no_sqrt)/R_R(2*nS-stage_idx+1,2*nS-stage_idx+1);
                else
                    bound_lower=(temp_sqrt-temp_no_sqrt)/R_R(2*nS-stage_idx+1,2*nS-stage_idx+1);
                    bound_upper=(-temp_sqrt-temp_no_sqrt)/R_R(2*nS-stage_idx+1,2*nS-stage_idx+1);
                end
                [~,B]=find((bound_lower<delta_x_set)&(delta_x_set<bound_upper));
                if stage_idx==1
                    B = B(ceil(length(B)/2):end); 
                end %if stage_idx is 1, delta_x_set->delta_x_set+
                x_list(2*nS-stage_idx+1,1:length(B))=delta_x_set(B);
            end
            if isempty(find(x_list(2*nS-stage_idx+1,:)~=100)) 
                %if the stage does not contain any candidate symbol
                if stage_idx==1
                    stage_idx=termination_idx; 
                else
                    Is_New_Stage=NO; 
                    stage_idx=stage_idx-1; 
                end
            else % candidates exist
                x_now(2*nS-stage_idx+1)=x_list(2*nS-stage_idx+1,1);
                x_list(2*nS-stage_idx+1,:)=[x_list(2*nS-stage_idx+1,[2:end]) 100];
                Is_New_Stage=YES; 
                stage_idx=stage_idx+1;
            end % if length(find(x_list(2*nT-stage_idx+1,:)~=100))==0
        end % if stage_idx==2*nT+1
    end % while
    if min_radius > max_of_min_radius
        max_of_min_radius = min_radius;
        opt_precoder_idx = precoder_idx;
    end
    precoder_score(i) =  min_radius;
end % for
disp(precoder_score);
if nS == 4
    F = F4(:,:,opt_precoder_idx);
    F_idx = opt_precoder_idx;
end

% radius_scores = precoder_score;