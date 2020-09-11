function [X_hat]=Sphere_Decoder_Complex(y,H,nT,L,var_n)
% Input parameters
%     y: received signal nRx1, H: Channel matrix, nRxnT, nT: number of Tx antennas
% Output parameter
%    X_hat: estimated signal, nTx1
if nargin < 1
    rng('default');
    rng(1111);
    disp('[Message] Sphere Decoding : Default set H, K, delta');
    nT = 2; nR = 4; 
    L = 4;
    y = randn(nR, 1);
    H = 1/sqrt(L)*sqrt(1/2)*(randn(nR, nT)+1i*randn(nR,nT));
    var_n=10^(-17/10)/L;
elseif nargin < 2
    disp('Sphere Decoding Error');
elseif nargin < 3
    disp('Sphere Decoding Error');
end

% 2020 Woosung

j=sqrt(-1);

switch L
    case 2
        real_const=[-1:2:1]/sqrt(2);
    case 4
        complex_const = zeros(1, L^2); curr_idx = 1;
        real_const=(-3:2:3)/sqrt(10);
        for xi = real_const
            for xy = real_const
                complex_const(curr_idx) = xi + xy * j;
                curr_idx = curr_idx + 1;
            end
        end
    case 6
        real_const=[-7:2:7]/sqrt(42); 
end

x_list = zeros(nT,length(complex_const));
x_now=zeros(nT,1); 

[~,R_R] = qr(H); 
size_H = size(H);
if size_H(1) ~= size_H(2)
    R_R = R_R(1:min(size_H),:);
end
R_R = complex(R_R);
R_R = diag(sign(diag(R_R))) * R_R;
% make diagonal elements of R matrix to positive
R_R = diag(sign(diag(R_R))) * R_R;
x_hat=H\y; % unconstrained ML solution (=ZF solution)
x_MMSE=(H'*H+var_n/2*eye(nT))\H'*y; 

x_pre=transpose(QAM_slicer(x_MMSE, L)); % using the complex slicer
radius_squared=norm(R_R*(x_pre-x_hat))^2;
stage_idx=1; Is_New_Stage=true;
searchOver = false;
candidate_count = zeros(nT, 1);
candidate_idx = ones(nT, 1);

while ~searchOver
   if stage_idx==nT + 1 %if a vector of length 2xnT is found
       metric_temp = norm(R_R*(x_now-x_hat))^2; %calculate the radius
       if metric_temp < radius_squared %new candidate vector has smaller metric, then restart
           x_pre=x_now; 
           radius_squared=metric_temp;  
           stage_idx=1; 
           Is_New_Stage=true; 
           x_list=zeros(nT,length(complex_const)); 
           x_now=zeros(nT,1); %initialization
       else
           stage_idx=stage_idx-1; 
           Is_New_Stage = false; 
       end %the found vector is the x_pre
   else %if the length of the vector x_now is shorter than 2nT
       if Is_New_Stage % In case of a new stage, we choose candidate symbols for the stage
           R_value_stage = R_R(nT-stage_idx+1,nT-stage_idx+1);
           RR_mul_x = R_R(end-stage_idx+1:end, end-stage_idx+2:end) * (x_now(end-stage_idx+2:end)-x_hat(end-stage_idx+2:end));
           
           if stage_idx == 1
               middle = radius_squared;
               coord = 0;
           else
              middle = radius_squared-sum(abs(RR_mul_x(2:end)).^2);
              coord = RR_mul_x(1) / R_value_stage;   
           end
           squared_const = abs(complex_const - x_hat(end-stage_idx + 1) + coord).^2 ;
           range = middle / (R_value_stage^2);
           const_cond = squared_const <= range;

           len_B = sum(const_cond);
           x_list(nT-stage_idx+1,1:len_B)=complex_const(const_cond);   
           candidate_count(stage_idx) = len_B;
           candidate_idx(stage_idx) = 1;
       end % if Is_New_Stage==YES
       if candidate_count(stage_idx) == 0 %if the stage does not contain any candidate symbol
           if stage_idx==1
               searchOver = true;
           else 
               Is_New_Stage=false; 
               stage_idx=stage_idx-1; 
           end 
       else % candidates exist
           x_now(end-stage_idx+1)=x_list(end-stage_idx+1, candidate_idx(stage_idx));
            Is_New_Stage=true; 
            candidate_idx(stage_idx) = candidate_idx(stage_idx) + 1;
            candidate_count(stage_idx) = candidate_count(stage_idx) - 1;
            stage_idx=stage_idx+1;
       end % if length(find(x_list(2*nT-stage_idx+1,:)~=0))==0
   end % if stage_idx==2*nT+1, else
end %while
X_hat=x_pre;