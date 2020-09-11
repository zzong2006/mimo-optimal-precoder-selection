function [X_hat]=Sphere_Decoder(y,H,nT, L, var_n)
% Input parameters
%     y: received signal nRx1, H: Channel matrix, nRxnT, nT: number of Tx
%     antennas, L : 4 => 16 QAM, 2 => 8 QAM 
% Output parameter
%    X_hat: estimated signal, nTx1
if nargin < 1
    disp('[Message] Sphere Decoding : Default set H, K, delta');
    nT = 2; nR = 2; 
    L = 4;
    y = randn(nR, 1);
    H = 1/sqrt(L)*sqrt(1/2)*(randn(nR, nT)+1i*randn(nR,nT));
    var_n=10^(-17/10)/L;
elseif nargin < 2
    disp('Sphere Decoding Error');
elseif nargin < 3
    disp('Sphere Decoding Error');
end
%MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd
switch L
    case 2, real_const=(-1:2:1)/sqrt(2);
    case 4, real_const=(-3:2:3)/sqrt(10);
    case 6, real_const=(-7:2:7)/sqrt(42); 
end
YES=1; NO=0; termination_idx=100; j=sqrt(-1);
y_R =[real(y); imag(y)]; 
H_R =[real(H) -(imag(H)); imag(H) real(H)]; % complex system -> real system
x_list=zeros(2*nT,length(real_const));
x_now=zeros(2*nT,1); 

[~,R_R] = qr(H_R); 
x_hat=pinv(H_R)*y_R; % unconstrained ML solution (=ZF solution)
x_MMSE=pinv(H_R'*H_R+var_n/2*eye(2*nT))*H_R'*y_R; 

temp_sliced=transpose(QAM_slicer(x_MMSE(1:nT)+j*x_MMSE(nT+1:2*nT),L)); % using the complex slicer
x_pre=[real(temp_sliced); imag(temp_sliced)]; 
radius_squared=norm(R_R*(x_pre-x_hat))^2;
stage_idx=1; Is_New_Stage=YES;

while (stage_idx~=termination_idx)
   if stage_idx==2*nT+1 %if a vector of length 2xnT is found
       metric_temp=norm(R_R*(x_now-x_hat))^2; %calculate the radius
       if metric_temp<radius_squared %new candidate vector has smaller metric, then restart
           x_pre=x_now; 
           radius_squared=metric_temp;  
           stage_idx=1; 
           Is_New_Stage=YES; 
           x_list=zeros(2*nT,length(real_const)); 
           x_now=zeros(2*nT,1); %initialization
       else stage_idx=stage_idx-1; Is_New_Stage=NO; end %the found vector is the x_pre
   else %if the length of the vector x_now is shorter than 2nT
       if Is_New_Stage==YES % In case of a new stage, we choose candidate symbols for the stage
           temp_sqrt=sqrt(radius_squared-norm(R_R(2*nT-stage_idx+2:end,2*nT-stage_idx+2:end)*...
               [x_now(2*nT-stage_idx+2:end)-x_hat(2*nT-stage_idx+2:end)])^2);
           temp_no_sqrt=R_R(2*nT-stage_idx+1,2*nT-stage_idx+2:end)*...
               (x_now(2*nT-stage_idx+2:end)-x_hat(2*nT-stage_idx+2:end));
           if R_R(2*nT-stage_idx+1,2*nT-stage_idx+1)>0
               bound_lower=(-temp_sqrt-temp_no_sqrt)/R_R(2*nT-stage_idx+1,2*nT-stage_idx+1)+x_hat(2*nT-stage_idx+1);
               bound_upper=(temp_sqrt-temp_no_sqrt)/R_R(2*nT-stage_idx+1,2*nT-stage_idx+1)+x_hat(2*nT-stage_idx+1);
           else
               bound_lower=(temp_sqrt-temp_no_sqrt)/R_R(2*nT-stage_idx+1,2*nT-stage_idx+1)+x_hat(2*nT-stage_idx+1);
               bound_upper=(-temp_sqrt-temp_no_sqrt)/R_R(2*nT-stage_idx+1,2*nT-stage_idx+1)+x_hat(2*nT-stage_idx+1);
           end
           [~,B]=find((bound_lower<real_const)&(real_const<bound_upper));
           x_list(2*nT-stage_idx+1,1:length(B))=real_const(B);   
       end % if Is_New_Stage==YES
       if length(find(x_list(2*nT-stage_idx+1,:)~=0))==0 %if the stage does not contain any candidate symbol
           if stage_idx==1; stage_idx=termination_idx; else Is_New_Stage=NO; stage_idx=stage_idx-1; end 
       else % candidates exist
           x_now(2*nT-stage_idx+1)=x_list(2*nT-stage_idx+1,1); 
           x_list(2*nT-stage_idx+1,:)=[x_list(2*nT-stage_idx+1,[2:end]) 0]; 
           Is_New_Stage=YES; stage_idx=stage_idx+1;
       end % if length(find(x_list(2*nT-stage_idx+1,:)~=0))==0
   end % if stage_idx==2*nT+1, else
end %while
X_hat=x_pre(1:nT)+j*x_pre(nT+1:2*nT);