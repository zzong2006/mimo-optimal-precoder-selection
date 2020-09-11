% Implemented by Woosung 2020.Jan.14th
% Complex Lenstra-Lenstra Lovasz Algorithm (CLLL)

function [Q, R, T] = Complex_LLL(Hp, delta)
    %% Default Value
    if nargin < 1
        disp('[Message] Complex_LLL : Default set Hp, delta');
        % Default value : nR = nT = nS = 4
        Hp = 1/sqrt(4)*sqrt(1/2)*(randn(4,4)+1i*randn(4,4));
        delta = 0.99;
    elseif nargin < 2
        disp('[Message] Complex_LLL : Default set only delta');
        delta = 0.99;
    end
    
    [Q, R] = qr(Hp);
    M = size(Hp, 2);
    T = eye(M);
    m = 2;
    
    while m <= M
        for q = (m - 1): -1: 1
            u = round(R(q,m) / R(q,q));
            if u ~= 0
                R(1:q, m) = R(1:q, m) - u * R(1:q, q);
                T(:, m) = T(:, m) - u * T(:, q);
            end
        end
        if delta * abs(R(m-1, m-1))^2 > abs(R(m,m))^2 + abs(R(m-1, m))^2
            % Swap the (m-1)th and mth columns in R and T
            temp_col = R(:, m); R(:, m) = R(:, m-1); R(:, m-1) = temp_col;
            temp_col = T(:, m); T(:, m) = T(:, m-1); T(:, m-1) = temp_col;
            
            % Set Updelta
            alpha = R(m-1, m-1) / norm(R(m-1:m, m-1));
            beta = R(m, m-1) / norm(R(m-1:m, m-1));
            Uptheta = [ alpha' beta ; -beta alpha ];
            
            R(m-1:m, m-1:M) = Uptheta * R(m-1:m, m-1:M);
            Q(:, m-1:m) = Q(:, m-1:m) * Uptheta';
            m = max(m - 1, 2);
        else
            m = m + 1;
        end
    end
end