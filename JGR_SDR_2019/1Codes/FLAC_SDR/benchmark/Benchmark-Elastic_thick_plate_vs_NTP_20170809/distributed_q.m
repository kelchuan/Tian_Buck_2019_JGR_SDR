function [ L ] = distributed_q( L,N,wn,g,dx,rho_c,delta_rho_c, D)
% update the L matrix base on new wn
% distributed_q takes the wn as an argument and calculate the updated 
% wn 
% update distributed load q(x) in the coefficients based on calculated wn

%infill_flag get value of 1 at where wn(x)>=0 and get 0 when wn(x)<0
    infill_flag = (wn >= 0);
    infill_flag(N-3:N) = 0; % for boundary condition 1~4
    %R = infill_flag * 
    num_i = N;  
    num_j = N - 4;
%{
coef_i_minus_2 = D / dx^4;               % coefficient for w(i-2) (j == i)
coef_i_minus_1 = -4 * D / dx^4;          % coefficient for w(i-1) (i = j+1)
coef_i = 6 * D / dx^4 + delta_rho_c * g; % coeeficient for w(i)   (i = j+2)
coef_i_plus_1 = -4 * D / dx^4;           % coefficient for w(i+1) (i = j+3)
coef_i_plus_2 = D / dx^4;                % coefficient for w(i+2) (i = j+4)
%}
coef_i_minus_2 = 1;                         % coefficient for w(i-2) (j == i)
coef_i_minus_1 = -4;                        % coefficient for w(i-1) (i = j+1)
coef_i = 6 + delta_rho_c * g / (D / dx^4);  % coeficient for w(i)   (i = j+2)
coef_i_wn_positive = 6 + rho_c * g / (D / dx^4); % coef for wn(x) >= 0
coef_i_plus_1 = -4;                         % coefficient for w(i+1) (i = j+3)
coef_i_plus_2 = 1;                          % coefficient for w(i+2) (i = j+4)

for j = 1:1:num_j
    for i = 1:1:num_i
        if(i == j)
           L(j,i) = coef_i_minus_2;
        elseif(i==j+1)
            L(j,i) = coef_i_minus_1;
        elseif(i==j+2)
            if(infill_flag(i) == 0)
                L(j,i) = coef_i;
            elseif(infill_flag(i) == 1)
                L(j,i) = coef_i_wn_positive;
            end
        elseif(i==j+3)
            L(j,i) = coef_i_plus_1;
        elseif(i==j+4)
            L(j,i) = coef_i_plus_2;
        end
    end
end

end

