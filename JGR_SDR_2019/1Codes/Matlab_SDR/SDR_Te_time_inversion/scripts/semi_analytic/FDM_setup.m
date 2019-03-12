%----------------------------------------------------------------
% Finite difference solution to D * d4w/dx4 + q(x) = 0
% Implicit solution L * W  = R
% L are the corresponding coefficients of w(i-2 ~ i+2)
% W are all the w
% R are the right hand side which are the extra loads
%----------------------------------------------------------------

% Return the LHS Lm matrix and RHS R matrix and whether lava sea or not
function [Lm,R,index_lava] = FDM_setup(N,D,dx,dV_0)

%----------------------------------------------------------------
% Lava sea or not
%----------------------------------------------------------------
index_lava = 0;
% index_lava = 1 means lava sea where:
                        % q(x) = (rho_i - rho_c) * g * w(x) for any w(x)
% index_lava = 0 means:
                        % q(x) = (rho_i - rho_c) * g at w(x) < 0
                        % q(x) = - rho_c * g at w(x) >= 0
%----------------------------------------------------------------

%----------------------------------------------------------------
% Setup the R matrix (N-5+1+4, 1)  % THe added 4 is for 4 BCs
%----------------------------------------------------------------
R = zeros(N-5+1+4, 1);

%----------------------------------------------------------------
% setup the L matrix (coefficients of w)
% D/dx^4 * (w(i+2) - 4*w(i+1) + 6*w(i) - 4*w(i-1) + w(i-2)) + ...
% delta_rho_c * g * w(i) = 0
%----------------------------------------------------------------
num_i = N;  
num_j = N - 4;
Lm = zeros(N, N);

%------------------------------------------
% Boundary conditions
%------------------------------------------
index_broken = 1; % if index_broken = 1, then it is calculating broken plate
                  % if index_broken != 1, it is for continuous plate

% BC1 approximate line load due to denser dike at the center 
% R(1) = -dV_0 / 2 / (D / dx^4);
% Rather than applying on the RHS of the first equation, we use BCs of 
% V = dV_0 at x = 0

if (index_broken == 1)
    % The remaining four BCs can be added to the L and R matrices
    % The row number of the four BCs is random and interchangeable
    %---------------------
    % BC1: wn(N) = 0;
    %---------------------
    Lm(num_j+1,num_i) = 1;
    %---------------------
    % BC2: w'(inf) = 0;
    % wn(N-1) = wn(N);
    %---------------------
    Lm(num_j+2,num_i) = 1;
    Lm(num_j+2,num_i-1) = -1;
    %---------------------
    % BC3: M = 0 at x = 0  --> d2w/dx2|(x=0) = 0
    % d2w/dx2 = (wn(i+1) - 2wn(i) + wn(i-1)) / (dx)^2;
    % M = -D * d2w/dx2 = -D * (wn(i+1) - 2wn(i) + wn(i-1)) / (dx)^2;
    %---------------------
    Lm(num_j+3, 1) = 1;
    Lm(num_j+3, 2) = -2;
    Lm(num_j+3, 3) = 1;
    %---------------------
    % BC4 V = dM/dx = -D * d3w/dx3 = V_0/2 = dV_0 at x = 0
    % d3w/dx3|(i-0.5) = (d2w/dx2|i - d2w/dx2|i-1) / dx
    %               = (wn(i+1) -3*wn(i)+3*wn(i-1)-wn(i-2)) / (dx)^3
    % -D * (-wn(1) + 3 * wn(2) - 3 * wn(3) + wn(4)) / dx^3;
    %---------------------
    Lm(num_j+4, 1) = -1;
    Lm(num_j+4, 2) = 3;
    Lm(num_j+4, 3) = -3;         
    Lm(num_j+4, 4) = 1;
    R(num_j+4,1) = dV_0 / (-D/dx^3); % value at the RHS
    
else % Continuous plate (index_broken != 1)
    %---------------------
    %BC1
    %wn(N) = 0;
    %---------------------
    Lm(num_j+1,num_i) = 1;
    %---------------------
    %BC2
    %wn(N-1) = wn(N);
    %---------------------
    Lm(num_j+2,num_i) = 1;
    Lm(num_j+2,num_i-1) = -1;
    %---------------------
    % BC3 wn' = 0 at x = 0
    % wn(1) = wn(2);
    %---------------------
    Lm(num_j+3, 1) = 1;
    Lm(num_j+3, 2) = -1;
    %---------------------
    % BC4 V = 0 at x = 0
    %wn(4) = wn(1) - 3 * wn(2) + 3 * wn(3);
    %---------------------
    Lm(num_j+4, 1) = -1;
    Lm(num_j+4, 2) = 3;
    Lm(num_j+4, 3) = -3;
    Lm(num_j+4, 4) = 1;
    R(num_j+4,1) = dV_0 / (-D/dx^3); % value at the RHS
end




