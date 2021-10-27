function [A, rhs, u_initial,P,L,x,y] = Make_Diffusion_Operator(m,n,Prol_Type,distributionType)

global T_FLAG B_FLAG L_FLAG R_FLAG   % 0 for Dirichlet conditions (function given on boundaries),
T_FLAG = 0;
B_FLAG = 0;
L_FLAG = 0;
R_FLAG = 0;
A = zeros(m,n,9);

% bilinear finite element - nine points.
if strcmp(distributionType,'lognormal')
    D_mat = exp(randn(m+1,n+1));
elseif strcmp(distributionType, 'uniform')
    D_mat = rand(m+1,n+1);%
end
A(:,:,2) = -1/6*(D_mat(1:end-1,1:end-1)+D_mat(1:end-1,2:end));
A(:,:,3) = -1/3*D_mat(1:end-1,2:end);
A(:,:,4) = -1/6*(D_mat(1:end-1,2:end)+D_mat(2:end,2:end));
A(:,:,5) = -1/3*D_mat(2:end,2:end);
A(:,:,6) = -1/6*(D_mat(2:end,2:end)+D_mat(2:end,1:end-1));
A(:,:,7) = -1/3*D_mat(2:end,1:end-1);
A(:,:,8) = -1/6*(D_mat(2:end,1:end-1)+D_mat(1:end-1,1:end-1));
A(:,:,9) = -1/3*D_mat(1:end-1,1:end-1);
A(:,:,1) = 2/3*(D_mat(1:end-1,1:end-1)+D_mat(1:end-1,2:end)+D_mat(2:end,1:end-1)+D_mat(2:end,2:end));

rhs = zeros(m,n);
u_initial = rand(m,n);


% Homogeneous Neumann boundary conditions
if (T_FLAG == 1)
    A(1,:,1) = A(1,:,1) + A(1,:,2);
end % if
if (B_FLAG == 1)
    A(end,:,1) = A(end,:,1) + A(end,:,6);
end % if
if (L_FLAG == 1)
    A(:,1,1) = A(:,1,1) + A(:,1,8);
end % if
if (R_FLAG == 1)
    A(:,end,1) = A(:,end,1) + A(:,end,4);
end % if

% This yields homogeneous boundary conditions
A(1,:,9) = 0;
A(1,:,2:3) = 0;
A(end,:,5:7) = 0;
A(:,1,7:9) = 0;
A(:,end,3:5) = 0;

if nargout>3
    % buid the prolongation and stencil via Galerkin formulation.
    %Prol_Type = 0;
    y = m;
    x = n;
    while((y(end) > 1) && (x(end) > 1)) % design coarse levels.
        y = [y floor(y(end)/2)];
        x = [x floor(x(end)/2)];
    end % while
    levels = length(y)-1;
    L = cell(levels,1); % Operator
    P = cell(levels,1); % Prolongation matrices
    % ############################################
    % Initialize; add "padding" of a line of zeros around the domain, so that we
    % don't need to deal with the boundary regions separately.
    L{1,1} = zeros(m+2,n+2,9);
    L{1,1}(:,:,1) = 1;
    L{1,1}(2:end-1,2:end-1,:) = A;
    
    for grid = 1:levels-1
        if (Prol_Type == -1)
            P{grid,1} = Bilin_Prol(y(grid), x(grid), y(grid+1), x(grid+1)); % Bi-linear interpolation
        elseif (Prol_Type == 0)
            P{grid,1} = Bbox_Prol(y(grid), x(grid), y(grid+1), x(grid+1), L{grid,1}); % Blackbox prolongation
        elseif (Prol_Type == 1)
            P{grid,1} = Triangle_Prol_half(y(grid), x(grid), y(grid+1), x(grid+1), L{grid,1}); % Triangle-Grid prolongation
        elseif (Prol_Type == 2)
            P{grid,1} = Triangle_Prol_one(y(grid), x(grid), y(grid+1), x(grid+1), L{grid,1}); % Triangle-Grid prolongation
        else
            P{grid,1} = Triangle_Prol_zero(y(grid), x(grid), y(grid+1), x(grid+1), L{grid,1}); % Triangle-Grid prolongation
        end
        
        L{grid+1,1} = Galerkin(P{grid,1}, L{grid,1}, y(grid), x(grid), ...
            y(grid+1), x(grid+1)); % Galerkin coarsening
    end
    
end
return

function P = Bbox_Prol(m, n, m_c, n_c, L) % The black-box prolongation operator
% The prolongation operator is a matrix of size  (m_c,n_c,9) plus pad.
% P(15,23,1), for example, represents the coefficient multiplying the corresponding
% coarse-grid variable (15,23), when interpolating to the fine-grid
% variable that coincides with it, while P(15,23,3) is the coefficient
% multiplying this variable when interpolating to the fine-grid variable
% that is shifted -1 in the i direction (up), and +1 in the j direction (right).
P = zeros(m_c+2,n_c+2,9);
P(2:end-1,2:end-1,1) = 1;
P(2:end-1,2:end-1,2) = ...
    -(L(2:2:end-2,3:2:end-1,5) + L(2:2:end-2,3:2:end-1,6) + L(2:2:end-2,3:2:end-1,7)) ./ ...
    (L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,1) + L(2:2:end-2,3:2:end-1,8));
P(2:end-1,2:end-1,4) = ...
    -(L(3:2:end-1,4:2:end,7) + L(3:2:end-1,4:2:end,8) + L(3:2:end-1,4:2:end,9)) ./ ...
    (L(3:2:end-1,4:2:end,6) + L(3:2:end-1,4:2:end,1) + L(3:2:end-1,4:2:end,2));
P(2:end-1,2:end-1,6) = ...
    -(L(4:2:end,3:2:end-1,9) + L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,2)) ./ ...
    (L(4:2:end,3:2:end-1,8) + L(4:2:end,3:2:end-1,1) + L(4:2:end,3:2:end-1,4));
P(2:end-1,2:end-1,8) = ...
    -(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,4) + L(3:2:end-1,2:2:end-2,5)) ./ ...
    (L(3:2:end-1,2:2:end-2,2) + L(3:2:end-1,2:2:end-2,1) + L(3:2:end-1,2:2:end-2,6));

P(2:end-1,2:end-1,3) = -(L(2:2:end-2,4:2:end,7) + P(2:end-1,2:end-1,2) .* ...
    L(2:2:end-2,4:2:end,8) + P(2:end-1,2:end-1,4) .* L(2:2:end-2,4:2:end,6)) ./ ...
    L(2:2:end-2,4:2:end,1);
P(2:end-1,2:end-1,5) = -(L(4:2:end,4:2:end,9) + P(2:end-1,2:end-1,4) .* ...
    L(4:2:end,4:2:end,2) + P(2:end-1,2:end-1,6) .* L(4:2:end,4:2:end,8)) ./ ...
    L(4:2:end,4:2:end,1);
P(2:end-1,2:end-1,7) = -(L(4:2:end,2:2:end-2,3) + P(2:end-1,2:end-1,6) .* ...
    L(4:2:end,2:2:end-2,4) + P(2:end-1,2:end-1,8) .* L(4:2:end,2:2:end-2,2)) ./ ...
    L(4:2:end,2:2:end-2,1);
P(2:end-1,2:end-1,9) = -(L(2:2:end-2,2:2:end-2,5) + P(2:end-1,2:end-1,2) .* ...
    L(2:2:end-2,2:2:end-2,4) + P(2:end-1,2:end-1,8) .* L(2:2:end-2,2:2:end-2,6)) ./ ...
    L(2:2:end-2,2:2:end-2,1);

% Refrain from interpolating to pads when m and/or n are even.
if (mod(m,2) == 0)
    P(end-1,2:end-1,5:7) = 0;
end % if
if (mod(n,2) == 0)
    P(2:end-1,end-1,3:5) = 0;
end % if

return

function P = Bilin_Prol(m, n, m_c, n_c) % Bilinear interpolation

global T_FLAG B_FLAG L_FLAG R_FLAG  % 0 for Dirichlet conditions (function given on boundaries),
% 1 for Neumann (outward derivatives).
P = zeros(m_c+2,n_c+2,9);
P(2:end-1,2:end-1,1) = 1;
P(2:end-1,2:end-1,2) = 0.5;
P(2:end-1,2:end-1,3) = 0.25;
P(2:end-1,2:end-1,4) = 0.5;
P(2:end-1,2:end-1,5) = 0.25;
P(2:end-1,2:end-1,6) = 0.5;
P(2:end-1,2:end-1,7) = 0.25;
P(2:end-1,2:end-1,8) = 0.5;
P(2:end-1,2:end-1,9) = 0.25;

% Modify near boundaries in case of Neumann boundary conditions
if (T_FLAG == 1)
    P(2,2:end-1,2:3) = 2*P(2,2:end-1,2:3);
    P(2,2:end-1,9) = 2*P(2,2:end-1,9);
end % if
if (B_FLAG == 1)
    P(end-1,2:end-1,5:7) = 2*P(end-1,2:end-1,5:7);
end % if
if (L_FLAG == 1)
    P(2:end-1,2,7:9) = 2*P(2:end-1,2,7:9);
end % if
if (R_FLAG == 1)
    P(2:end-1,end-1,3:5) = 2*P(2:end-1,end-1,3:5);
end % if

% Refrain from interpolating to boundary "pads"
if (mod(m,2) == 0)
    P(end-1,2:end-1,5:7) = 0;
end % if
if (mod(n,2) == 0)
    P(2:end-1,end-1,3:5) = 0;
end % if

return

function L_c = Galerkin(P, L, m, n, m_c, n_c) % Galerkin coarsening

L_c = zeros(m_c+2,n_c+2,9);
L_c(:,:,1) = 1;

% We have to do this nine times to separate the contributions. We put 1 at
% some coarse-grid variables, prolong, apply L on the fine grid, and restrict.
% The coarse-grid values thus obtained are the required coefficients.
for i = 2:4
    for j = 2:4
        u_c = zeros(m_c+2, n_c+2);
        if ((m_c >= i-1) && (n_c >= j-1))
            u_c(i:3:end-1,j:3:end-1) = 1;
            u = Prolong(P, u_c, m, n);
            u = Lu(L, u);
            u_c = Restrict(P, u, m_c, n_c);
            L_c(i:3:end-1,j:3:end-1,1) = u_c(i:3:end-1,j:3:end-1);
            L_c(i-1:3:end,j:3:end,6) = u_c(i-1:3:end,j:3:end);
            L_c(i:3:end,j-1:3:end,4) = u_c(i:3:end,j-1:3:end);
            L_c(i-1:3:end,j-1:3:end,5) = u_c(i-1:3:end,j-1:3:end);
            if (m_c >= i) % Avoid stepping out of bounds
                L_c(i+1:3:end,j:3:end,2) = u_c(i+1:3:end,j:3:end);
                L_c(i+1:3:end,j-1:3:end,3) = u_c(i+1:3:end,j-1:3:end);
            end % if
            if (n_c >= j) % Avoid stepping out of bounds
                L_c(i-1:3:end,j+1:3:end,7) = u_c(i-1:3:end,j+1:3:end);
                L_c(i:3:end,j+1:3:end,8) = u_c(i:3:end,j+1:3:end);
            end % if
            if ((n_c >= j) && (m_c >= i)) % Avoid stepping out of bounds
                L_c(i+1:3:end,j+1:3:end,9) = u_c(i+1:3:end,j+1:3:end);
            end % if
        end % if
    end % for
end % for

return;

function P = Triangle_Prol_half(m, n, m_c, n_c, L) % Prolongation from two near neighbors -
% suitable for the case of triangular grids

% The prolongation operator is a matrix of size  (m_c,n_c,9) plus pad.
% P(15,23,1), for example, represents the coefficient multiplying the corresponding
% coarse-grid variable (15,23), when interpolating to the fine-grid
% variable that coincides with it, while P(15,23,3) is the coefficient
% multiplying this variable when interpolating to the fine-grid variable
% that is shifted -1 in the i direction (up), and +1 in the j direction (right).
P = zeros(m_c+2,n_c+2,9);
P(2:end-1,2:end-1,1) = 1;
P(2:end-1,2:end-1,2) = ...
    -(L(2:2:end-2,3:2:end-1,6) + 0.5*(L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,7))) ./ ...
    (L(2:2:end-2,3:2:end-1,1)  + 0.5*(L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,7) + ...
    L(2:2:end-2,3:2:end-1,3) + L(2:2:end-2,3:2:end-1,8)));
P(2:end-1,2:end-1,3) = ...
    -(L(2:2:end-2,4:2:end,7) + 0.5*(L(2:2:end-2,4:2:end,6) + L(2:2:end-2,4:2:end,8))) ./ ...
    (L(2:2:end-2,4:2:end,1)  + 0.5*(L(2:2:end-2,4:2:end,6) + L(2:2:end-2,4:2:end,8) + ...
    L(2:2:end-2,4:2:end,2) + L(2:2:end-2,4:2:end,4)));
P(2:end-1,2:end-1,4) = ...
    -(L(3:2:end-1,4:2:end,8) + 0.5*(L(3:2:end-1,4:2:end,2) + L(3:2:end-1,4:2:end,7))) ./ ...
    (L(3:2:end-1,4:2:end,1)  + 0.5*(L(3:2:end-1,4:2:end,2) + L(3:2:end-1,4:2:end,7) + ...
    L(3:2:end-1,4:2:end,3) + L(3:2:end-1,4:2:end,6)));
P(2:end-1,2:end-1,6) = ...
    -(L(4:2:end,3:2:end-1,2) + 0.5*(L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,8))) ./ ...
    (L(4:2:end,3:2:end-1,1)  + 0.5*(L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,8) + ...
    L(4:2:end,3:2:end-1,4) + L(4:2:end,3:2:end-1,7)));
P(2:end-1,2:end-1,7) = ...
    -(L(4:2:end,2:2:end-2,3) + 0.5*(L(4:2:end,2:2:end-2,2) + L(4:2:end,2:2:end-2,4))) ./ ...
    (L(4:2:end,2:2:end-2,1)  + 0.5*(L(4:2:end,2:2:end-2,2) + L(4:2:end,2:2:end-2,4) + ...
    L(4:2:end,2:2:end-2,6) + L(4:2:end,2:2:end-2,8)));
P(2:end-1,2:end-1,8) = ...
    -(L(3:2:end-1,2:2:end-2,4) + 0.5*(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,6))) ./ ...
    (L(3:2:end-1,2:2:end-2,1)  + 0.5*(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,6) + ...
    L(3:2:end-1,2:2:end-2,2) + L(3:2:end-1,2:2:end-2,7)));
% Refrain from interpolating to pads when m and/or n are even.
if (mod(m,2) == 0)
    P(end-1,2:end-1,5:7) = 0;
end % if
if (mod(n,2) == 0)
    P(2:end-1,end-1,3:5) = 0;
end % if

return

function P = Triangle_Prol_one(m, n, m_c, n_c, L) % Prolongation from two near neighbors -
% suitable for the case of triangular grids
% The prolongation operator is a matrix of size  (m_c,n_c,9) plus pad.
% P(15,23,1), for example, represents the coefficient multiplying the corresponding
% coarse-grid variable (15,23), when interpolating to the fine-grid
% variable that coincides with it, while P(15,23,3) is the coefficient
% multiplying this variable when interpolating to the fine-grid variable
% that is shifted -1 in the i direction (up), and +1 in the j direction (right).
P = zeros(m_c+2,n_c+2,9);
P(2:end-1,2:end-1,1) = 1;
P(2:end-1,2:end-1,2) = ...
    -(L(2:2:end-2,3:2:end-1,6) + 1*(L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,7))) ./ ...
    (L(2:2:end-2,3:2:end-1,1)  + 0*(L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,7) + ...
    L(2:2:end-2,3:2:end-1,3) + L(2:2:end-2,3:2:end-1,8)));
P(2:end-1,2:end-1,3) = ...
    -(L(2:2:end-2,4:2:end,7) + 1*(L(2:2:end-2,4:2:end,6) + L(2:2:end-2,4:2:end,8))) ./ ...
    (L(2:2:end-2,4:2:end,1)  + 0*(L(2:2:end-2,4:2:end,6) + L(2:2:end-2,4:2:end,8) + ...
    L(2:2:end-2,4:2:end,2) + L(2:2:end-2,4:2:end,4)));
P(2:end-1,2:end-1,4) = ...
    -(L(3:2:end-1,4:2:end,8) + 1*(L(3:2:end-1,4:2:end,2) + L(3:2:end-1,4:2:end,7))) ./ ...
    (L(3:2:end-1,4:2:end,1)  + 0*(L(3:2:end-1,4:2:end,2) + L(3:2:end-1,4:2:end,7) + ...
    L(3:2:end-1,4:2:end,3) + L(3:2:end-1,4:2:end,6)));
P(2:end-1,2:end-1,6) = ...
    -(L(4:2:end,3:2:end-1,2) + 1*(L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,8))) ./ ...
    (L(4:2:end,3:2:end-1,1)  + 0*(L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,8) + ...
    L(4:2:end,3:2:end-1,4) + L(4:2:end,3:2:end-1,7)));
P(2:end-1,2:end-1,7) = ...
    -(L(4:2:end,2:2:end-2,3) + 1*(L(4:2:end,2:2:end-2,2) + L(4:2:end,2:2:end-2,4))) ./ ...
    (L(4:2:end,2:2:end-2,1)  + 0*(L(4:2:end,2:2:end-2,2) + L(4:2:end,2:2:end-2,4) + ...
    L(4:2:end,2:2:end-2,6) + L(4:2:end,2:2:end-2,8)));
P(2:end-1,2:end-1,8) = ...
    -(L(3:2:end-1,2:2:end-2,4) + 1*(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,6))) ./ ...
    (L(3:2:end-1,2:2:end-2,1)  + 0*(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,6) + ...
    L(3:2:end-1,2:2:end-2,2) + L(3:2:end-1,2:2:end-2,7)));

% Refrain from interpolating to pads when m and/or n are even.
if (mod(m,2) == 0)
    P(end-1,2:end-1,5:7) = 0;
end % if
if (mod(n,2) == 0)
    P(2:end-1,end-1,3:5) = 0;
end % if

return

function P = Triangle_Prol_zero(m, n, m_c, n_c, L) % Prolongation from two near neighbors -
% suitable for the case of triangular grids

% The prolongation operator is a matrix of size  (m_c,n_c,9) plus pad.
% P(15,23,1), for example, represents the coefficient multiplying the corresponding
% coarse-grid variable (15,23), when interpolating to the fine-grid
% variable that coincides with it, while P(15,23,3) is the coefficient
% multiplying this variable when interpolating to the fine-grid variable
% that is shifted -1 in the i direction (up), and +1 in the j direction (right).
P = zeros(m_c+2,n_c+2,9);
P(2:end-1,2:end-1,1) = 1;
P(2:end-1,2:end-1,2) = ...
    -(L(2:2:end-2,3:2:end-1,6) + 0*(L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,7))) ./ ...
    (L(2:2:end-2,3:2:end-1,1)  + 1*(L(2:2:end-2,3:2:end-1,4) + L(2:2:end-2,3:2:end-1,7) + ...
    L(2:2:end-2,3:2:end-1,3) + L(2:2:end-2,3:2:end-1,8)));
P(2:end-1,2:end-1,3) = ...
    -(L(2:2:end-2,4:2:end,7) + 0*(L(2:2:end-2,4:2:end,6) + L(2:2:end-2,4:2:end,8))) ./ ...
    (L(2:2:end-2,4:2:end,1)  + 1*(L(2:2:end-2,4:2:end,6) + L(2:2:end-2,4:2:end,8) + ...
    L(2:2:end-2,4:2:end,2) + L(2:2:end-2,4:2:end,4)));
P(2:end-1,2:end-1,4) = ...
    -(L(3:2:end-1,4:2:end,8) + 0*(L(3:2:end-1,4:2:end,2) + L(3:2:end-1,4:2:end,7))) ./ ...
    (L(3:2:end-1,4:2:end,1)  + 1*(L(3:2:end-1,4:2:end,2) + L(3:2:end-1,4:2:end,7) + ...
    L(3:2:end-1,4:2:end,3) + L(3:2:end-1,4:2:end,6)));
P(2:end-1,2:end-1,6) = ...
    -(L(4:2:end,3:2:end-1,2) + 0*(L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,8))) ./ ...
    (L(4:2:end,3:2:end-1,1)  + 1*(L(4:2:end,3:2:end-1,3) + L(4:2:end,3:2:end-1,8) + ...
    L(4:2:end,3:2:end-1,4) + L(4:2:end,3:2:end-1,7)));
P(2:end-1,2:end-1,7) = ...
    -(L(4:2:end,2:2:end-2,3) + 0*(L(4:2:end,2:2:end-2,2) + L(4:2:end,2:2:end-2,4))) ./ ...
    (L(4:2:end,2:2:end-2,1)  + 1*(L(4:2:end,2:2:end-2,2) + L(4:2:end,2:2:end-2,4) + ...
    L(4:2:end,2:2:end-2,6) + L(4:2:end,2:2:end-2,8)));
P(2:end-1,2:end-1,8) = ...
    -(L(3:2:end-1,2:2:end-2,4) + 0*(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,6))) ./ ...
    (L(3:2:end-1,2:2:end-2,1)  + 1*(L(3:2:end-1,2:2:end-2,3) + L(3:2:end-1,2:2:end-2,6) + ...
    L(3:2:end-1,2:2:end-2,2) + L(3:2:end-1,2:2:end-2,7)));


% Refrain from interpolating to pads when m and/or n are even.
if (mod(m,2) == 0)
    P(end-1,2:end-1,5:7) = 0;
end % if
if (mod(n,2) == 0)
    P(2:end-1,end-1,3:5) = 0;
end % if

return
function f = Restrict(P, r, m_c, n_c) % Restrict residuals

f = zeros(m_c+2,n_c+2);
f(2:end-1,2:end-1) = ...
    P(2:end-1,2:end-1,1) .* r(3:2:end-1,3:2:end-1) + ...
    P(2:end-1,2:end-1,2) .* r(2:2:end-2,3:2:end-1) + ...
    P(2:end-1,2:end-1,3) .* r(2:2:end-2,4:2:end) + ...
    P(2:end-1,2:end-1,4) .* r(3:2:end-1,4:2:end) + ...
    P(2:end-1,2:end-1,5) .* r(4:2:end,4:2:end) + ...
    P(2:end-1,2:end-1,6) .* r(4:2:end,3:2:end-1) + ...
    P(2:end-1,2:end-1,7) .* r(4:2:end,2:2:end-2) + ...
    P(2:end-1,2:end-1,8) .* r(3:2:end-1,2:2:end-2) + ...
    P(2:end-1,2:end-1,9) .* r(2:2:end-2,2:2:end-2);

return

function u = Prolong(P, u_c, m, n) % Interpolate u_c

u = zeros(m+2,n+2);

u(3:2:end-1,3:2:end-1) = P(2:end-1,2:end-1,1) .* u_c(2:end-1,2:end-1);
u(2:2:end-2,3:2:end-1) = u(2:2:end-2,3:2:end-1) + ...
    P(2:end-1,2:end-1,2) .* u_c(2:end-1,2:end-1);
u(2:2:end-2,4:2:end) = u(2:2:end-2,4:2:end) + ...
    P(2:end-1,2:end-1,3) .* u_c(2:end-1,2:end-1);
u(3:2:end-1,4:2:end) = u(3:2:end-1,4:2:end) + ...
    P(2:end-1,2:end-1,4) .* u_c(2:end-1,2:end-1);
u(4:2:end,4:2:end) = u(4:2:end,4:2:end) + ...
    P(2:end-1,2:end-1,5) .* u_c(2:end-1,2:end-1);
u(4:2:end,3:2:end-1) = u(4:2:end,3:2:end-1) + ...
    P(2:end-1,2:end-1,6) .* u_c(2:end-1,2:end-1);
u(4:2:end,2:2:end-2) = u(4:2:end,2:2:end-2) + ...
    P(2:end-1,2:end-1,7) .* u_c(2:end-1,2:end-1);
u(3:2:end-1,2:2:end-2) = u(3:2:end-1,2:2:end-2) + ...
    P(2:end-1,2:end-1,8) .* u_c(2:end-1,2:end-1);
u(2:2:end-2,2:2:end-2) = u(2:2:end-2,2:2:end-2) + ...
    P(2:end-1,2:end-1,9) .* u_c(2:end-1,2:end-1);

return

function v = Lu(L,u) % Compute L(u)
v = u;
v(2:end-1,2:end-1) = ...
    L(2:end-1,2:end-1,1) .* u(2:end-1,2:end-1) + ...
    L(2:end-1,2:end-1,2) .* u(1:end-2,2:end-1) + ...
    L(2:end-1,2:end-1,3) .* u(1:end-2,3:end) + ...
    L(2:end-1,2:end-1,4) .* u(2:end-1,3:end) + ...
    L(2:end-1,2:end-1,5) .* u(3:end,3:end) + ...
    L(2:end-1,2:end-1,6) .* u(3:end,2:end-1) + ...
    L(2:end-1,2:end-1,7) .* u(3:end,1:end-2) + ...
    L(2:end-1,2:end-1,8) .* u(2:end-1,1:end-2) + ...
    L(2:end-1,2:end-1,9) .* u(1:end-2,1:end-2);

return
