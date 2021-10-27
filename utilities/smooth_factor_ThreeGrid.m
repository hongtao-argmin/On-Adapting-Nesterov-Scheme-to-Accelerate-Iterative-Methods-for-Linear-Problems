% This function is written to analysis the smooth factor of three-grid with
% Jacobi relaxation.
%
% L: the stencil
% N: the number of points
% v1: the number of pre-relaxation. The default is one.
% v2: the number of post-relaxation. The default is zero.
% param: a param to control different formulations of coarse and the coarse
% correction
% ==========================================================================

function [mu_opt,eig_set] = smooth_factor_ThreeGrid(L,N,...
    param,isShow)

if nargin<=3
    isShow = true;
end

if isfield(param,'CycType')
    CycType = param.CycType;
else
    CycType = 1;
end

if isfield(param,'Omega_opt') % the optimally damped for the Jacobi method.
    Omega_opt = param.Omega_opt;
else
    Omega_opt = 1;
end

if isfield(param,'v1')
    v1 = param.v1;
else
    v1 = 1;
end

if isfield(param,'v2')
    v2 = param.v2;
else
    v2 = 0;
end

% isPi_2H decides whether contains the boundary part
if isfield(param,'isPi_2H')
    isPi_2H = param.isPi_2H;
else
    isPi_2H = false;
end

if isfield(param,'RestrictionMethod')
    RestrictionMethod = param.RestrictionMethod;
else
    RestrictionMethod = {'Bilinear','Bilinear'};
end
if isfield(param,'RelaxationMethod')
    RelaxationMethod = param.RelaxationMethod;
else
    RelaxationMethod = 'Jacobi';
end

epsilon = 1e-8; % tolerance to avoid numerical error
matr_size = 4;  % 4\times 4 matrix - 2D
i = sqrt(-1);   % the complex unit
m = size(L,1);

% if sum(sum(abs(L-L')))<epsilon
%     isSymmtric = true;
% else
%     isSymmtric = false;
% end

% the size of stencil is 3*3
% formulate the matrix for calculating the symbol
global smooth_operator smooth_operat_Jacobi
L_hat = @(angle_1,angle_2)(reshape([exp(-i*angle_1-i*angle_2) exp(-i*angle_1) exp(-i*angle_1+i*angle_2);...
    exp(-i*angle_2) 1 exp(i*angle_2);exp(i*angle_1-i*angle_2) exp(i*angle_1) exp(i*angle_1+i*angle_2)],9,1));
smooth_operator = @(L_kernel,angle_1,angle_2)(L_kernel(:)'*L_hat(angle_1,angle_2));

switch RelaxationMethod
    case 'Jacobi'
        smooth_operat_Jacobi = @(L_kernel,angle_1,angle_2)(1-...
            Omega_opt*((L_kernel(:)'*L_hat(angle_1,angle_2))/L_kernel(2,2)));
    case 'RB'
        smooth_operat_Jacobi = @(L_kernel,angle_1,angle_2)(1-...
            Omega_opt*((L_kernel(:)'*L_hat(angle_1,angle_2))/L_kernel(2,2)));
    case 'GS'
        L_low = zeros(m,m);
        m_ = (m-1)/2;
        L_low(1,:) = L(1,:);
        L_low(2,1:m_+1) = L(2,1:m_+1);
        smooth_operat_Jacobi = @(L_kernel,angle_1,angle_2)(1-Omega_opt*smooth_operator(L_kernel,angle_1,angle_2)/smooth_operator(L_low,angle_1,angle_2));
    case 'SymGS'
        L_low = zeros(m,m);
        L_diag = zeros(m,m);
        m_ = (m-1)/2;
        L_diag(m_+1,m_+1) = L(m_+1,m_+1);
        L_low(1,:) = L(1,:);
        L_low(2,1:m_+1) = L(2,1:m_+1);
        L_pos = L-L_low+L_diag;
        smooth_operat_Jacobi = @(L_kernel,angle_1,angle_2)(1-Omega_opt*smooth_operator(L_diag,angle_1,angle_2)*smooth_operator(L_kernel,angle_1,angle_2)/smooth_operator(L_low,angle_1,angle_2)/smooth_operator(L_pos,angle_1,angle_2));
    otherwise
        error([LocalSmooth ' Smoother is not defined\n'])
end
% focus on low-frequency so the number of samples reduces to N/4 since this
% is three-grid analysis.
theta = -pi/4:2*pi/N:pi/4;
theta = theta(:);

mu_set = zeros(N/4+1,N/4+1);
if nargout>1
    eig_set = [];
end
for index_i = 1:length(theta)
    theta_1 = theta(index_i);
    if isPi_2H
        if abs(theta_1)>=pi/4
            continue;
        end
    else
        if abs(theta_1)>pi/4
            continue;
        end
    end
    
    for index_j = 1:length(theta)
        theta_2 = theta(index_j);
        if isPi_2H
            if  abs(theta_2)>=pi/4
                continue;
            end
        else
            if  abs(theta_2)>pi/4
                continue;
            end
        end
        
        %         if isSymmtric&&theta_2<0
        %             continue;
        %         end
        
        % build analysis on pi/2
        [theta_1_hat,theta_2_hat] = AngleAlias(theta_1,theta_2,'Pi_2');
        
        % pay attention to the 2\theta and 2h
        theta_1_2g = 2*theta_1;
        theta_2_2g = 2*theta_2;
        theta_1_hat_2g = 2*theta_1_hat;
        theta_2_hat_2g = 2*theta_2_hat;
        [L_2g_h_1,L_2g_h_2,L_2g_h_3,L_2g_h_4] = Smooth_L_h_Pi_2(0.25*L,theta_1_2g,theta_2_2g,theta_1_hat_2g,theta_2_hat_2g);
        
        if abs(L_2g_h_1)<epsilon || abs(L_2g_h_2)<epsilon || ...
                abs(L_2g_h_3)<epsilon || abs(L_2g_h_4)<epsilon
            continue;
        end
        
        % formulate the matrix on the second level first
        I_2g = speye(matr_size);
        
        [R_2gMat,P_2gMat] = FormualteRP(theta_1_2g,theta_2_2g,theta_1_hat_2g,theta_2_hat_2g,RestrictionMethod);
        
        % formulate smooth matrix
        if strcmp(RelaxationMethod,'RB')
            [S_2g_R,S_2g_B] = RBSmoothOperator(L,theta_1_2g,theta_2_2g,theta_1_hat_2g,theta_2_hat_2g);
            S_2g = S_2g_B*S_2g_R;
        else
            S_2g = spdiags([smooth_operat_Jacobi(L,theta_1_2g,theta_2_2g);...
                smooth_operat_Jacobi(L,theta_1_hat_2g,theta_2_hat_2g);...
                smooth_operat_Jacobi(L,theta_1_hat_2g,theta_2_2g);...
                smooth_operat_Jacobi(L,theta_1_2g,theta_2_hat_2g)],0,matr_size,matr_size);
        end
        L_2g_h = spdiags([L_2g_h_1;L_2g_h_2;L_2g_h_3;L_2g_h_4],0,matr_size,matr_size);
        L_2g_H = smooth_operator((0.25^2)*L,2*theta_1_2g,2*theta_2_2g);
        if abs(L_2g_H)<epsilon
            continue;
        end
        
        M_2g_h_H = (S_2g^v2)*(I_2g-P_2gMat*(1/L_2g_H)*R_2gMat*L_2g_h)*(S_2g^v1);
        L_3g_H_inv = (I_2g-M_2g_h_H^CycType)/L_2g_h; % main part of the iteration matrix on the second level
        
        %         formulate the matrix on the finest level first
        
        [theta_1_hat_3g,theta_2_hat_3g] = AngleAlias(theta_1,theta_2,'Pi');
        [L_h_1_temp,L_h_2_temp,L_h_3_temp,L_h_4_temp] = Smooth_L_h_Pi(L,theta_1,theta_2,theta_1_hat_3g,theta_2_hat_3g);
        if abs(L_h_1_temp)<epsilon || abs(L_h_2_temp)<epsilon || ...
                abs(L_h_3_temp)<epsilon || abs(L_h_2_temp)<epsilon
            continue;
        end
        [R_2gMat_temp_1,P_2gMat_temp_1] = FormualteRP(theta_1,theta_2,theta_1_hat_3g,theta_2_hat_3g,RestrictionMethod);
        L_h_3g_temp1 = spdiags([L_h_1_temp;L_h_2_temp;L_h_3_temp;L_h_4_temp],0,matr_size,matr_size);
        if strcmp(RelaxationMethod,'RB')
            [S_3g_R,S_3g_B] = RBSmoothOperator(L,theta_1,theta_2,theta_1_hat_3g,theta_2_hat_3g);
            S_3g_temp1 = S_3g_B*S_3g_R;
        else
            S_3g_temp1 = spdiags([smooth_operat_Jacobi(L,theta_1,theta_2);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2_hat_3g);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2);...
                smooth_operat_Jacobi(L,theta_1,theta_2_hat_3g)],0,matr_size,matr_size);
        end
        
        [theta_1_hat_3g,theta_2_hat_3g] = AngleAlias(theta_1_hat,theta_2_hat,'Pi');
        [L_h_1_temp,L_h_2_temp,L_h_3_temp,L_h_4_temp] = Smooth_L_h_Pi(L,theta_1_hat,theta_2_hat,theta_1_hat_3g,theta_2_hat_3g);
        if abs(L_h_1_temp)<epsilon || abs(L_h_2_temp)<epsilon || ...
                abs(L_h_3_temp)<epsilon || abs(L_h_2_temp)<epsilon
            continue;
        end
        [R_2gMat_temp_2,P_2gMat_temp_2] = FormualteRP(theta_1_hat,theta_2_hat,theta_1_hat_3g,theta_2_hat_3g,RestrictionMethod);
        L_h_3g_temp2 = spdiags([L_h_1_temp;L_h_2_temp;L_h_3_temp;L_h_4_temp],0,matr_size,matr_size);
        if strcmp(RelaxationMethod,'RB')
            [S_3g_R,S_3g_B] = RBSmoothOperator(L,theta_1_hat,theta_2_hat,theta_1_hat_3g,theta_2_hat_3g);
            S_3g_temp2 = S_3g_B*S_3g_R;
        else
            S_3g_temp2 = spdiags([smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2_hat_3g);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2_hat);...
                smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat_3g)],0,matr_size,matr_size);
        end
        
        
        [theta_1_hat_3g,theta_2_hat_3g] = AngleAlias(theta_1_hat,theta_2,'Pi');
        [L_h_1_temp,L_h_2_temp,L_h_3_temp,L_h_4_temp] = Smooth_L_h_Pi(L,theta_1_hat,theta_2,theta_1_hat_3g,theta_2_hat_3g);
        if abs(L_h_1_temp)<epsilon || abs(L_h_2_temp)<epsilon || ...
                abs(L_h_3_temp)<epsilon || abs(L_h_2_temp)<epsilon
            continue;
        end
        [R_2gMat_temp_3,P_2gMat_temp_3] = FormualteRP(theta_1_hat,theta_2,theta_1_hat_3g,theta_2_hat_3g,RestrictionMethod);
        L_h_3g_temp3 = spdiags([L_h_1_temp;L_h_2_temp;L_h_3_temp;L_h_4_temp],0,matr_size,matr_size);
        if strcmp(RelaxationMethod,'RB')
            [S_3g_R,S_3g_B] = RBSmoothOperator(L,theta_1_hat,theta_2,theta_1_hat_3g,theta_2_hat_3g);
            S_3g_temp3 = S_3g_B*S_3g_R;
        else
            S_3g_temp3 = spdiags([smooth_operat_Jacobi(L,theta_1_hat,theta_2);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2_hat_3g);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2);...
                smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat_3g)],0,matr_size,matr_size);
        end
        
        [theta_1_hat_3g,theta_2_hat_3g] = AngleAlias(theta_1,theta_2_hat,'Pi');
        [L_h_1_temp,L_h_2_temp,L_h_3_temp,L_h_4_temp] = Smooth_L_h_Pi(L,theta_1,theta_2_hat,theta_1_hat_3g,theta_2_hat_3g);
        if abs(L_h_1_temp)<epsilon || abs(L_h_2_temp)<epsilon || ...
                abs(L_h_3_temp)<epsilon || abs(L_h_2_temp)<epsilon
            continue;
        end
        [R_2gMat_temp_4,P_2gMat_temp_4] = FormualteRP(theta_1,theta_2_hat,theta_1_hat_3g,theta_2_hat_3g,RestrictionMethod);
        L_h_3g_temp4 = spdiags([L_h_1_temp;L_h_2_temp;L_h_3_temp;L_h_4_temp],0,matr_size,matr_size);
        
        if strcmp(RelaxationMethod,'RB')
            [S_3g_R,S_3g_B] = RBSmoothOperator(L,theta_1,theta_2_hat,theta_1_hat_3g,theta_2_hat_3g);
            S_3g_temp4 = S_3g_B*S_3g_R;
        else
            S_3g_temp4 = spdiags([smooth_operat_Jacobi(L,theta_1,theta_2_hat);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2_hat_3g);...
                smooth_operat_Jacobi(L,theta_1_hat_3g,theta_2_hat);...
                smooth_operat_Jacobi(L,theta_1,theta_2_hat_3g)],0,matr_size,matr_size);
        end
        
        L_3g_h = sparse(blkdiag(L_h_3g_temp1,L_h_3g_temp2,L_h_3g_temp3,L_h_3g_temp4));
        R_3gMat = sparse(blkdiag(R_2gMat_temp_1,R_2gMat_temp_2,R_2gMat_temp_3,R_2gMat_temp_4));
        P_3gMat = sparse(blkdiag(P_2gMat_temp_1,P_2gMat_temp_2,P_2gMat_temp_3,P_2gMat_temp_4));
        S_3gMat = sparse(blkdiag(S_3g_temp1,S_3g_temp2,S_3g_temp3,S_3g_temp4));
        
        E_h = S_3gMat^v2*(speye(2^4)-P_3gMat*L_3g_H_inv*R_3gMat*L_3g_h)*S_3gMat^v1;
        mu_set(index_i,index_j) = vrho(E_h);%nthroot(vrho(E_h),v1+v2);
        if nargout>1
            eig_set = [eig_set;eig(E_h)]; % real(E_h)
        end
    end
end
mu_opt = max(mu_set(:));

if isShow
    fprintf('The optimal smooth factor is %f \n',mu_opt)
end
clear  global
return
function [S_R,S_B] = RBSmoothOperator(L,theta_1,theta_2,theta_1_hat,theta_2_hat)
global smooth_operat_Jacobi

S_R = 0.5*blkdiag([smooth_operat_Jacobi(L,theta_1,theta_2)+1 smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat)-1;...
    smooth_operat_Jacobi(L,theta_1,theta_2)-1 smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat)+1],[smooth_operat_Jacobi(L,theta_1,theta_2_hat)+1 smooth_operat_Jacobi(L,theta_1_hat,theta_2)-1;...
    smooth_operat_Jacobi(L,theta_1,theta_2_hat)-1 smooth_operat_Jacobi(L,theta_1_hat,theta_2)+1]);
S_B = 0.5*blkdiag([smooth_operat_Jacobi(L,theta_1,theta_2)+1 -smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat)+1;...
    -smooth_operat_Jacobi(L,theta_1,theta_2)+1 smooth_operat_Jacobi(L,theta_1_hat,theta_2_hat)+1],[smooth_operat_Jacobi(L,theta_1,theta_2_hat)+1 -smooth_operat_Jacobi(L,theta_1_hat,theta_2)+1;...
    -smooth_operat_Jacobi(L,theta_1,theta_2_hat)+1 smooth_operat_Jacobi(L,theta_1_hat,theta_2)+1]);


return
function [L_h_1,L_h_2,L_h_3,L_h_4] = Smooth_L_h_Pi_2(L,theta_1,theta_2,theta_1_hat,theta_2_hat)
global  smooth_operator
if nargin<4
    Type = 'Pi_2';
    [theta_1_hat,theta_2_hat] = AngleAlias(theta_1,theta_2,Type);
end

L_h_1 = smooth_operator(L,theta_1,theta_2);
L_h_2 = smooth_operator(L,theta_1_hat,theta_2_hat);
L_h_3 = smooth_operator(L,theta_1_hat,theta_2);
L_h_4 = smooth_operator(L,theta_1,theta_2_hat);

return

function [L_h_1,L_h_2,L_h_3,L_h_4,L_h_deno1,L_h_deno2] = Smooth_L_h_Pi(L,theta_1,theta_2,theta_1_hat,theta_2_hat)
global  smooth_operator

if nargin<4
    Type = 'Pi';
    [theta_1_hat,theta_2_hat] = AngleAlias(theta_1,theta_2,Type);
end

L_h_1 = smooth_operator(L,theta_1,theta_2);
L_h_2 = smooth_operator(L,theta_1_hat,theta_2_hat);
L_h_3 = smooth_operator(L,theta_1_hat,theta_2);
L_h_4 = smooth_operator(L,theta_1,theta_2_hat);

if nargout>4
    % -------------------------------------------------------------------------
    L_temp_1 = triu(L) - diag([L(1,1);0;0]);
    L_temp_1(1,2) = 0;
    L_temp_1(3,2) =  L(3,2);
    % L_temp_1_num = L - L_temp_1;
    
    L_h_1_deno1 = smooth_operator(L_temp_1,theta_1,theta_2);
    % L_h_1_num1 =  smooth_operator(L_temp_1_num,theta_1,theta_2);
    
    L_h_2_deno1 = smooth_operator(L_temp_1,theta_1_hat,theta_2_hat);
    % L_h_2_num1 = smooth_operator(L_temp_1_num,theta_1_hat,theta_2_hat);
    
    L_h_3_deno1 = smooth_operator(L_temp_1,theta_1_hat,theta_2);
    % L_h_3_num1 = smooth_operator(L_temp_1_num,theta_1_hat,theta_2);
    
    L_h_4_deno1 = smooth_operator(L_temp_1,theta_1,theta_2_hat);
    % L_h_4_num1 = smooth_operator(L_temp_1_num,theta_1,theta_2_hat);
    % L_h_backward = [L_h_1_num1;L_h_2_num1;L_h_3_num1;L_h_4_num1]./[L_h_1_deno1;L_h_2_deno1;L_h_3_deno1;L_h_4_deno1];
    % L_h_backward = diag(L_h_backward);
    L_h_deno1 = [L_h_1_deno1;L_h_2_deno1;L_h_3_deno1;L_h_4_deno1];
    L_h_deno1 = 1./L_h_deno1;
    L_h_deno1 = diag(L_h_deno1);
    % -------------------------------------------------------------------------
    
    L_temp_2 = tril(L) - diag([0;0;L(3,3)]);
    L_temp_2(3,2) = 0;
    L_temp_2(1,2) = L(1,2);
    % L_temp_2_num = L - L_temp_2;
    
    L_h_1_deno2 = smooth_operator(L_temp_2,theta_1,theta_2);
    % L_h_1_num2 = smooth_operator(L_temp_2_num,theta_1,theta_2);
    
    L_h_2_deno2 = smooth_operator(L_temp_2,theta_1_hat,theta_2_hat);
    % L_h_2_num2 = smooth_operator(L_temp_2_num,theta_1_hat,theta_2_hat);
    
    L_h_3_deno2 = smooth_operator(L_temp_2,theta_1_hat,theta_2);
    % L_h_3_num2 = smooth_operator(L_temp_2_num,theta_1_hat,theta_2);
    
    L_h_4_deno2 = smooth_operator(L_temp_2,theta_1,theta_2_hat);
    % L_h_4_num2 = smooth_operator(L_temp_2_num,theta_1,theta_2_hat);
    % L_h_forward = [L_h_1_num2;L_h_2_num2;L_h_3_num2;L_h_4_num2]./[L_h_1_deno2;L_h_2_deno2;L_h_3_deno2;L_h_4_deno2];
    % L_h_forward = diag(L_h_forward);
    L_h_deno2 = [L_h_1_deno2;L_h_2_deno2;L_h_3_deno2;L_h_4_deno2];
    L_h_deno2 = 1./L_h_deno2;
    L_h_deno2 = diag(L_h_deno2);
end
% -------------------------------------------------------------------------
return

function [theta_1_hat,theta_2_hat] = AngleAlias(theta_1,theta_2,Type)

if strcmp(Type,'Pi_2')
    if theta_1>=0&&theta_2>=0
        theta_1_hat = theta_1-pi/2;
        theta_2_hat = theta_2-pi/2;
    elseif theta_1>=0&&theta_2<0
        theta_1_hat = theta_1-pi/2;
        theta_2_hat = theta_2+pi/2;
    elseif theta_1<0&&theta_2>=0
        theta_1_hat = theta_1+pi/2;
        theta_2_hat = theta_2-pi/2;
    elseif theta_1<0&&theta_2<0
        theta_1_hat = theta_1+pi/2;
        theta_2_hat = theta_2+pi/2;
    else
        error('It has an error\n');
    end
elseif strcmp(Type,'Pi')
    if theta_1>=0&&theta_2>=0
        theta_1_hat = theta_1-pi;
        theta_2_hat = theta_2-pi;
    elseif theta_1>=0&&theta_2<0
        theta_1_hat = theta_1-pi;
        theta_2_hat = theta_2+pi;
    elseif theta_1<0&&theta_2>=0
        theta_1_hat = theta_1+pi;
        theta_2_hat = theta_2-pi;
    elseif theta_1<0&&theta_2<0
        theta_1_hat = theta_1+pi;
        theta_2_hat = theta_2+pi;
    else
        error('It has an error\n');
    end
else
    error('We do not have such an algle transform')
end

return

function [R_mat,P_mat] = FormualteRP(theta_1,theta_2,theta_1_hat,theta_2_hat,RestrictionMethod)

switch RestrictionMethod{1}
    case 'Bilinear'
        I_h_H_1 = 0.25*(1+cos(theta_1))*(1+cos(theta_2));
        I_h_H_2 = 0.25*(1+cos(theta_1_hat))*(1+cos(theta_2_hat));
        I_h_H_3 = 0.25*(1+cos(theta_1_hat))*(1+cos(theta_2));
        I_h_H_4 = 0.25*(1+cos(theta_1))*(1+cos(theta_2_hat));
    case 'Cubic'
        I_h_H_1 = 0.25*(1+9/8*cos(theta_1)-0.125*cos(3*theta_1))*(1+1.125*cos(theta_2)-0.125*cos(3*theta_2));
        I_h_H_2 = 0.25*(1+1.125*cos(theta_1_hat)-0.125*cos(3*theta_1_hat))*(1+1.125*cos(theta_2_hat)-0.125*cos(3*theta_2_hat));
        I_h_H_3 = 0.25*(1+1.125*cos(theta_1_hat)-1/8*cos(3*theta_1_hat))*(1+1.125*cos(theta_2)-0.125*cos(3*theta_2));
        I_h_H_4 = 0.25*(1+1.125*cos(theta_1)-0.125*cos(3*theta_1))*(1+1.125*cos(theta_2_hat)-0.125*cos(3*theta_2_hat));
    case '6th'
        a = 2*150/256;
        b = -2*25/256;
        c = 6/256;
        I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2));
        I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat));
        I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2));
        I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat));
    case '8th'
        a = 2*1225/2048;
        b = -2*245/2048;
        c = 2*49/2048;
        d = -10/2048;
        I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2));
        I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat));
        I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2));
        I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat));
    case '10th'
        x = [1 1 1 1 1;1 3^2 5^2 7^2 9^2;1 3^4 5^4 7^4 9^4;1 3^6 5^6 7^6 9^6;1 3^8 5^8 7^8 9^8]\[0.5;0;0;0;0];
        x = 2*x;
        a = x(1);
        b = x(2);
        c = x(3);
        d = x(4);
        e = x(5);
        I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1)+e*cos(9*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2)+e*cos(9*theta_2));
        I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat)+e*cos(9*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat)+e*cos(9*theta_2_hat));
        I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat)+e*cos(9*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2)+e*cos(9*theta_2));
        I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1)+e*cos(9*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat)+e*cos(9*theta_2_hat));
    otherwise
        error([RestrictionMethod ' is not defined \n'])
end

P_mat = [I_h_H_1 I_h_H_2 I_h_H_3 I_h_H_4]';

if strcmp(RestrictionMethod{2},RestrictionMethod{1})
    R_mat = P_mat';
else
    switch RestrictionMethod{2}
        case 'Bilinear'
            I_h_H_1 = 0.25*(1+cos(theta_1))*(1+cos(theta_2));
            I_h_H_2 = 0.25*(1+cos(theta_1_hat))*(1+cos(theta_2_hat));
            I_h_H_3 = 0.25*(1+cos(theta_1_hat))*(1+cos(theta_2));
            I_h_H_4 = 0.25*(1+cos(theta_1))*(1+cos(theta_2_hat));
        case 'Cubic'
            I_h_H_1 = 0.25*(1+9/8*cos(theta_1)-0.125*cos(3*theta_1))*(1+1.125*cos(theta_2)-0.125*cos(3*theta_2));
            I_h_H_2 = 0.25*(1+1.125*cos(theta_1_hat)-0.125*cos(3*theta_1_hat))*(1+1.125*cos(theta_2_hat)-0.125*cos(3*theta_2_hat));
            I_h_H_3 = 0.25*(1+1.125*cos(theta_1_hat)-1/8*cos(3*theta_1_hat))*(1+1.125*cos(theta_2)-0.125*cos(3*theta_2));
            I_h_H_4 = 0.25*(1+1.125*cos(theta_1)-0.125*cos(3*theta_1))*(1+1.125*cos(theta_2_hat)-0.125*cos(3*theta_2_hat));
        case '6th'
            a = 2*150/256;
            b = -2*25/256;
            c = 6/256;
            I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2));
            I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat));
            I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2));
            I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat));
        case '8th'
            a = 2*1225/2048;
            b = -2*245/2048;
            c = 2*49/2048;
            d = -10/2048;
            I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2));
            I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat));
            I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2));
            I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat));
        case '10th'
            x = [1 1 1 1 1;1 3^2 5^2 7^2 9^2;1 3^4 5^4 7^4 9^4;1 3^6 5^6 7^6 9^6;1 3^8 5^8 7^8 9^8]\[0.5;0;0;0;0];
            x = 2*x;
            a = x(1);
            b = x(2);
            c = x(3);
            d = x(4);
            e = x(5);
            I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1)+e*cos(9*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2)+e*cos(9*theta_2));
            I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat)+e*cos(9*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat)+e*cos(9*theta_2_hat));
            I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat)+e*cos(9*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2)+e*cos(9*theta_2));
            I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1)+e*cos(9*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat)+e*cos(9*theta_2_hat));
        otherwise
            error([RestrictionMethod ' is not defined \n'])
    end
    R_mat = [I_h_H_1 I_h_H_2 I_h_H_3 I_h_H_4];
end

return
