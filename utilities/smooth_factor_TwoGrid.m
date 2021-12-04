% This function is written to analysis the smooth factor of two grid with
% diff. relaxations.
%
% Inputs:
%       L: the given stencil.
%       N: the number of test grids
%       omega: the value of damping factor in Jacobi method
%       v: the number of relaxation
%       param: paras, 'isPDE' if it is true, rediscretization coarsing otherwise
%       Galerkin. PRmethod, a struct first one defines the porlongation,
%       second one for the restriction.

function [mu_opt,eig_set] = smooth_factor_TwoGrid(L_stencil,N,omega,param)

if isfield(param,'isPDE')
    isPDE = param.isPDE;
else
    isPDE = true;
end
if isfield(param,'PRmethod')
    PRmethod = param.PRmethod;
else
    PRmethod = {'Bilinear','Bilinear'};
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

if isfield(param,'RelaxationMethod')
    RelaxationMethod = param.RelaxationMethod;
else
    RelaxationMethod = 'Jacobi';
end
epsilon = 1e-8;
matr_size = 4;

theta = -pi/2:2*pi/N:pi/2;
theta = theta(:);

i = sqrt(-1); % the complex unit
m = size(L_stencil,1);
m_ = (m-1)/2; % take the central position

L_hat_left = @(x)(exp(i*(-m_:1:m_)*x));
L_hat_right = @(x)(exp(i*(m_:-1:-m_).'*x));

L_hat = @(x1,x2)(reshape(kron(L_hat_left(x1),L_hat_right(x2)),m^2,1));

smooth_operat = @(L,x1,x2)(L(:).'*L_hat(x1,x2));
mu_set = zeros(N,N);
eig_set = [];
critera_low = @(x)(x>=-pi/2&&x<pi/2);
for index_i = 1:length(theta)
    theta_1 = theta(index_i);
    for index_j = 1:length(theta)
        theta_2 = theta(index_j);
        if ~critera_low(theta_1)
            continue;
        end
        
        [S_theta_1,S_theta_2,S_theta_3,S_theta_4,...
            L_h_1,L_h_2,L_h_3,L_h_4] = Smooth_L_h_calculate(L_stencil,theta_1,theta_2,omega,RelaxationMethod,smooth_operat);
        L_h = spdiags([L_h_1;L_h_2;L_h_3;L_h_4],0,matr_size,matr_size);
        if abs(L_h_1)<epsilon || abs(L_h_2)<epsilon || ...
                abs(L_h_3)<epsilon || abs(L_h_4)<epsilon
            %fprintf('Exclude set\n')
            continue;
        end
        if strcmp(RelaxationMethod,'RB')
            S_R = 0.5*blkdiag([S_theta_1+1 S_theta_2-1;S_theta_1-1 S_theta_2+1],[S_theta_3+1 S_theta_4-1;S_theta_3-1 S_theta_4+1]);
            S_B = 0.5*blkdiag([S_theta_1+1 -S_theta_2+1;-S_theta_1+1 S_theta_2+1],[S_theta_3+1 -S_theta_4+1;-S_theta_3+1 S_theta_4+1]);
            S_mat = S_B*S_R;
        else
            S_mat = spdiags([S_theta_1;S_theta_2;S_theta_3;S_theta_4],0,matr_size,matr_size);
        end
        [I_h_H_1,I_h_H_2,I_h_H_3,I_h_H_4] = SymbolTransfer(theta_1,theta_2,PRmethod{1});
        P_mat = [I_h_H_1 I_h_H_2 I_h_H_3 I_h_H_4]';
        if strcmp(PRmethod{1},PRmethod{2})
            R_mat = P_mat';
        else
            [I_h_H_1,I_h_H_2,I_h_H_3,I_h_H_4] = SymbolTransfer(theta_1,theta_2,PRmethod{2});
            R_mat = [I_h_H_1 I_h_H_2 I_h_H_3 I_h_H_4];
        end
        
        if isPDE
            L_H = 0.25*smooth_operat(L_stencil,2*theta_1,2*theta_2);
            % build coarse correction
            if abs(L_H)>epsilon
                Q_h = (speye(matr_size)-P_mat*(1/L_H)*R_mat*L_h);
                %Q_h = spdiags([0;1;1;1],0,matr_size,matr_size);
            else
                continue;
            end
        else
            L_H = R_mat*L_h*P_mat;
            % coarse correction
            Q_h = (speye(matr_size)-P_mat*(1/L_H)*R_mat*...
                L_h);
        end
        
        E_h = S_mat^v2*Q_h*S_mat^v1;
        if (v1+v2)>1
            mu_set(index_i,index_j) = vrho(E_h);% nthroot(vrho(E_h),v1+v2);
        else
            mu_set(index_i,index_j) = vrho(E_h);
        end
        eig_set = [eig_set;eig(E_h)];%real(E_h)
    end
end
mu_opt = max(mu_set(:));
fprintf('The optimal smooth factor for true TG analysis is %f  \n',mu_opt)
clear smooth_operat
return

function [S_theta_1,S_theta_2,S_theta_3,S_theta_4,...
    L_h_1,L_h_2,L_h_3,L_h_4] = Smooth_L_h_calculate(L_stencil,theta_1,theta_2,omega,RelaxationMethod,smooth_operat)
% assume the size of L_stencil is square.
m = size(L_stencil,1);

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
    error('It has error\n');
end

L_h_1 = smooth_operat(L_stencil,theta_1,theta_2);
L_h_2 = smooth_operat(L_stencil,theta_1_hat,theta_2_hat);
L_h_3 = smooth_operat(L_stencil,theta_1_hat,theta_2);
L_h_4 = smooth_operat(L_stencil,theta_1,theta_2_hat);

switch RelaxationMethod
    case 'Jacobi'
        L_pos = zeros(m,m);
        m_ = (m-1)/2;
        L_pos(m_+1,m_+1) = L_stencil(m_+1,m_+1);
        S_theta_1 = 1-omega*smooth_operat(L_stencil,theta_1,theta_2)/smooth_operat(L_pos,theta_1,theta_2);
        S_theta_2 = 1-omega*smooth_operat(L_stencil,theta_1_hat,theta_2_hat)/smooth_operat(L_pos,theta_1_hat,theta_2_hat);
        S_theta_3 = 1-omega*smooth_operat(L_stencil,theta_1_hat,theta_2)/smooth_operat(L_pos,theta_1_hat,theta_2);
        S_theta_4 = 1-omega*smooth_operat(L_stencil,theta_1,theta_2_hat)/smooth_operat(L_pos,theta_1,theta_2_hat);
    case 'RB'
        L_pos = zeros(m,m);
        m_ = (m-1)/2;
        L_pos(m_+1,m_+1) = L_stencil(m_+1,m_+1);
        S_theta_1 = 1-omega*smooth_operat(L_stencil,theta_1,theta_2)/smooth_operat(L_pos,theta_1,theta_2);
        S_theta_2 = 1-omega*smooth_operat(L_stencil,theta_1_hat,theta_2_hat)/smooth_operat(L_pos,theta_1_hat,theta_2_hat);
        S_theta_3 = 1-omega*smooth_operat(L_stencil,theta_1_hat,theta_2)/smooth_operat(L_pos,theta_1_hat,theta_2);
        S_theta_4 = 1-omega*smooth_operat(L_stencil,theta_1,theta_2_hat)/smooth_operat(L_pos,theta_1,theta_2_hat);
    case 'GS'
        L_low = zeros(m,m);
        m_ = (m-1)/2;
        L_low(1,:) = L_stencil(1,:);
        L_low(2,1:m_+1) = L_stencil(2,1:m_+1);
        S_theta_1 = 1-omega*smooth_operat(L_stencil,theta_1,theta_2)/smooth_operat(L_low,theta_1,theta_2);
        S_theta_2 = 1-omega*smooth_operat(L_stencil,theta_1_hat,theta_2_hat)/smooth_operat(L_low,theta_1_hat,theta_2_hat);
        S_theta_3 = 1-omega*smooth_operat(L_stencil,theta_1_hat,theta_2)/smooth_operat(L_low,theta_1_hat,theta_2);
        S_theta_4 = 1-omega*smooth_operat(L_stencil,theta_1,theta_2_hat)/smooth_operat(L_low,theta_1,theta_2_hat);
    case 'SymGS'
        L_low = zeros(m,m);
        L_diag = zeros(m,m);
        m_ = (m-1)/2;
        L_diag(m_+1,m_+1) = L_stencil(m_+1,m_+1);
        L_low(1,:) = L_stencil(1,:);
        L_low(2,1:m_+1) = L_stencil(2,1:m_+1);
        L_pos = L_stencil-L_low+L_diag;
        S_theta_1 = 1-omega*smooth_operat(L_diag,theta_1,theta_2)*smooth_operat(L_stencil,theta_1,theta_2)/smooth_operat(L_low,theta_1,theta_2)/smooth_operat(L_pos,theta_1,theta_2);
        S_theta_2 = 1-omega*smooth_operat(L_diag,theta_1_hat,theta_2_hat)*smooth_operat(L_stencil,theta_1_hat,theta_2_hat)/smooth_operat(L_low,theta_1_hat,theta_2_hat)/smooth_operat(L_pos,theta_1_hat,theta_2_hat);
        S_theta_3 = 1-omega*smooth_operat(L_diag,theta_1_hat,theta_2)*smooth_operat(L_stencil,theta_1_hat,theta_2)/smooth_operat(L_low,theta_1_hat,theta_2)/smooth_operat(L_pos,theta_1_hat,theta_2);
        S_theta_4 = 1-omega*smooth_operat(L_diag,theta_1,theta_2_hat)*smooth_operat(L_stencil,theta_1,theta_2_hat)/smooth_operat(L_low,theta_1,theta_2_hat)/smooth_operat(L_pos,theta_1,theta_2_hat);
    otherwise
        error([LocalSmooth ' Smoother is not defined\n'])
end
return
