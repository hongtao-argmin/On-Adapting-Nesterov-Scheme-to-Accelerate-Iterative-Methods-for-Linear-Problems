% Test on Laplacian problem.
%%
% clc;
clear;
close all;
%%
addpath('./utilities')
%%
N = 1024;
nu1 = 1;
nu2 = 0;
n_cycles = 50;
epsi = 1; 
phi = 0;
CycType = 1;
N_setup = 64;
StopTolerance = 1e-8;
Magfy = 1;
%%
C = cos(phi);
S = sin(phi);
L_stencil = -[-0.5*(1-epsi)*C*S epsi*C^2+S^2 0.5*(1-epsi)*C*S;...
    C^2+epsi*S^2 -2*(1+epsi) C^2+epsi*S^2; ...
    0.5*(1-epsi)*C*S epsi*C^2+S^2 -0.5*(1-epsi)*C*S];
omega = omega_mu_opt(phi,epsi);
bndcond =  'Dirichlet'; 
ishomogeneous = true;  
RelaxationMethod = 'RB';
if strcmp(RelaxationMethod,'RB')||strcmp(RelaxationMethod,'GS')||strcmp(RelaxationMethod,'SymGS')
    omega = 1;
end
%%
param.CycType = CycType;
param.Omega_opt = omega;
param.v1 = nu1;
param.v2 = nu2;
param.RelaxationMethod = RelaxationMethod;
%[mu_opt,eig_set] = smooth_factor_TwoGrid(L_stencil,N_setup,omega,param);
[u_standard,u_set_standard,rnorm_standard,CPUTimestandard] =  MG_Base(L_stencil,N,nu1,nu2,n_cycles,CycType,omega,bndcond,ishomogeneous,RelaxationMethod,StopTolerance);
%%
param.CycType = CycType;
if strcmp(RelaxationMethod,'Jacobi')&&(nu1+nu2<=1)
    omega_modify = 8/13;
else
    omega_modify = omega;
end
if strcmp(RelaxationMethod,'RB')||strcmp(RelaxationMethod,'GS')||strcmp(RelaxationMethod,'SymGS')
    omega_modify = 1;
end
param.Omega_opt = omega_modify;
param.v1 = nu1;
param.v2 = nu2;
param.RelaxationMethod = RelaxationMethod;
[mu_opt,eig_set] = smooth_factor_ThreeGrid(L_stencil,N_setup,param);
%[mu_opt,eig_set] = smooth_factor_TwoGrid(L_stencil,N_setup,omega,param);
eig_set_abs = abs(eig_set);
if strcmp(RelaxationMethod,'RB')&&(nu1+nu2<=1)
    [c_opt,r_star] = CalculateOptimalFixedpara(0,Magfy*0.34); % 1/0 RB % 0.358
elseif strcmp(RelaxationMethod,'RB')&&(nu1+nu2>1)
    [c_opt,r_star] = CalculateOptimalFixedpara(0,Magfy*0.34^(nu1+nu2)); % 0.358^(nu1+nu2) 1/1 RB
else
    [c_opt,r_star] = CalculateOptimalFixedpara(min(eig_set_abs(:)),Magfy*max(eig_set_abs(:)));
end
[u_Acc,u_set_Acc,rnorm_Acc,CPUTimeAcc] = MG_BaseAcc(L_stencil,N,nu1,nu2,n_cycles,CycType,omega_modify,bndcond,c_opt,ishomogeneous,RelaxationMethod,StopTolerance);
%%
param.Omega_opt = omega;
% nu2 = 1;
% param.v2 = nu2;
[mu_opt,eig_set] = smooth_factor_ThreeGrid(L_stencil,N_setup,param);
[u_ChebAcc,u_set_ChebAcc,rnorm_ChebAcc,CPUTimeChebAcc] = ...
    MG_BaseChebAcc(L_stencil,N,nu1,nu2,n_cycles,CycType,param.Omega_opt,bndcond,ishomogeneous,RelaxationMethod,StopTolerance,min(real(eig_set)),Magfy*max(real(eig_set)));
%%
% PCG-MG or GMRES
if strcmp(RelaxationMethod,'RB')
    restart = n_cycles;
    isPre = true;
    [u_GMRESMG,rnorm_GMRESMG,CPUTime_GMRESMG] = gmres_wiki(L_stencil,N,restart,n_cycles,isPre,nu1,nu2,omega,CycType,bndcond,ishomogeneous,RelaxationMethod,StopTolerance);
else
    [u_PCGMG,u_set_PCGMG,rnorm_PCGMG,CPUTime_PCGMG] = PCG_MG(L_stencil,N,n_cycles,nu1,nu2,omega,CycType,bndcond,ishomogeneous,RelaxationMethod,StopTolerance);
end
%%
if strcmp(RelaxationMethod,'RB')
    
    figure
    semilogy(0:length(rnorm_standard)-1,rnorm_standard,'k-',0:length(rnorm_Acc)-1,rnorm_Acc,'r--',...
        0:length(rnorm_GMRESMG)-1,rnorm_GMRESMG,'b-.',0:length(rnorm_ChebAcc)-1,rnorm_ChebAcc,'g:','linewidth',2.2);
    set(gca,'fontsize',16)
    xlabel('Iterations','fontsize',20)
    ylabel('Residual Norm','fontsize',20)
    legend(['V(' num2str(nu1) ',' num2str(nu2) ')'],['Nesterov-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['GMRES-V(' num2str(nu1) ',' num2str(nu2) ')'],['Chebyshev-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        'location','best','fontsize',17)
    grid on
    
    figure
    semilogy(CPUTimestandard,rnorm_standard,'k-',CPUTimeAcc,rnorm_Acc,'r--',...
        CPUTime_GMRESMG,rnorm_GMRESMG,'b-.',CPUTimeChebAcc,rnorm_ChebAcc,'g:','linewidth',2.2);
    set(gca,'fontsize',16)
    xlabel('CPU Time (Seconds)','fontsize',20)
    ylabel('Residual Norm','fontsize',20)
    grid on
    legend(['V(' num2str(nu1) ',' num2str(nu2) ')'],['Nesterov-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['GMRES-V(' num2str(nu1) ',' num2str(nu2) ')'],['Chebyshev-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        'location','best','fontsize',17)
    
    figure
    plot(1:length(rnorm_standard)-1,rnorm_standard(2:end) ./ rnorm_standard(1:end-1),'k-',...
        1:length(rnorm_Acc)-1,rnorm_Acc(2:end) ./ rnorm_Acc(1:end-1),'r--',1:length(rnorm_GMRESMG)-1,...
        rnorm_GMRESMG(2:end)./rnorm_GMRESMG(1:end-1),'b-.',...
        1:length(rnorm_ChebAcc)-1,rnorm_ChebAcc(2:end) ./ rnorm_ChebAcc(1:end-1),'g:','linewidth',2.2);
    set(gca,'fontsize',14)
    xlabel('Iterations','fontsize',18)
    ylabel('Convergence Factor','fontsize',18)
    grid on
    legend(['V(' num2str(nu1) ',' num2str(nu2) ')'],['Nesterov-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['GMRES-V(' num2str(nu1) ',' num2str(nu2) ')'],['Chebyshev-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        'location','best','fontsize',17)
else
    
    figure
    semilogy(0:length(rnorm_standard)-1,rnorm_standard,'k-',0:length(rnorm_Acc)-1,rnorm_Acc,'r--',...
        0:length(rnorm_PCGMG)-1,rnorm_PCGMG,'b-.',0:length(rnorm_ChebAcc)-1,rnorm_ChebAcc,'g:',...
        'linewidth',2.2);
    set(gca,'fontsize',16)
    xlabel('Iterations','fontsize',20)
    ylabel('Residual Norm','fontsize',20)
    legend(['V(' num2str(nu1) ',' num2str(nu2) ')'],['Nesterov-V(' num2str(nu1) ',' num2str(nu2) ')'],['PCG-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['Chebyshev-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        'location','best','fontsize',17)
    grid on
    
    figure
    semilogy(CPUTimestandard,rnorm_standard,'k-',CPUTimeAcc,rnorm_Acc,'r--',CPUTime_PCGMG,rnorm_PCGMG,'b-.',...
        CPUTimeChebAcc,rnorm_ChebAcc,'g:','linewidth',2.2);
    set(gca,'fontsize',16)
    xlabel('CPU Time (Seconds)','fontsize',20)
    ylabel('Residual Norm','fontsize',20)
    grid on
    legend(['V(' num2str(nu1) ',' num2str(nu2) ')'],['Nesterov-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['PCG-V(' num2str(nu1) ',' num2str(nu2) ')'],['Chebyshev-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['PCG-V(' num2str(nu1) ',' num2str(nu2) ') - SymGS'],'location','best','fontsize',17)
    
    figure
    plot(1:length(rnorm_standard)-1,rnorm_standard(2:end) ./ rnorm_standard(1:end-1),'k-',...
        1:length(rnorm_Acc)-1,rnorm_Acc(2:end) ./ rnorm_Acc(1:end-1),'r--',...
        1:length(rnorm_PCGMG)-1,rnorm_PCGMG(2:end) ./ rnorm_PCGMG(1:end-1),'b-.',...
        1:length(rnorm_ChebAcc)-1,rnorm_ChebAcc(2:end) ./ rnorm_ChebAcc(1:end-1),'g:',...
        'linewidth',2.2);
    set(gca,'fontsize',16)
    xlabel('Iterations','fontsize',20)
    ylabel('Convergence Factor','fontsize',20)
    grid on
    legend(['V(' num2str(nu1) ',' num2str(nu2) ')'],['Nesterov-V(' num2str(nu1) ',' num2str(nu2) ')'],...
        ['PCG-V(' num2str(nu1) ',' num2str(nu2) ')'],['Chebyshev-V(' num2str(nu1) ',' num2str(nu2) ')'],...
            'location','best','fontsize',17)
    
end