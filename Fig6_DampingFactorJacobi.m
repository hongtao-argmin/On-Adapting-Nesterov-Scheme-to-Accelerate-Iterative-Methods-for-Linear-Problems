%
% Check optimal damping factor for V(1,0) with damped Jacobi as the
% relaxation.
%
%%
clc;
clear all;
close all;
%%
addpath('./utilities')
%%
r_base_set = 0.6;
omega_set = 0.01:0.01:0.99;%0.5:0.01:0.9;%
ratio_set = zeros(length(omega_set),1);
r_opt_set_Nesterov = zeros(length(omega_set),1);
r_opt_set_base = zeros(length(omega_set),1);
r_opt_set_Nesterov_practice = zeros(length(omega_set),1);
r_opt_set_base_practice = zeros(length(omega_set),1);
%%
N = 256;
nu1 = 1;
nu2 = 0;
n_cycles = 2000;
epsi = 1;
phi = 0;
CycType = 1;
StopTolerance = 1e-8;
Magfy = 1;
%%
C = cos(phi);
S = sin(phi);
L_stencil = -[-0.5*(1-epsi)*C*S epsi*C^2+S^2 0.5*(1-epsi)*C*S;...
    C^2+epsi*S^2 -2*(1+epsi) C^2+epsi*S^2; ...
    0.5*(1-epsi)*C*S epsi*C^2+S^2 -0.5*(1-epsi)*C*S];
bndcond =  'Dirichlet';  
ishomogeneous = true;  
RelaxationMethod = 'Jacobi'; 
%%
for index_i = 1:length(omega_set)
    b_1 = 1-2*omega_set(index_i);
    b_N = 1-0.5*omega_set(index_i);
    ratio_set(index_i) = b_1/b_N;
    r_opt_set_base(index_i) = max(abs(b_1),abs(b_N));
    [c_opt,r_opt_set_Nesterov(index_i)] = CalculateOptimalFixedpara(b_1,b_N);
    [u_standard,u_set_standard,rnorm_standard,CPUTimestandard] =  MG_Base(L_stencil,N,nu1,nu2,n_cycles,CycType,omega_set(index_i),bndcond,ishomogeneous,RelaxationMethod,StopTolerance);
    temp = rnorm_standard(2:end) ./ rnorm_standard(1:end-1);
    r_opt_set_base_practice(index_i) = geomean(temp(end-5:end));
    [u_Acc,u_set_Acc,rnorm_Acc,CPUTimeAcc] = MG_BaseAcc(L_stencil,N,nu1,nu2,n_cycles,CycType,omega_set(index_i),bndcond,c_opt,ishomogeneous,RelaxationMethod,StopTolerance);
    temp = rnorm_Acc(2:end) ./ rnorm_Acc(1:end-1);
    r_opt_set_Nesterov_practice(index_i) = geomean(temp(end-5:end));
end
%%
figure
colororder({'k','b'})
yyaxis left
plot(omega_set,r_opt_set_base_practice,'-x',omega_set,r_opt_set_base,'-',omega_set,r_opt_set_Nesterov_practice,'-o',omega_set,r_opt_set_Nesterov,'--','linewidth',2)
xlabel('Damping Factor')
axis tight
grid on
set(gca,'fontsize',11)
xlabel('Damping Factor','fontsize',20)
ylabel('$r^*$','Interpreter','latex','fontsize',19)
ylim([min(r_opt_set_Nesterov) max(r_opt_set_base)])
yyaxis right
plot(omega_set,ratio_set,'-.','linewidth',2)
ylabel('$\frac{b_1}{b_N}$','Interpreter','latex','fontsize',19)
ylim([min(ratio_set) max(ratio_set)])
legend('V(1,0)-Practice','V(1,0)','Nesterov-V(1,0)-Practice','Nesterov-V(1,0)','$\frac{b_1}{b_N}$','Interpreter','latex','location','best','fontsize',17)