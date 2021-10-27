%%
clc
clear all;
close all;
%%
addpath('./utilities')
%%
isRev = false;
b_N = 0.9;  % 0.34;
b_1 = -0.3; %-0.5; %-0.5; % -0.5; % 0;%
if isRev
    temp = b_1;
    b_1 = -b_N;
    b_N = -temp;
end
% -------------------------------------------------------------------------
% if this is true then we test when the ACF is better than the base alg.
% otherwise we obtain the regions which are not worse than the accelerated ACF.
isStable = false;
epsilon = 0;
%%
if isStable
    c_opt = CalculateOptimalFixedpara(b_1,b_N);
    r_star_Nest = b_N;
else
    [c_opt,r_star_Nest] = CalculateOptimalFixedpara(b_1,b_N);
end
r_star_cheby = ChebyshevIter(b_1,b_1,b_N);
%%
% Compare diff. complex regions of Chebyshev acceleration and Nesterov's scheme.
magify = 1;
N = 3e2;
x_modu = 0:magify*b_N/N:1-magify*b_N/N;
y_arg = -pi:pi/N:pi;
Nesterovregion = [];
Chebyregion = [];
ChebyregionworseNest = [];
for index_i = 1:length(x_modu)
    for index_j = 1:length(y_arg)
        b_cmp = x_modu(index_i)*exp(1i*y_arg(index_j));
        if abs(imag(b_cmp))<epsilon
            b_cmp = real(b_cmp);
        end
        r_cheby = ChebyshevIter(b_cmp,b_1,b_N);
        r_nes = NesterovIter(b_cmp,c_opt,[],[]);
        if r_cheby<=r_star_Nest+epsilon
            Chebyregion = [Chebyregion b_cmp];
        else
            if r_nes<=r_star_Nest+epsilon
                ChebyregionworseNest = [ChebyregionworseNest b_cmp];
            end
        end
        if r_nes<=r_star_Nest+epsilon
            Nesterovregion = [Nesterovregion b_cmp];
        end
    end
end
%%
figure
plot(real(Nesterovregion), imag(Nesterovregion), 'b.')
zoom on
grid on
hold on
plot( [-1 1],[0 0], '-k', 'linewidth', 2.2)
plot( [0 0],[-1 1], '-k', 'linewidth', 2.2)
if isRev
    if b_N<=-1/3*b_1-epsilon
        plot(abs(-1/3*b_1)*cos(y_arg),abs(-1/3*b_1)*sin(y_arg),'-r','linewidth',2.2)
    else
        plot(abs(b_N)*cos(y_arg),abs(b_N)*sin(y_arg),'-r','linewidth',2.2)
    end
else
    if b_1>=-1/3*b_N-epsilon
        plot(abs(-1/3*b_N)*cos(y_arg),abs(-1/3*b_N)*sin(y_arg),'-r','linewidth',2.2)
    else
        plot(abs(b_1)*cos(y_arg),abs(b_1)*sin(y_arg),'-r','linewidth',2.2)
    end
end
xticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xticklabels({'-1' '-0.75' '-0.5' '-0.25' '0' '0.25' '0.5' '0.75' '1'})
yticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
yticklabels({'-1' '-0.75' '-0.5' '-0.25' '0' '0.25' '0.5' '0.75' '1'})
set(gca,'FontSize',13)
axis square
xlabel('Real','FontSize',18)
ylabel('Imaginary','FontSize',18)
%%
figure
plot(real(ChebyregionworseNest), imag(ChebyregionworseNest), 'b.')
zoom on
grid on
hold on
plot(real(Chebyregion), imag(Chebyregion), 'c.')
plot([-1 1],[0 0], '-k', 'linewidth', 2.2)
plot([0 0],[-1 1], '-k', 'linewidth', 2.2)
if isRev
    if b_N<=-1/3*b_1-epsilon
        plot(abs(-1/3*b_1)*cos(y_arg),abs(-1/3*b_1)*sin(y_arg),'-r','linewidth',2.2)
    else
        plot(abs(b_N)*cos(y_arg),abs(b_N)*sin(y_arg),'-r','linewidth',2.2)
    end
else
    if b_1>=-1/3*b_N-epsilon
        plot(abs(-1/3*b_N)*cos(y_arg),abs(-1/3*b_N)*sin(y_arg),'-r','linewidth',2.2)
    else
        plot(abs(b_1)*cos(y_arg),abs(b_1)*sin(y_arg),'-r','linewidth',2.2)
    end
end
xticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
xticklabels({'-1' '-0.75' '-0.5' '-0.25' '0' '0.25' '0.5' '0.75' '1'})
yticks([-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1])
yticklabels({'-1' '-0.75' '-0.5' '-0.25' '0' '0.25' '0.5' '0.75' '1'})
set(gca,'FontSize',13)
axis square
xlabel('Real','FontSize',18)
ylabel('Imaginary','FontSize',18)
%%
% test the change of ACF with respect to arguments for a given modulus
% Test on Nesterov and Chebyshev to see their difference.
b_N = 0.9;
b_1 = -0.5;
N = 3e2; 
[c_opt,r_star_Nest] = CalculateOptimalFixedpara(b_1,b_N);
epsilon = 0;
if abs(b_1)<=b_N
    if b_1>=-1/3*b_N-epsilon
        X_modus_set = [1/3*b_N;0.5*b_N;0.8*b_N];
    else
        X_modus_set = [abs(b_1);0.7*b_N;0.9*b_N];
    end
else
    if b_N<=-1/3*b_1+epsilon
        X_modus_set = [-1/3*b_1;-0.5*b_1;-0.8*b_1];
    else
        X_modus_set = [-1/3*b_1;-0.7*b_1;-0.9*b_1];
    end
end
y_arg = 0:pi/N:pi;
r_set_nest = zeros(length(y_arg),length(X_modus_set));
r_set_cheb = zeros(length(y_arg),length(X_modus_set));
for index_i = 1:length(X_modus_set)
    for index_j = 1:length(y_arg)
        b_cmp = X_modus_set(index_i)*exp(1i*y_arg(index_j));
        r_set_nest(index_j,index_i) = NesterovIter(b_cmp,c_opt,[],[]);
        r_set_cheb(index_j,index_i) = ChebyshevIter(b_cmp,b_1,b_N);
    end
end
%% 
Star_rate = r_star_Nest*ones(length(y_arg),1);
figure
if b_1>=-b_N/3 || b_N<=-b_1/3
    plot(y_arg,Star_rate, '-k', 'linewidth', 1.5,'DisplayName',['Nesterov ACF: ' '$r^* = $' ' ' num2str(Star_rate(1))],'linewidth',2.2)
else
    plot(y_arg,Star_rate, '-k', 'linewidth', 1.5,'DisplayName',['Nesterov ACF: ' '$r^* = $' ' ' num2str(Star_rate(1))],'linewidth',2.2)
end
hold on
for kk = 1:length(X_modus_set)
    plot(y_arg,r_set_nest(:,kk),'--','linewidth',1.5,'DisplayName',['Nesterov: ' '$\bar b^c = $' ' ' num2str(X_modus_set(kk))],'linewidth',2.2)
    plot(y_arg,r_set_cheb(:,kk),':','linewidth',1.5,'DisplayName',['Chebyshev: ' '$\bar b^c = $' ' ' num2str(X_modus_set(kk))],'linewidth',2.2)
end
set(gca,'FontSize',13)
set(legend,'FontSize',14)
axis([0 3.5 0.35 2*max(abs(b_1),abs(b_N))])

grid on
hl = legend('location','best');
set(hl, 'Interpreter','latex');
xlabel('$\theta$','interpret','latex','fontsize',20)
if b_1>=-b_N/3 || b_N<=-b_1/3
    ylabel('$r_{c_{top}}(b^c)$','interpret','latex','fontsize',20)
else
    ylabel('$r_{c_{mid}}(b^c)$','interpret','latex','fontsize',20)
end