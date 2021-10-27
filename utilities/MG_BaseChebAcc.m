
function [u,u_set,rnorm,CPUTime] = MG_BaseChebAcc(L_stencil,n,nu1,nu2,n_cycles,CycType,w,bndcond,ishomogeneous,RelaxationMethod,StopTolerance,b_1,b_N)

h = 1/n;
[u,f] = initialdata(n,bndcond,ishomogeneous,L_stencil);
r = f - Lu(u,L_stencil,h,bndcond);                   % The initial residual.
rnorm = norm(r(:))/n;                                % This will store the properly normalized
gamma = 2/(2-b_1-b_N);
sigma = (b_N-b_1)/(2-b_1-b_N);
%rho_k = 1;
rho_k_1 = 1/(1-0.5*sigma^2);
u_set{1} = u;
u_k = u;
CPUTime = zeros(n_cycles+1,1);
for iters = 1:n_cycles
    StartTime = tic;
    y_k = V_cycle(u_k,f,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);
    
    if iters==1
        u_k_1 = gamma*y_k+(1-gamma)*u_k;
        u_k_min_k = u_k;
        u_k = u_k_1;
    else
        u_k_1 = rho_k_1*(gamma*y_k+(1-gamma)*u_k)+(1-rho_k_1)*u_k_min_k;
        u_k_min_k = u_k;
        u_k = u_k_1;
        rho_k_1 = 1/(1-0.25*sigma^2*rho_k_1);
    end
    
    CPUTime(iters+1) = toc(StartTime);
    r = f - Lu(u_k,L_stencil,h,bndcond);
    rnorm = [rnorm, norm(r(:))/n];
    u_set{iters+1} = u_k;
    if rnorm(end)<StopTolerance
        break;
    end
end
u = u_k;
n_cycles = length(rnorm)-1;
CPUTime = cumsum(CPUTime);
CPUTime = CPUTime(1:n_cycles+1);
temp = rnorm(2:end) ./ rnorm(1:end-1);
fprintf('%f\n',geomean(temp(end-5:end)))
% figure
% semilogy([0:n_cycles],rnorm,'k-');
% title('Residual Norm History')
% xlabel('Iterations')
% ylabel('Residual Norm')
% % Plot residual norm convergence factors.
% figure
% plot([1:n_cycles],rnorm(2:end) ./ rnorm(1:end-1),'k-');
% title('Residual Convergence Factor Per Cycle')
% xlabel('Iterations')
% ylabel('Convergence Factor')

end