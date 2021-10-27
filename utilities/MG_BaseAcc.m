
function [u,u_set,rnorm,CPUTime] = MG_BaseAcc(L_stencil,n,nu1,nu2,n_cycles,CycType,w,bndcond,accpara,ishomogeneous,RelaxationMethod,StopTolerance)

h = 1/n;
[u,f] = initialdata(n,bndcond,ishomogeneous,L_stencil);
r = f - Lu(u,L_stencil,h,bndcond);                   % The initial residual.
rnorm = norm(r(:))/n;                                % This will store the properly normalized
%t_k = 1;
u_set{1} = u;
u_k = u;
y_k = u;
CPUTime = zeros(n_cycles+1,1);
for iters = 1:n_cycles
    StartTime = tic;
    u_k_1 = V_cycle(y_k,f,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);
    %     if mod(iters,100)==0
    %        y_k = u_k_1;
    %     else
    %t_k_1 = (1+sqrt(1+4*t_k^2))/2;
    %y_k = u_k_1 + (iters-2)/(iters+1)*(u_k_1-u_k);%(t_k-1)/t_k_1
    %t_k = t_k_1;
    y_k = u_k_1 + accpara*(u_k_1-u_k);
    
    %     end
    u_k = u_k_1;
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