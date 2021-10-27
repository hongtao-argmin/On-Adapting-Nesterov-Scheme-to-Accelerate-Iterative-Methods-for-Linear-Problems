% PCG with MG as preconditioner.
%---------------------------------------------------------------------------------------
function [u,u_set,rnorm,CPUTime] = PCG_MG(L_stencil,n,MaxIter,nu1,nu2,w,CycType,bndcond,ishomogeneous,RelaxationMethod,StopTolerance)

h = 1/n;                               % Mesh-size (the domain is the unit square)

[u,f] = initialdata(n,bndcond,ishomogeneous,L_stencil); % Define right-hand side, bounday conditions,
r_k_1 = u;
p_k_1 = u;
% and initial guess for the solution
r = f - Lu(u,L_stencil,h,bndcond);        % The initial residual
rnorm = norm(r(:))/n;             % This will store the residual vector norm after each cycle
u_set{1} = u;
CPUTime = zeros(MaxIter+1,1);
for iters = 1:MaxIter
    StartTime = tic;
    % CG
    if strcmp(bndcond,'Dirichlet')
        if iters ==1
            r_k = f-Lu(u,L_stencil,h,bndcond);
            z_k = V_cycle(zeros(size(u)),r_k,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);    % Performs a V cycle; recursive.
            r_k_z_k = sum(sum((r_k(2:end-1,2:end-1).*z_k(2:end-1,2:end-1))));
            p_k = z_k;
        else
            r_k_1(2:end-1,2:end-1)  = r_k(2:end-1,2:end-1) -alpha_k*Ap_k(2:end-1,2:end-1);
            z_k_1 = V_cycle(zeros(size(u)),r_k_1,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);
            r_k_1_z_k_1 = sum(sum(z_k_1(2:end-1,2:end-1).*r_k_1(2:end-1,2:end-1)));
            beta_k = r_k_1_z_k_1/r_k_z_k;
            p_k_1(2:end-1,2:end-1) = z_k_1(2:end-1,2:end-1)+beta_k*p_k(2:end-1,2:end-1);
            r_k_z_k = r_k_1_z_k_1;
            p_k = p_k_1;
            r_k = r_k_1;
        end
        Ap_k = Lu(p_k,L_stencil,h,bndcond);
        alpha_k = r_k_z_k/sum(sum(p_k(2:end-1,2:end-1).*Ap_k(2:end-1,2:end-1)));
        u(2:end-1,2:end-1) = u(2:end-1,2:end-1)+alpha_k*p_k(2:end-1,2:end-1);
    elseif strcmp(bndcond,'Periodic')
        if iters ==1
            r_k = f-Lu(u,L_stencil,h,bndcond);
            z_k = V_cycle(zeros(size(u)),r_k,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);    % Performs a V cycle; recursive.
            r_k_z_k = sum(sum((r_k.*z_k)));
            p_k = z_k;
        else
            r_k_1  = r_k -alpha_k*Ap_k;
            z_k_1 = V_cycle(zeros(size(u)),r_k_1,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);
            r_k_1_z_k_1 = sum(sum(z_k_1.*r_k_1));
            beta_k = r_k_1_z_k_1/r_k_z_k;
            p_k_1 = z_k_1+beta_k*p_k;
            r_k_z_k = r_k_1_z_k_1;
            p_k = p_k_1;
            r_k = r_k_1;
        end
        Ap_k = Lu(p_k,L_stencil,h,bndcond);
        alpha_k = r_k_z_k/sum(sum(p_k.*Ap_k));
        u = u+alpha_k*p_k;
    else
        error(['No such boundary condition - ' bndcond '\n'])
    end
    CPUTime(iters+1) = toc(StartTime);
    u_set{iters+1} = u;
    r = f - Lu(u,L_stencil,h,bndcond);
    rnorm = [rnorm; norm(r(:))/n];
    if rnorm(end)<StopTolerance
        CPUTime = CPUTime(1:iters+1);
        break;
    end
end

MaxIter = length(rnorm)-1;
CPUTime = cumsum(CPUTime);
CPUTime = CPUTime(1:MaxIter+1);
temp = rnorm(2:end) ./ rnorm(1:end-1);
fprintf('%f\n',geomean(temp(end-5:end)))
% figure
% semilogy([0:MaxIter],rnorm,'k-');
% title('Residual Norm History')
% xlabel('Iterations')
% ylabel('Residual Norm')
% % Plot residual norm convergence factors.
% figure
% plot([1:MaxIter],rnorm(2:end) ./ rnorm(1:end-1),'k-');
% title('Residual Convergence Factor Per Cycle')
% xlabel('Iterations')
% ylabel('Convergence Factor')
return
