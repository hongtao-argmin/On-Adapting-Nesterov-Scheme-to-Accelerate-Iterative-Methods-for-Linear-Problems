function [u,rnorm,CPUTime] = gmres_wiki(L_stencil,n,restart,MaxIter,isPre,nu1,nu2,w,CycType,bndcond,ishomogeneous,RelaxationMethod,StopTolerance)
% implement GMRES with preconditioner.
% Especially for PDE-based problem. 
% Inputs:
%   L_stencil: the stencil of PDE.
%   n: the size of grids.
%   restart: the number of restart.
%   MaxIter: maximal iteration number.
%   isPre: whether use precondition, if it is true, call V or W-Cycle.
%   nu1,nu2: the number of pre- and post-relaxation.
%   w: damping factor.
%   CycType: the type of Cycle. 1: V-Cycle; 2: W-Cycle.
%   bndcond: boundary condition, 'Periodic' and 'Dirichlet'.
%   ishomogeneous: true: zero f otherwise use a random initialization with a
%   corresponding f.
%   RelaxationMethod and StopTolerance: relaxation method and stop
%   tolerance.
% 
% Outputs:
%   u: the final solution.
%   rnorm: residual norm VS iteration.
%   CPUTime: CPU time VS iteration.
%
% Reference: check wikipedia.
% Tao.
% -------------------------------------------------------------------------
h = 1/n;                                % Mesh-size (the domain is the unit square)
[u,b] = initialdata(n,bndcond,ishomogeneous,L_stencil); % Define right-hand side, bounday conditions

if strcmp(bndcond,'Periodic')
    Ax = @(x)Lu(x,L_stencil,h,bndcond);
    total_A = n^2;
    size_A = [n n];
    vec2max = @(x)reshape(x,size_A);
    x = u;
elseif strcmp(bndcond,'Dirichlet')
    Ax = @(x)Lu([zeros(1,n+1);zeros(n-1,1) x zeros(n-1,1);zeros(1,n+1)],L_stencil,h,bndcond);
    total_A = (n-1)^2;
    size_A = [n-1,n-1];
    x = u(2:end-1,2:end-1);
    vec2max = @(x)[zeros(1,n+1);zeros(n-1,1) reshape(x,size_A) zeros(n-1,1);zeros(1,n+1)];
else
    error(['No such boundary condition - ' bndcond '\n']);
end
if isPre
    % if use preconditioned, MG method is called.
    % Performs a V cycle; recursive.
    % if more cycles are prefered, then run it many times.
    if strcmp(bndcond,'Periodic')
        M = @(x)(V_cycle(zeros(size_A),vec2max(x),h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod));
    elseif strcmp(bndcond,'Dirichlet')
        M = @(x)(V_cycle(zeros(size_A+2),vec2max(x),h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod));
    else
        error(['No such boundary condition - ' bndcond '\n']);
    end
    xx = x;
else
    M = @(x)(vec2max(x));
end

% set parameters
if isempty(restart)
    m = MaxIter;
else
    m = restart;
end
MaxIter = round(MaxIter/restart);

CPUTime = zeros(MaxIter*restart+1,1);
iters = 1;
r = b - Ax(x);
rnorm = norm(r(:))/n;
isbreak = false;
StartTime = tic;
for iter = 1:MaxIter
    % use x as the initial vector
    r = b - Ax(x);
    if strcmp(bndcond,'Dirichlet')
        r = r(2:end-1,2:end-1);
        r = M(r(:));
        r = r(2:end-1,2:end-1);
    elseif strcmp(bndcond,'Periodic')
        r = M(r(:));
    else
        error(['No such boundary condition - ' bndcond '\n']);
    end
    
    % initialize the 1D vectors
    sn = zeros(m, 1);
    cs = zeros(m, 1);
    e1 = zeros(total_A, 1);
    e1(1) = 1;
    r_norm = norm(r(:));
    Q(:,1) = r(:)/r_norm;
    beta = r_norm*e1;
    
    for k = 1:m
        
        % run arnoldi
        [H(1:k+1, k),Q(:, k+1)] = arnoldi(Ax,Q,k,M,size_A,bndcond);
        
        % eliminate the last element in H ith row and update the rotation matrix
        [H(1:k+1, k),cs(k),sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k);
        
        % update the residual vector
        beta(k + 1) = -sn(k) * beta(k);
        beta(k)     = cs(k) * beta(k);
        errorII       = abs(beta(k + 1))/n;
        
        % save the error
        CPUTime(iters+1) = toc(StartTime);
        iters = iters+1;
        if isPre
            % calculate the result
            yy = H(1:k, 1:k)\beta(1:k);
            xx(:) = x(:) + Q(:, 1:k) * yy;
            rr = b-Ax(xx);
            errorII = norm(rr(:))/n;
        end
        rnorm = [rnorm; errorII];
        StartTime = tic;
        if (errorII <= StopTolerance)
            isbreak = true;
            break;
        end
    end
    % calculate the result
    y = H(1:k, 1:k)\beta(1:k);
    x(:) = x(:) + Q(:, 1:k) * y;
    if isbreak
        CPUTime(iters+1) = toc(StartTime);
        break;
    end
end
if ~isbreak
    CPUTime(iters+1) = toc(StartTime);
end
CPUTime(iters+1) = CPUTime(iters)+CPUTime(iters+1);
CPUTime(iters+1) = 0;
CPUTime = cumsum(CPUTime(1:iters));
if strcmp(bndcond,'Dirichlet')
    u(2:end-1,2:end-1) = x;
else
    u = x;
end

return

%----------------------------------------------------%
%                  Arnoldi Function                  %
%----------------------------------------------------%

function [h, q] = arnoldi(Ax, Q, k,M,size_A,bndcond)

q = Ax(reshape(Q(:,k),size_A));

if strcmp(bndcond,'Periodic')
    q = M(q);
    q = q(:);
elseif strcmp(bndcond,'Dirichlet')
    q = q(2:end-1,2:end-1);
    q = M(q);
    q = q(2:end-1,2:end-1);
    q = q(:);
else
    error(['No such boundary condition - ' bndcond '\n']);
end

for i = 1:k
    h(i) = q' * Q(:, i);
    q = q - h(i) * Q(:, i);
end
h(k + 1) = norm(q);
q = q / h(k + 1);

return

%---------------------------------------------------------------------%
%                  Applying Givens Rotation to H col                  %
%---------------------------------------------------------------------%
function [h, cs_k, sn_k] = apply_givens_rotation(h, cs, sn, k)
% apply for ith column
for i = 1:k-1
    temp   =  cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
end

% update the next sin cos values for rotation
[cs_k,sn_k] = givens_rotation(h(k), h(k + 1));

% eliminate H(i + 1, i)
h(k) = cs_k * h(k) + sn_k * h(k + 1);
h(k + 1) = 0.0;
return

%%----Calculate the Given rotation matrix----%%
function [cs, sn] = givens_rotation(v1, v2)
if (v1 == 0)
    cs = 0;
    sn = 1;
else
    t = sqrt(v1^2 + v2^2);
    cs = abs(v1) / t;
    sn = cs * v2 / v1;
end
return