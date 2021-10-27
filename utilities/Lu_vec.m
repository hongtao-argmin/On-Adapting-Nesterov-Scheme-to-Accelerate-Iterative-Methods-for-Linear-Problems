function f = Lu_vec(u,M,N,L_stencil,h,bndcond) % Computes L(u)

if strcmp(bndcond,'Dirichlet')
    u = reshape(u,M-1,N-1);
    f = zeros(size(u)+2);
    f(2:end-1,2:end-1) = conv2(u,L_stencil/h^2,'same');
elseif strcmp(bndcond,'Periodic')
    u = reshape(u,M,N);
    f = imfilter(u,L_stencil/h^2,'circular','conv');
else
    error(['No such boundary condition - ' bndcond '\n'])
end
f = f(:);
return

