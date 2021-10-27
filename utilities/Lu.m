
function f = Lu(u,L_stencil,h,bndcond) % Computes L(u)

if strcmp(bndcond,'Dirichlet')
    f = zeros(size(u));
    f(2:end-1,2:end-1) = conv2(u(2:end-1,2:end-1),L_stencil/h^2,'same');
elseif strcmp(bndcond,'Periodic')
    f = imfilter(u,L_stencil/h^2,'circular','conv');
else
    error(['No such boundary condition - ' bndcond '\n'])
end

return
