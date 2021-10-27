function u = Jacobi(u,f,h,w,L_stencil,bndcond)            % Jacobi relaxation with damping w

if sum(u(:))==0
    u = u + (h^2*w/L_stencil(2,2))*f;  % Assumes f = L(u) on boundaries
else
    u = u + (h^2*w/L_stencil(2,2))*(f - Lu(u,L_stencil,h,bndcond));  % Assumes f = L(u) on boundaries
end

return