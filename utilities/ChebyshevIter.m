% *** Chebyshev polynomial ***
%
% Estimate the Conv. Rate for a given symbol and number of iterations.
% symbol can be a complex number.
%
function r_final = ChebyshevIter(symbol,b_1,b_N,isPrint)% niter,

% if niter<5
%     % we need enough larger iteration to obtain a more accurate solution.
%     niter = 50;
% end
if nargin<=3
    isPrint = false;
end

gamma = 2/(2-b_N-b_1);
sigma = (b_N-b_1)/(2-b_N-b_1);

rho_final = 2/(1+sqrt(1-sigma^2));
r_final = vrho([rho_final*(gamma*symbol+1-gamma) 1-rho_final;1 0]);

if isPrint
    fprintf('The final conv. factor for symbol %f%+fj is: %f \n',real(symbol),imag(symbol),r_final);
end
