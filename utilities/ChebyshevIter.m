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

% % stra. I
% omega = (2*symbol-b_N-b_1)/(b_N-b_1);
% omega_ground = (2-b_N-b_1)/(b_N-b_1);
%
% T_k = 1;
% T_k_1 = omega;
%
% T_k_ground = 1;
% T_k_1_ground = omega_ground;
% for iter = 2:niter
%     T_k_1_post = 2*omega*T_k_1-T_k;
%     T_k = T_k_1;
%     T_k_1 = T_k_1_post;
%
%     T_k_1_post_ground = 2*omega_ground*T_k_1_ground-T_k_ground;
%     T_k_ground = T_k_1_ground;
%     T_k_1_ground = T_k_1_post_ground;
%
% end
% r_final = max(abs(rootsc(niter,abs(T_k_1)/abs(T_k_1_post_ground))));

% stra. II -- use matrix form
% niter = 50;
% gamma = 2/(2-b_N-b_1);
% sigma = (b_N-b_1)/(2-b_N-b_1);
% rho_k = 1/(1-0.5*sigma^2);
% r_final_set = zeros(niter-1,1);
% for iter = 2:niter
%     rho_k_1 = 1/(1-0.25*sigma^2*rho_k);
%     r_final_set(iter-1) = vrho([rho_k_1*(gamma*symbol+1-gamma) 1-rho_k_1;1 0]);
%     rho_k = rho_k_1;
% end
% 
% r_final = max(r_final_set);

% stra. III -- use extreme rho form

gamma = 2/(2-b_N-b_1);
sigma = (b_N-b_1)/(2-b_N-b_1);

rho_final = 2/(1+sqrt(1-sigma^2));
r_final = vrho([rho_final*(gamma*symbol+1-gamma) 1-rho_final;1 0]);

% stra. IV -- use matrix form iteration
% niter = 50;
% gamma = 2/(2-b_N-b_1);
% sigma = (b_N-b_1)/(2-b_N-b_1);
% rho_k = 1/(1-0.5*sigma^2);
% MM = [rho_k*(gamma*symbol+1-gamma) 1-rho_k;1 0];
% for iter = 2:niter
%     rho_k_1 = 1/(1-0.25*sigma^2*rho_k);
%     MM = MM*[rho_k_1*(gamma*symbol+1-gamma) 1-rho_k_1;1 0];
%     rho_k = rho_k_1;
% end
% r_final = (vrho(MM))^(1/niter);

if isPrint
    fprintf('The final conv. factor for symbol %f%+fj is: %f \n',real(symbol),imag(symbol),r_final);
end
