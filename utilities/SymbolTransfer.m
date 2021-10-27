% This function is written to obtain the symbol of the restriction with different orders in 2D.
%
% LF and HF are 2th: 0.5[0.5 1 0.5].
% LF and HF are 4th: 0.5[-1/16 0 9/16 1 9/16 0 -1/16].
% LF and HF are 6th: 0.5[3/256 0 -25/256 0 150/256 1 150/256 0 -25/256 0 3/256].
% LF and HF are 8th: 0.5[-5/2048 0 49/2048 0 -245/2048 0 1225/2048 1 1225/2048 0 -245/2048 0 49/2048 0 -5/2048].
% ...
% [a b c d...] = A\[0.5;0;0;0;...].
% A = [1 1 1 1 1...;1 3^2 5^2 7^2 8^2...;1 3^4....;1 3^6...;1 3^8...]
% Author: Tao Hong, contact: hongtao@cs.technion.ac.il.
% ========================================================================
function [I_h_H_1,I_h_H_2,I_h_H_3,I_h_H_4] = SymbolTransfer(theta_1,theta_2,RestrictionMethod,Type)

if nargin<=3
    Type = 'Pi';
end

[theta_1_hat,theta_2_hat] = AngleAlias(theta_1,theta_2,Type);

switch RestrictionMethod
    case 'Bilinear'
        I_h_H_1 = 0.25*(1+cos(theta_1))*(1+cos(theta_2));
        I_h_H_2 = 0.25*(1+cos(theta_1_hat))*(1+cos(theta_2_hat));
        I_h_H_3 = 0.25*(1+cos(theta_1_hat))*(1+cos(theta_2));
        I_h_H_4 = 0.25*(1+cos(theta_1))*(1+cos(theta_2_hat));
    case 'Cubic'
        I_h_H_1 = 0.25*(1+9/8*cos(theta_1)-0.125*cos(3*theta_1))*(1+1.125*cos(theta_2)-0.125*cos(3*theta_2));
        I_h_H_2 = 0.25*(1+1.125*cos(theta_1_hat)-0.125*cos(3*theta_1_hat))*(1+1.125*cos(theta_2_hat)-0.125*cos(3*theta_2_hat));
        I_h_H_3 = 0.25*(1+1.125*cos(theta_1_hat)-1/8*cos(3*theta_1_hat))*(1+1.125*cos(theta_2)-0.125*cos(3*theta_2));
        I_h_H_4 = 0.25*(1+1.125*cos(theta_1)-0.125*cos(3*theta_1))*(1+1.125*cos(theta_2_hat)-0.125*cos(3*theta_2_hat));
    case '6th'
        a = 2*150/256;
        b = -2*25/256;
        c = 6/256;
        I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2));
        I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat));
        I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2));
        I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat));
    case '8th'
        a = 2*1225/2048;
        b = -2*245/2048;
        c = 2*49/2048;
        d = -10/2048;
        I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2));
        I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat));
        I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2));
        I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat));
    case '10th'
        x = [1 1 1 1 1;1 3^2 5^2 7^2 9^2;1 3^4 5^4 7^4 9^4;1 3^6 5^6 7^6 9^6;1 3^8 5^8 7^8 9^8]\[1;0;0;0;0];
        a = x(1);
        b = x(2);
        c = x(3);
        d = x(4);
        e = x(5);
        I_h_H_1 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1)+e*cos(9*theta_1))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2)+e*cos(9*theta_2));
        I_h_H_2 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat)+e*cos(9*theta_1_hat))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat)+e*cos(9*theta_2_hat));
        I_h_H_3 = 0.25*(1+a*cos(theta_1_hat)+b*cos(3*theta_1_hat)+c*cos(5*theta_1_hat)+d*cos(7*theta_1_hat)+e*cos(9*theta_1_hat))*(1+a*cos(theta_2)+b*cos(3*theta_2)+c*cos(5*theta_2)+d*cos(7*theta_2)+e*cos(9*theta_2));
        I_h_H_4 = 0.25*(1+a*cos(theta_1)+b*cos(3*theta_1)+c*cos(5*theta_1)+d*cos(7*theta_1)+e*cos(9*theta_1))*(1+a*cos(theta_2_hat)+b*cos(3*theta_2_hat)+c*cos(5*theta_2_hat)+d*cos(7*theta_2_hat)+e*cos(9*theta_2_hat));
    
    case {'2','4','6','8','10','12','14','16','18','20','22','24','26','28','30','32','34','36','38','40','42','44','46','50'}
        
        [I_h_H_1,I_h_H_2,I_h_H_3,I_h_H_4] = HighOrderSymbol(str2double(RestrictionMethod),theta_1,theta_2,theta_1_hat,theta_2_hat);
        
    case 'LF4thHF2th'
        I_h_H_1 = (0.125)^2*(5+4*cos(theta_1)-cos(2*theta_1))*(5+4*cos(theta_2)-cos(2*theta_2));
        I_h_H_2 = (0.125)^2*(5+4*cos(theta_1_hat)-cos(2*theta_1_hat))*(5+4*cos(theta_2_hat)-cos(2*theta_2_hat));
        I_h_H_3 = (0.125)^2*(5+4*cos(theta_1_hat)-cos(2*theta_1_hat))*(5+4*cos(theta_2)-cos(2*theta_2));
        I_h_H_4 = (0.125)^2*(5+4*cos(theta_1)-cos(2*theta_1))*(5+4*cos(theta_2_hat)-cos(2*theta_2_hat));
    case 'LF2thHF4th'
        I_h_H_1 = (0.125)^2*(3+4*cos(theta_1)+cos(2*theta_1))*(3+4*cos(theta_2)+cos(2*theta_2));
        I_h_H_2 = (0.125)^2*(3+4*cos(theta_1_hat)+cos(2*theta_1_hat))*(3+4*cos(theta_2_hat)+cos(2*theta_2_hat));
        I_h_H_3 = (0.125)^2*(3+4*cos(theta_1_hat)+cos(2*theta_1_hat))*(3+4*cos(theta_2)+cos(2*theta_2));
        I_h_H_4 = (0.125)^2*(3+4*cos(theta_1)+cos(2*theta_1))*(3+4*cos(theta_2_hat)+cos(2*theta_2_hat));
    case 'Ideal'
        I_h_H_1 = 1;
        I_h_H_2 = 0;
        I_h_H_3 = 0;
        I_h_H_4 = 0;
    otherwise
        error([RestrictionMethod ' is not defined \n'])          
end

return

function [theta_1_hat,theta_2_hat] = AngleAlias(theta_1,theta_2,Type)

if strcmp(Type,'Pi_2')
    if theta_1>=0&&theta_2>=0
        theta_1_hat = theta_1-pi/2;
        theta_2_hat = theta_2-pi/2;
    elseif theta_1>=0&&theta_2<0
        theta_1_hat = theta_1-pi/2;
        theta_2_hat = theta_2+pi/2;
    elseif theta_1<0&&theta_2>=0
        theta_1_hat = theta_1+pi/2;
        theta_2_hat = theta_2-pi/2;
    elseif theta_1<0&&theta_2<0
        theta_1_hat = theta_1+pi/2;
        theta_2_hat = theta_2+pi/2;
    else
        error('It has an error\n');
    end
elseif strcmp(Type,'Pi')
    if theta_1>=0&&theta_2>=0
        theta_1_hat = theta_1-pi;
        theta_2_hat = theta_2-pi;
    elseif theta_1>=0&&theta_2<0
        theta_1_hat = theta_1-pi;
        theta_2_hat = theta_2+pi;
    elseif theta_1<0&&theta_2>=0
        theta_1_hat = theta_1+pi;
        theta_2_hat = theta_2-pi;
    elseif theta_1<0&&theta_2<0
        theta_1_hat = theta_1+pi;
        theta_2_hat = theta_2+pi;
    else
        error('It has an error\n');
    end
else
    error('We do not have such an algle transform')
end

return

function [I_h_H_1,I_h_H_2,I_h_H_3,I_h_H_4] = HighOrderSymbol(order,theta_1,theta_2,theta_1_hat,theta_2_hat)

% A = ones(1,order/2);
% for index_i = 1:order/2-1
%     A = [A;(1:2:order-1).^(2*index_i)];
% end
A = fliplr(vander((1:2:order-1).^2))';
for index = 2:size(A,1)
    A(index,:) = A(index,:)/max(A(index,:));
end
b = [1;zeros(order/2-1,1)];
x = A\b;%(A+1e-8*eye(size(A)))\b;
term1 = 1;
term2 = 1;
term3 = 1;
term4 = 1;
        
for index_i = 1:length(x)
    term1 = term1 + x(index_i)*cos((2*index_i-1)*theta_1);
    term2 = term2 + x(index_i)*cos((2*index_i-1)*theta_2);
    term3 = term3 + x(index_i)*cos((2*index_i-1)*theta_1_hat);
    term4 = term4 + x(index_i)*cos((2*index_i-1)*theta_2_hat);
end
I_h_H_1 = 0.25*term1*term2;
I_h_H_2 = 0.25*term3*term4;
I_h_H_3 = 0.25*term3*term2;
I_h_H_4 = 0.25*term1*term4;

return