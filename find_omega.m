function [Phi, Phin, index] = find_omega(Omega,X,S,D,N,M,K,H)

rho(1) = 11.08; 
rho(2) = 2.71; 
rho(3) = 7.4;

Cp(1) = 0.031;
Cp(2) = 0.181;
Cp(3) = 0.054;

[~, G1]=transcendental(Omega,M,S,K,H,D);
% options = optimset('TolX',1e-6);
% Omega = fzero(f,x0,options);
% result = subs(simplify(det(Gamma)),Omega)
% G1 = double(subs(Gamma,omega,Omega));

[U,E,V] = svd(G1);

sol = V(:,end);

A = [sol(1), sol(3), sol(5)];
B = [sol(2), sol(4), sol(6)];

for j = 1:N
    x = X(j);
    [Phi(j), index(j)] = phi_i(x,Omega,A,B,D,S);
end

I1 = find(index==1);
I2 = find(index==2);
I3 = find(index==3);

Nn = rho(1)*Cp(1)*trapz(X(I1),Phi(I1).^2)+rho(2)*Cp(2)*trapz(X(I2),Phi(I2).^2)+rho(3)*Cp(3)*trapz(X(I3),Phi(I3).^2);


% Phin(I1) = (1/sqrt(3))*Phi(I1)/sqrt(trapz(X(I1),Phi(I1).^2));
% Phin(I2) = (1/sqrt(3))*Phi(I2)/sqrt(trapz(X(I2),Phi(I2).^2));
% Phin(I3) = (1/sqrt(3))*Phi(I3)/sqrt(trapz(X(I3),Phi(I3).^2));
% % % 
 Phin(I1) = sqrt(rho(1)*Cp(1)).*Phi(I1)./sqrt(Nn); % Phi/sqrt(trapz(X,Phi.^2));
 Phin(I2) = sqrt(rho(2)*Cp(2)).*Phi(I2)./sqrt(Nn); % Phi/sqrt(trapz(X,Phi.^2));
 Phin(I3) = sqrt(rho(3)*Cp(3)).*Phi(I3)./sqrt(Nn); % Phi/sqrt(trapz(X,Phi.^2));
 
% % %  Phin(I1) = Phi(I1)./sqrt(Nn); % Phi/sqrt(trapz(X,Phi.^2));
% % %  Phin(I2) = Phi(I2)./sqrt(Nn); % Phi/sqrt(trapz(X,Phi.^2));
% % %  Phin(I3) = Phi(I3)./sqrt(Nn); % Phi/sqrt(trapz(X,Phi.^2));

% Phin = Phi/sqrt(Nn);

end

