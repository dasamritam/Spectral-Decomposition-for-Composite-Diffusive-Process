function [d_Gamma, Gamma_n] = transcendental(omega_n,no_layer,thick,K,h,D)
%%% Solving the boundary conditions for the spatial frequency response U_i(x,t). 
%%% U_i(x,t) = a_in cos(w_n x/D_i) + b_in sin(w_n x/D_i).
%%% Can be written as Gamma_n (w_n) c_n = 0. c_n = col(a_in,...b_in).
%%% w_n can be obtianed from solving det (Gamma_n(w_n)) = 0. n = 1,2,...N.
%%% c_n is obtained from the nullspace of the re-constructed matrix Gamma_n(w_n).

%___________ Construct Gamma_n(w_n) as function of omega_n ________________

Gamma_n = zeros(2*no_layer,2*no_layer);% Initialize a matrix filled with zeros

Gamma_n(1,1) = -(omega_n/D(1))*sin(omega_n*thick(1)/D(1))...
               - (h(1)/K(1))*cos(omega_n*thick(1)/D(1)); % filling up Gamma_n [Row 1]
Gamma_n(1,2) = (omega_n/D(1))*cos(omega_n*thick(1)/D(1))...
               - (h(1)/K(1))*sin(omega_n*thick(1)/D(1));

Gamma_n(2*no_layer,2*no_layer-1) = -(omega_n/D(no_layer))*sin(omega_n*thick(no_layer+1)/D(no_layer))...
               + (h(no_layer+1)/K(no_layer))*cos(omega_n*thick(no_layer+1)/D(no_layer)); % filling up Gamma_n [Row 2M]
Gamma_n(2*no_layer,2*no_layer) = (omega_n/D(no_layer))*cos(omega_n*thick(no_layer+1)/D(no_layer))...
               + (h(no_layer+1)/K(no_layer))*sin(omega_n*thick(no_layer+1)/D(no_layer));

for i = 1:1:no_layer-1 % filling rows
    Gamma_n(2*i,2*i-1) = sqrt(K(i)/(D(i)*D(i)))*(-(K(i)/h(i+1))*(omega_n/D(i))*sin(omega_n*thick(i+1)/D(i))...
               + cos(omega_n*thick(i+1)/D(i)));
    Gamma_n(2*i,2*i) =  sqrt(K(i)/(D(i)*D(i)))*((K(i)/h(i+1))*(omega_n/D(i))*cos(omega_n*thick(i+1)/D(i))...
               + sin(omega_n*thick(i+1)/D(i)));       
    Gamma_n(2*i,2*i+1) = sqrt(K(i+1)/(D(i+1)*D(i+1)))*(-cos(omega_n*thick(i+1)/D(i+1)));    
    Gamma_n(2*i,2*i+2) = sqrt(K(i+1)/(D(i+1)*D(i+1)))*(-sin(omega_n*thick(i+1)/D(i+1)));
    
    Gamma_n(2*i+1,2*i-1) = sqrt(K(i)/(D(i)*D(i)))*(-(K(i)/K(i+1))*(omega_n/D(i))*sin(omega_n*thick(i+1)/D(i)));
    Gamma_n(2*i+1,2*i) =  sqrt(K(i)/(D(i)*D(i)))*((K(i)/K(i+1))*(omega_n/D(i))*cos(omega_n*thick(i+1)/D(i)));       
    Gamma_n(2*i+1,2*i+1) = sqrt(K(i+1)/(D(i+1)*D(i+1)))*((omega_n/D(i+1))*sin(omega_n*thick(i+1)/D(i+1)));    
    Gamma_n(2*i+1,2*i+2) = sqrt(K(i+1)/(D(i+1)*D(i+1)))*(-(omega_n/D(i+1))*cos(omega_n*thick(i+1)/D(i+1)));
end
d_Gamma = det(Gamma_n);