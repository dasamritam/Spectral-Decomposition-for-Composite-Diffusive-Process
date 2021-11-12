function G = g_i1(x,no_layer,K,H,thick)
h = H;
G = NaN;

% P_1 is a matrix of [2M X 2M]. X_1 is unknown as 
% X_1:= [A_11;B_11;A_21;B_21;...;A_M1;B_M1].

P_1 = zeros(2*no_layer,2*no_layer); % Initialize a matrix filled with zeros
Y_1 = zeros(2*no_layer,1);

P_1(1,1) = K(1) - h(1)*thick(1); % filling up P_1 [Row 1]
P_1(1,2) = -h(1);

P_1(2*no_layer,2*no_layer-1) = K(no_layer) + h(no_layer+1)*thick(no_layer+1); % filling up P_1 [Row 2M]
P_1(2*no_layer,2*no_layer) = h(no_layer+1);

for index = 1:1:no_layer-1 % filling rows from 2 till 2m-1
    P_1(2*index,2*index-1) = (K(index)/h(index+1)) + thick(index+1);
    P_1(2*index,2*index) = 1;
    P_1(2*index,2*index+1) = -thick(index+1);
    P_1(2*index,2*index+2) = -1; 
    
    P_1(2*index+1,2*index-1) = K(index)/K(index+1);
    P_1(2*index+1,2*index+1) = -1;
end

Y_1(1) = -h(1);

steady_top = inv(P_1)*Y_1;

for layer = 1:no_layer
     if (thick(layer) <= x) && (x <=  thick(layer+1))
          G = steady_top(2*layer-1)*x+steady_top(2*layer);
     end 
end
  
end

