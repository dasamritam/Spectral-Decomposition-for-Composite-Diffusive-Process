function G = g_i2(x,no_layer,K,H,thick)
h = H;
G = NaN;

P_2 = zeros(2*no_layer,2*no_layer);
Y_2 = zeros(2*no_layer,1);

P_2(1,1) = K(1) - h(1)*thick(1); % filling up P_1 [Row 1]
P_2(1,2) = -h(1);

P_2(2*no_layer,2*no_layer-1) = K(no_layer) + h(no_layer+1)*thick(no_layer+1); % filling up P_1 [Row 2M]
P_2(2*no_layer,2*no_layer) = h(no_layer+1);

for index = 1:1:no_layer-1 % filling rows from 2 till 2m-1
    P_2(2*index,2*index-1) = (K(index)/h(index+1)) + thick(index+1);
    P_2(2*index,2*index) = 1;
    P_2(2*index,2*index+1) = -thick(index+1);
    P_2(2*index,2*index+2) = -1; 
    
    P_2(2*index+1,2*index-1) = K(index)/K(index+1);
    P_2(2*index+1,2*index+1) = -1;
end


Y_2(2*no_layer) = h(no_layer+1);

steady_bottom = inv(P_2)*Y_2;

for layer = 1:no_layer
     if (thick(layer) <= x) && (x <=  thick(layer+1))
          G = steady_bottom(2*layer-1)*x+steady_bottom(2*layer);
     end 
end
  
end

