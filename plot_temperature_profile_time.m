function plot_temperature_profile_time(x,S,T,input1, input2,T_st,TM,UR1,UR2,winx,winy)

M = winx*winy;
for i = 1:M
subplot(winx,winy,i);
plot_wall(S,75,0); 
hold on;
plot_temperature_profile(x,T_st+UR1(:,T(i)).','b--');
plot_temperature_profile(x,T_st+UR2(:,T(i)).','r-');
plot_labels(70,TM,S);

text = ['Time = ',num2str(T(i),'%2d'), '(s)'];
%text =['$q_1^{int}$ = ',num2str(input1(T(i))/1000),' (watt$/$m$^2$)','; ',...
 %   '$q_3^{int}$ = ',num2str(input2(T(i))/1000),' (watt$/$m$^2$)'];
title(text,'interpreter','latex');


text_x = ['$x$(m)'];
text_y = ['$T$($^{\circ}$ C)'];
xlabel(text_x,'interpreter','latex');
ylabel(text_y,'interpreter','latex')
ylim ([0 75])



end
