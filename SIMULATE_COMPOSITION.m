%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%                                   Main Script             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Amritam Das
% Date: 12-02-2017

%% Purpose of the software %%
% Modeling and simulation of 1D heat diffusion in Composite media. For
% composite madia the physical properties are piece-wise constant over the
% space. This is an example of "Spectral Method" where explicit fourier
% bases are used to describe the modal functions. In stead of using
% discretized space (mesh) this method uses basis functions to be defined
% on the whole domain. We are using SOBOLEV SPACE to define new piecewise
% continuous basis functions. 
%% How to use %%
% 1. Run this code if you want to just see an example. 
% 2. If you want to change the physical properties, see the readme.txt file.

%% Setting up
clear;
clc;
close all;

%% set parameters
set_parameters;

%% Set ambient temperature.
% So far we are considering only constant temperature
T0 = 30;
TM = 15;

%% set spatial domain
N = 1000;
R_find = 50;
R_use = 50;
R_POD = 5;
Modal = 20;
X = linspace(S(1),S(end),N); dx = X(2) - X(1);
T = linspace(0,1e-3,200);

%% set input
kcal = 5e4;
input1=kcal*square(5000*pi*T);
input2=-kcal*square(5000*pi*T);

u = [input1; input2];

%% construct the intial condition and steady-state solution

T_0 = zeros(1,N);
T_st = zeros(1,N);

for i = 1:N
    x = X(i);
    T_0(i) = T_i0(x,S);
    T_st(i) = T0*g_i1(x,M,K,H,S) + TM*g_i2(x,M,K,H,S);
end
T_0 = T_st;
plot_wall(S,100,10); 
hold on;
plot_temperature_profile(X,T_0,'k:');
plot_temperature_profile(X,T_st,'g--');
plot_labels(85,TM,S);
axis tight
%matlab2tikz('ic.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%% construct Gamma
% syms omega
% [Gamma,sym] = construct_GAMMA();
% 
% SYM = [S(1), S(2), S(3), S(4),...
%        H(1), H(2), H(3), H(4),...
%        K(1), K(2), K(3),...
%        D(1), D(2), D(3), omega];
%    
% Gamma = subs(Gamma,sym,SYM);

%% find omega_n and function U_i

% % % %  omega_n = [10:0.00001:1000];   % Time Vector
% % % %  for j =1:length(omega_n)
% % % %   y(j) = transcendental(omega_n(j),M,S,K,H,D); 
% % % %  end
% % % %  zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
% % % %  zx = zci(y);  
% % % %  Omega = omega_n(zx);

% F = simplify(det(Gamma));
% f = matlabFunction(F);
% 
% x0 = 0; 
% Omega = [0];
Phi = zeros(R_find,N);
Phin = zeros(R_find,N);


% create progress bar
% text = ['Finding ', num2str(R_find), ' natural modes']; h = waitbar(0,text);

load Omega

for i = 1:R_find
    
[Phi_, Phin_, Index] = find_omega(Omega(i),X,S,D,N,M,K,H);
    

Phi(i,:) = Phi_;
Phin(i,:) = Phin_;
 
 
end

I1 = find(Index==1);
I2 = find(Index==2);
I3 = find(Index==3);

% sort the natural frequencies
% [Omega, ind] = sort(Omega);
% Phi = Phi(ind,:);

% plot first 18 natural frequencies
figure;
plot_natural_modes(Phi,Omega,X,S,3,4,'k-');
figure;
plot_natural_modes(Phin,Omega,X,S,3,4,'k-');

%% check orthogonality
R = R_use;

for i = 1:R
    for j = 1:R
        A(i,j) = rho(1)*Cp(1)*trapz(X(I1),Phi(i,I1).*Phi(j,I1))+rho(2)*Cp(2)*trapz(X(I2),Phi(i,I2).*Phi(j,I2))...
            + rho(3)*Cp(3)*trapz(X(I3),Phi(i,I3).*Phi(j,I3));%trapz(X,Phi(i,:).*Phi(j,:));
        An(i,j) = trapz(X,Phin(i,:).*Phin(j,:)); %trapz(X(I1),Phin(i,I1).*Phin(j,I1))+trapz(X(I2),Phin(i,I2).*Phin(j,I2))...
           %+ trapz(X(I3),Phin(i,I3).*Phin(j,I3)); %trapz(X,Phin(i,:).*Phin(j,:)); %
    end
end
% for i = 1:R
%     for j = 1:R
% asd1(i,j)=sum(Phin(i,:).*Phin(j,:))*(X(end)-X(1))/length(I1);
%     end
% end
figure
subplot(1,2,1);
surf(A); axis tight; view(45,30);
title('$<\phi_{i,n}(x),\phi_{i,m}(x)>$','interpreter','latex')
subplot(1,2,2);
surf(An); axis tight; view(45,30);
title('$<\psi_n(x),psi_m(x)>$','interpreter','latex')

%% fit intial condition
w_0 = T_0 - T_st;

c = zeros(R,N);
cn = zeros(R,N);
C = zeros(R,R);
Cn = zeros(R,R);

for i = 1:R
    for j = 1:R
        C(i,j) = rho(1)*Cp(1)*trapz(X(I1),Phi(i,I1).*Phi(j,I1))+rho(2)*Cp(2)*trapz(X(I2),Phi(i,I2).*Phi(j,I2))...
            + rho(3)*Cp(3)*trapz(X(I3),Phi(i,I3).*Phi(j,I3));%trapz(X,Phi(i,:).*Phi(j,:));
        
        Cn(i,j) =trapz(X,Phin(i,:).*Phin(j,:));%trapz(X(I1),Phin(i,I1).*Phin(j,I1))+trapz(X(I2),Phin(i,I2).*Phin(j,I2))...
            %+ trapz(X(I3),Phin(i,I3).*Phin(j,I3));%trapz(X,Phi(i,:).*Phi(j,:))
    end
        B(i,1) = rho(1)*Cp(1)*trapz(X(I1),w_0(I1).*Phi(i,I1))+rho(2)*Cp(2)*trapz(X(I2),w_0(I2).*Phi(i,I2))...
            + rho(3)*Cp(3)*trapz(X(I3),w_0(I3).*Phi(i,I3));%trapz(X,w_0.*Phi(i,:));
        
        Bn(i,1) =trapz(X,w_0.*Phin(i,:));% %trapz(X(I1),w_0(I1).*Phin(i,I1))+trapz(X(I2),w_0(I2).*Phin(i,I2))...
            %+ trapz(X(I3),w_0(I3).*Phin(i,I3));%trapz(X,w_0.*Phin(i,:));
end
c = linsolve(C,B); 
cn = linsolve(Cn,Bn); 

U = (repmat(c,1,N).*Phi).'*ones(R,1);
Un = (repmat(cn,1,N).*Phin).'*ones(R,1);

%% fit reference data
w_f = T_0 + T_st;

d = zeros(R,N);
dn = zeros(R,N);
dd = zeros(R,R);
ddn = zeros(R,R);

for i = 1:R
    for j = 1:R
        dd(i,j) = rho(1)*Cp(1)*trapz(X(I1),Phi(i,I1).*Phi(j,I1))+rho(2)*Cp(2)*trapz(X(I2),Phi(i,I2).*Phi(j,I2))...
            + rho(3)*Cp(3)*trapz(X(I3),Phi(i,I3).*Phi(j,I3));%trapz(X,Phi(i,:).*Phi(j,:));
        
        ddn(i,j) =trapz(X,Phin(i,:).*Phin(j,:));%trapz(X(I1),Phin(i,I1).*Phin(j,I1))+trapz(X(I2),Phin(i,I2).*Phin(j,I2))...
            %+ trapz(X(I3),Phin(i,I3).*Phin(j,I3));%trapz(X,Phi(i,:).*Phi(j,:))
    end
        Bb(i,1) = rho(1)*Cp(1)*trapz(X(I1),w_f(I1).*Phi(i,I1))+rho(2)*Cp(2)*trapz(X(I2),w_f(I2).*Phi(i,I2))...
            + rho(3)*Cp(3)*trapz(X(I3),w_f(I3).*Phi(i,I3));%trapz(X,w_0.*Phi(i,:));
        
        Bbn(i,1) =trapz(X,w_0.*Phin(i,:));% %trapz(X(I1),w_0(I1).*Phin(i,I1))+trapz(X(I2),w_0(I2).*Phin(i,I2))...
            %+ trapz(X(I3),w_0(I3).*Phin(i,I3));%trapz(X,w_0.*Phin(i,:));
end
d = linsolve(dd,Bb); 
dn = linsolve(ddn,Bbn); 

Uu = (repmat(d,1,N).*Phi).'*ones(R,1);
Uun = (repmat(dn,1,N).*Phin).'*ones(R,1);

figure(1);
plot_temperature_profile(X,Uun.' + T_st,'b-');
hold on;
plot_temperature_profile(X,Un.' - T_st,'k-');

%% construct state space

l1 = ones(size(X)); l1(I2) = 0; l1(I3) = 0;
l2 = ones(size(X)); l2(I1) = 0; l2(I3) = 0;
l3 = ones(size(X)); l3(I2) = 0; l3(I1) = 0;
I{1} = I1; I{2} = I2; I{3} = I3; 

A = []; B = []; Norm = [];
for i = 1:R
    Norm (i,i) = rho(1)*Cp(1)*trapz(X(I1),Phi(i,I1).^2)+rho(2)*Cp(2)*trapz(X(I2),Phi(i,I2).^2)...
        +rho(3)*Cp(3)*trapz(X(I3),Phi(i,I3).^2);

    A(i,i) = -Omega(i)^2;
%%%%    B(i,1) = rho(1)*Cp(1)*trapz(X(I1),l1(I1).*Phi(i,I1))+rho(2)*Cp(2)*trapz(X(I2),l1(I2).*Phi(i,I2))...
%%%%             + rho(3)*Cp(3)*trapz(X(I3),l1(I3).*Phi(i,I3));% trapz(X,l1.*Phin(i,:));
%%%%    B(i,2) = rho(1)*Cp(1)*trapz(X(I1),l3(I1).*Phi(i,I1))+rho(2)*Cp(2)*trapz(X(I2),l3(I2).*Phi(i,I2))...
%%%%             + rho(3)*Cp(3)*trapz(X(I3),l3(I3).*Phi(i,I3));%trapz(X,l3.*Phin(i,:));
    B(i,1) = (D(1)^2/K(1))*trapz(X,l1.*Phin(i,:));
    B(i,2) = (D(3)^2/K(3))*trapz(X,l3.*Phin(i,:));

end

C = eye(R,R);

%%%% sys= ss(A,inv(Norm)*B,C,zeros(R,2)); 
sys= ss(A,B,C,zeros(R,2)); 

%% simulate state space


%%%% [y,t,q] = lsim(sys,u,T,c);
[y,t,q] = lsim(sys,u,T,cn);

%% Temperature Alone
for i = 1:length(t)
   UR1(:,i) = (repmat(y(i,:).',1,N).*Phin).'*ones(R,1);
   UR2(:,i) = (repmat(y(i,1:Modal).',1,N).*Phin(1:Modal,:)).'*ones(Modal,1);%3 %20
end
%% Temperature + Moisture [Same Basis]
% % % % for i = 1:length(t)
% % % %    UR1(:,i) = (repmat(y(i,:).',1,N).*Phin).'*ones(R,1);
% % % %    UR2(:,i) = (repmat(y(i,1:Modal).',1,N).*Phin(1:Modal,:)).'*ones(Modal,1);%3 %20
% % % % end
% % % % for i = 1:length(t)
% % % %    US1(:,i) = (repmat(y(i,:).',1,N).*Phin).'*ones(R,1);
% % % %    US2(:,i) = (repmat(y(i,1:Modal).',1,N).*Phin(1:Modal,:)).'*ones(Modal,1);%3 %20
% % % % end
%% Temperature + Moisture [diff Basis, same coeff]
% % % % for i = 1:length(t)
% % % %    UR1(:,i) = (repmat(y(i,:).',1,N).*Phin).'*ones(R,1);
% % % %    UR2(:,i) = (repmat(y(i,1:Modal).',1,N).*Phin(1:Modal,:)).'*ones(Modal,1);%3 %20% Modal
% % % %    %UR1(:,i) = (repmat(y(i,1:Modal).',1,N).*Phin(1:Modal,:)).'*ones(Modal,1);%3 %20
% % % %     UR3(:,i) = (repmat(y(i,5:10).',1,N).*Phin(5:10,:)).'*ones(6,1);%3 %20
% % % % end

%disp('press any key to start the simulation video')
%pause;

% r=10;
% % video name with current time stamp
% VidName = ['r' num2str(r) '_2d_' datestr(now,'ddmmyyHHMMSS') '.mp4'];
% writerObj = VideoWriter(['C:\Users\AmDas\Dropbox\PhD\MATLAB\MODRED' VidName],'MPEG-4');
% % number of frames per second
% writerObj.FrameRate = 5;
% writerObj.Quality = 100;
% open(writerObj);


% set(gcf, 'Position', [100, 100, 1100, 600]); % give figure a certain size (can comment it out if u want)
% set(gcf,'color','w');	% make background white	
% % figure()
% % plot_temperature_profile_time(X,S,[1:10:31],input1, input2, T_st,TM,UR1, UR2,2,2); %[50:50:length(T)]
% % %% For savepdf
% % set(gcf,'Renderer','painters');
% % figureHandle = gcf;
% % set(findall(figureHandle,'type','text'),'fontSize',12,'interpreter','latex');
% % set(findall(figureHandle,'type','axes'),'Color','w','fontSize',12);
% % set(findall(figureHandle,'type','axes'),'TickLabelInterpreter','latex');
% % 
% % %axis tight
% % %savePDF(figureHandle.Position([3]),figureHandle.Position([4]))
% % 
% % %%
% % 
% % figure()
% % plot_temperature_profile_diff(X,S,[1:10:31],input1, input2, T_st,TM,UR1, UR2,2,2); %[50:50:length(T)]
% % %% For savepdf
% % set(gcf,'Renderer','painters');
% % figureHandle = gcf;
% % set(findall(figureHandle,'type','text'),'fontSize',12,'interpreter','latex');
% % set(findall(figureHandle,'type','axes'),'Color','w','fontSize',12);
% % set(findall(figureHandle,'type','axes'),'TickLabelInterpreter','latex');
% % %axis tight
% % savePDF(figureHandle.Position([3]),figureHandle.Position([4]))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure(8);
% % % for i = 1:length(T)
% % %     cla
% % %     plot_wall(S,100,10); 
% % %     hold on;
% % %    plot_temperature_profile(X,T_0,'k:');
% % %    plot_temperature_profile(X,T_st,'g:');
% % %     %% Temperature Alone
% % %    plot_temperature_profile(X,T_st+UR1(:,i).','b--');
% % %    plot_temperature_profile(X,T_st+UR2(:,i).','r-');
% % %    %% Temperature + Moisture [Same Basis]
% % % %    plot_temperature_profile(X,T_st+UR1(:,i).','b--');
% % % %    plot_temperature_profile(X,T_st+UR2(:,i).','r-')
% % % %    plot_temperature_profile(X,T_st+US1(:,i).'-2*ones(1,1000),'g--');
% % % %    plot_temperature_profile(X,T_st+US2(:,i).'-2*ones(1,1000),'k-');
% % % % % % %      plot_temperature_profile(X,0.78*(T_st+UR1(:,i).'),'k-');
% % % % % % %      plot_temperature_profile(X,0.78*(T_st+UR2(:,i).'),'c-');
% % % % % % %      plot_temperature_profile(X,T_st+UR3(:,i).','r-');
% % %     %%
% % %     %plot_temperature_profile(X,T_st+UR1(:,i).'- 0.25*UR2(:,i).','b-');
% % %     text_x = ['$x$(m)'];
% % %     text_y = ['$T$(degree C)',' $u$(mg/m$^3$)'];
% % %    % title(text,'interpreter','latex');
% % %     %text = ['FO order: ', num2str(R), '$\;$(blue); $\;$', 'FO order: ', num2str(Modal), '$\;$(red)'];
% % %     xlabel(text_x,'interpreter','latex');
% % %     ylabel(text_y,'interpreter','latex')
% % %     
% % %     plot_labels(85,TM,S);
% % %     axis([S(1) S(end) 10 100]);
% % %     pause(1/60);
% % %     %pause;
% % % % %     drawnow
% % % % % 
% % % % %     set(gcf,'Renderer','painters');
% % % % %     frame = getframe(gcf);
% % % % %  
% % % % %     writeVideo(writerObj,frame);
% % %     %Image(i) = getframe;
% % % end
% % % 
% % % close(writerObj);
% % % 
% % % 
% % % %matlab2tikz('sim2.tikz', 'height', '\figureheight', 'width', '\figurewidth');
% % %   %matlab2tikz('sim1.tikz', 'height', '\figureheight', 'width', '\figurewidth');
% % % 
% % % 
% % % 
% % % 
% % % for i=1:numel(time)
% % % 	% plot ur stuff and put legend etc.
% % % 	plot(x,y(t))
% % % 	
% % % 	drawnow
% % % 
% % %     	set(gcf,'Renderer','painters');
% % %     	frame = getframe(gcf);
% % %  
% % %     	writeVideo(writerObj,frame);
% % % end
% % % 
% % % close(writerObj);
% % % 
% % % 
% % % %% POD
% % % 
% % % W = UR1*UR1.';
% % % 
% % % [U,E,V] = svd(W);
% % % 
% % % Phi_POD = (1/(sqrt(dx))*U(:,1:R_POD).');
% % % %Phi_POD = U(:,1:R_POD).';
% % % A = [];
% % % for i = 1:R_POD
% % %     for j = 1:R_POD;
% % %         A(i,j) = trapz(X,Phi_POD(i,:).*Phi_POD(j,:));
% % %     end
% % %     
% % %     Phi_POD(i,:) = Phi_POD(i,:)/sqrt(trapz(X,Phi_POD(i,:).^2));
% % % end
% % % 
% % % figure(5)
% % % subplot(1,2,1)
% % % surf(A); axis([0 R_POD 0 R_POD -0.1 1.1]); axis tight; view(45,30);
% % % title('$<\phi_{pod,n}(x),\phi_{pod,k}(x)>$','interpreter','latex')
% % % subplot(1,2,2)
% % % plot_wall(S,2,-2); hold on; 
% % % plot(X,Phi_POD);
% % % title('Natural modes for the POD','interpreter','latex')
% % % 
% % % %%
% % % 
% % % for i = 1:R_POD
% % %     for j = 1:R_POD
% % %         Cp(i,j) = (trapz(X,Phi_POD(i,:).*Phi_POD(j,:)));
% % %     end
% % %     
% % %     Bp(i,1) = trapz(X,w_0.*Phi_POD(i,:));
% % % end
% % % 
% % % cp = linsolve(Cp.',Bp); 
% % % 
% % % Phi_POD_ddot1=[]; Phi_POD_ddot2=[]; Phi_POD_ddot3=[]; Phi_POD_ddot=[];
% % % for i = 1:R_POD
% % %     %Phi_POD_ddot = [Phi_POD_ddot; gradient(gradient(Phi_POD(i,:),dx),dx)];
% % %     Phi_POD_ddot1 = [Phi_POD_ddot1; gradient(gradient(Phi_POD(i,I1),dx),dx)];
% % %     Phi_POD_ddot2 = [Phi_POD_ddot2; gradient(gradient(Phi_POD(i,I2),dx),dx)];
% % %     Phi_POD_ddot3 = [Phi_POD_ddot3; gradient(gradient(Phi_POD(i,I3),dx),dx)];
% % %     %Phi_POD_ddot1 = [Phi_POD_ddot1; num_diffxx(Phi_POD(i,I1),dx)];
% % %     %Phi_POD_ddot2 = [Phi_POD_ddot2; num_diffxx(Phi_POD(i,I2),dx)];
% % %     %Phi_POD_ddot3 = [Phi_POD_ddot3; num_diffxx(Phi_POD(i,I3),dx)];
% % % end
% % % 
% % % Phi_POD_ddot = [Phi_POD_ddot1,Phi_POD_ddot2,Phi_POD_ddot3];
% % % 
% % % A_POD = []; B = [];Norm_POD = [];
% % % for i = 1:R_POD
% % %     for j = 1:R_POD
% % %         Norm_POD(i,j) = trapz(X,Phi_POD(i,:).*Phi_POD(j,:));
% % %         A_POD(i,j) = (D(1)^2)*trapz(X(I1),Phi_POD_ddot(i,I1).*Phi_POD(j,I1)) + (D(2)^2)*trapz(X(I2),Phi_POD_ddot(i,I2).*Phi_POD(j,I2)) + (D(3)^2)*trapz(X(I3),Phi_POD_ddot(i,I3).*Phi_POD(j,I3));
% % %         %A(i,j) = (D(1)^2)*trapz(X,Phi_POD_ddot(i,:).*Phi_POD(j,:)) + (D(2)^2)*trapz(X,Phi_POD_ddot(i,:).*Phi_POD(j,:)) + (D(3)^2)*trapz(X,Phi_POD_ddot(i,:).*Phi_POD(j,:));
% % %     end
% % %     
% % %     B(i,1) = (D(1)^2/K(1))*trapz(X,l1.*Phi_POD(i,:));
% % %     B(i,2) = (D(3)^2/K(3))*trapz(X,l3.*Phi_POD(i,:));
% % % end
% % % 
% % % C = eye(R_POD,R_POD);
% % % D = zeros(R_POD,2);
% % % 
% % % %%%% sys = ss(A_POD.',(Norm(1:R_POD,1:R_POD))*B,C,D);
% % % sys = ss(A_POD.',B,C,D);
% % % 
% % % [y,t,q] = lsim(sys,u,T,cp);
% % % 
% % % for i = 1:length(t)
% % %    UPOD(:,i) = (repmat(y(i,:).',1,N).*Phi_POD).'*ones(R_POD,1);
% % % end
% % % 
% % % disp('press any key to start the simulation video for POD')
% % % pause;
% % % figure(1);
% % % for i = 1:length(T)
% % %     cla
% % %     plot_wall(S,100,10); 
% % %     hold on;
% % %     plot_temperature_profile(X,T_0,'k:');
% % %     plot_temperature_profile(X,T_st,'k--');
% % %     plot_temperature_profile(X,T_st+UR1(:,i).','b-');
% % %     plot_temperature_profile(X,T_st+UR2(:,i).','r-');
% % %     plot_temperature_profile(X,T_st+UPOD(:,i).','g-');
% % %     
% % %     text = ['$q_1$ = ',num2str(input1(i)/1000),' (kcal)','; ','$q_2$ = ',num2str(input2(i)/1000),' (kcal)'];
% % %     %'time = ',num2str(1000*T(i),'%2.3f'),' (ms)','; ',
% % %     title(text,'interpreter','latex');
% % %     text = ['FO order: ', num2str(R), '$\;$(blue); $\;$', 'POD order: ', num2str(R_POD), '$\;$(green); $\;$', 'FO order: ', num2str(Modal), '$\;$(red)'];
% % %     xlabel(text,'interpreter','latex')
% % %     plot_labels(85,TM,S);
% % %     axis([S(1) S(end) 10 100]);
% % %     pause(1/60);
% % %    % pause;
% % % %     drawnow
% % % 
% % % end
% % % % matlab2tikz('reference.tikz', 'height', '\figureheight', 'width', '\figurewidth');

