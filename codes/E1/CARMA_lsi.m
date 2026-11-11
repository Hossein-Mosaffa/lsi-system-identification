clc
clear 
close all

%_____________________________________________
% CARMA\LSI %
%_____________________________________________


L=3000;    % data length
sigma_2=0.25;    % noise variance
e=0.001;   % error tolerance

u1=normrnd(0,1,[1,L]);
u2=normrnd(0,1,[1,L]);
u=[u1;u2];       %system input

v1=normrnd(0,sigma_2,[1,L]);
v2=normrnd(0,sigma_2,[1,L]);
v=[v1;v2];       %white noise

A1=[0.6 0.5;-0.8 1];                       %n_a=1
B1=[1.5 -0.4;-0.5 1.1];                    %n_b=1
D1=[0.2 -0.1;-0.1 0.6];                    %n_c=1
teta=[A1';B1';D1'];      %true teta

y(:,1)=v(:,1);      %y(1)=v(1)    if t<=0  : u(t),y(t),v(t)=0
Y(:,1)=y(:,1);

for a=2:L
    phi(:,a)=[-y(:,a-1);u(:,a-1);v(:,a-1)];
    y(:,a)=teta'*phi(:,a)+v(:,a);
    Y(:,a)=y(:,a);    %stacked output matrix
end

% initial value:
vhat1=rand([1,L]);
vhat2=rand([1,L]);
vhat=[vhat1;vhat2];       %noise estimate: v_hat(0)=random vector

% % iterative alghorithm:
for k=1:20
    phi_hat(:,1)=zeros(6,1);
    PHI_hat(:,1)=phi_hat(:,1);
    for b=2:L
        phi_hat(:,b)=[-y(:,b-1);u(:,b-1);vhat(:,b-1)];
        PHI_hat(:,b)=phi_hat(:,b);
    end
    
    % estimation of system parameters:
    teta_hat(:,:,k)=((PHI_hat*PHI_hat')^(-1))*PHI_hat*Y'
    
    vhat=y-teta_hat(:,:,k)'*phi_hat;
    
    % % estimation error:
    error(k)=norm(teta_hat(:,:,k)-teta)/norm(teta)*100
    
    %-------------------------------------------------------------
    %  % stop condition:
    if (k>=2)   &&   (norm(teta_hat(:,:,k)-teta_hat(:,:,k-1))<=e)   
    break;
    else
    end
    %-------------------------------------------------------------
end