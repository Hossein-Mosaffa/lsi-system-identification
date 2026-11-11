clc
clear 
close all

%======================================
%======================================
%/box jenkins_recursive
%======================================
%======================================


L=1000;            %data length
sigma_2=0.25;         %noise variance


%----------------------------------
%----------------------------------
%n_a=2
%n_b=2
%n_c=1
%n_d=1
%----------------------------------
%----------------------------------


p0=10^6;      %p0
P0=p0*eye(6);    % % initial value for covariance matrix
teta_hat_0=1/p0*ones(6,1);   %  %  initial value for teta_hat


teta_s=[-0.4;0.44;-1.8;0.24];       %about system
teta_n=[-0.48;0.85];                 %about noise
teta=[teta_s;teta_n];     %true teta


u=normrnd(0,1,[1,L]);            %system input
v=normrnd(0,sigma_2,[1,L]);         %white noise


%   % first and second initial values for collecting the output of system:
% % if t<=0   ,, u(t),v(t)=0
phi_s(:,1)=[0;0;0;0];
x(1)=(phi_s(:,1))'*teta_s;
phi_s(:,2)=[-x(1);0;u(1);0];
x(2)=(phi_s(:,2))'*teta_s;

phi_n(:,1)=[0;0];
w(1)=(phi_n(:,1))'*teta_n+v(1);
phi_n(:,2)=[-w(1);v(1)];
w(2)=(phi_n(:,2))'*teta_n+v(2);

phi(:,1)=[phi_s(:,1);phi_n(:,1)];
phi(:,2)=[phi_s(:,2);phi_n(:,2)];

y(1)=(phi(:,1))'*teta+v(1);
y(2)=(phi(:,2))'*teta+v(2);


% % collecting system output:
for a=3:L
    phi_n(:,a)=[-w(a-1);v(a-1)];
    w(a)=(phi_n(:,a))'*teta_n+v(a);
    
    phi_s(:,a)=[-x(a-1);-x(a-2);u(a-1);u(a-2)];
    x(a)=(phi_s(:,a))'*teta_s;
    
    phi(:,a)=[phi_s(:,a);phi_n(:,a)];
    
    %  % outpuy of system:
    y(a)=(phi(:,a))'*teta+v(a);
end

%   % first initial values:
phi_hat_s(:,1)=[-1/p0;-1/p0;0;0];
phi_hat_n(:,1)=[-1/p0;1/p0];
phi_hat(:,1)=[phi_hat_s(:,1);phi_hat_n(:,1)];

P(:,:,1)=P0-P0*phi_hat(:,1)*(phi_hat(:,1))'*P0/(1+(phi_hat(:,1))'*P0*phi_hat(:,1));

teta_hat(:,1)=teta_hat_0+P(:,:,1)*phi_hat(:,1)*(y(1)-(phi_hat(:,1))'*teta_hat_0);
b=teta_hat(:,1);
teta_hat_s(:,1)=b(1:4,:);
teta_hat_n(:,1)=b(5:6,:);

x_hat(1)=(phi_hat_s(:,1))'*teta_hat_s(:,1);
w_hat(1)=y(1)-x_hat(1);
v_hat(1)=w_hat(1)-(phi_hat_n(:,1))'*teta_hat_n(:,1);



% % second initiel values:
phi_hat_s(:,2)=[-x_hat(1);-1/p0;u(1);0];
phi_hat_n(:,2)=[-w_hat(1);v_hat(1)];
phi_hat(:,2)=[phi_hat_s(:,2);phi_hat_n(:,2)];

P(:,:,2)=P(:,:,1)-P(:,:,1)*phi_hat(:,2)*(phi_hat(:,2))'*P(:,:,1)/(1+(phi_hat(:,2))'*P(:,:,1)*phi_hat(:,2));

teta_hat(:,2)=teta_hat(:,1)+P(:,:,2)*phi_hat(:,2)*(y(2)-(phi_hat(:,2))'*teta_hat(:,1));
c=teta_hat(:,2);
teta_hat_s(:,2)=c(1:4,:);
teta_hat_n(:,2)=c(5:6,:);

x_hat(2)=(phi_hat_s(:,2))'*teta_hat_s(:,2);
w_hat(2)=y(2)-x_hat(2);
v_hat(2)=w_hat(2)-(phi_hat_n(:,2))'*teta_hat_n(:,2);



for t=3:L
    % % information matrix:
    phi_hat_s(:,t)=[-x_hat(t-1);-x_hat(t-2);u(t-1);u(t-2)];
    phi_hat_n(:,t)=[-w_hat(t-1);v_hat(t-1)];
    phi_hat(:,t)=[phi_hat_s(:,t);phi_hat_n(:,t)];
    
    % % covariance matrix:
    P(:,:,t)=P(:,:,t-1)-P(:,:,t-1)*phi_hat(:,t)*(phi_hat(:,t))'*P(:,:,t-1)/(1+(phi_hat(:,t))'*P(:,:,t-1)*phi_hat(:,t));
    
    % % online estimation: 
    teta_hat(:,t)=teta_hat(:,t-1)+P(:,:,t)*phi_hat(:,t)*(y(t)-(phi_hat(:,t))'*teta_hat(:,t-1));
    d=teta_hat(:,t);
    teta_hat_s(:,t)=d(1:4,:);
    teta_hat_n(:,t)=d(5:6,:);
    
    % % inner variables:
    x_hat(t)=(phi_hat_s(:,t))'*teta_hat_s(:,t);
    w_hat(t)=y(t)-x_hat(t);
    v_hat(t)=w_hat(t)-(phi_hat_n(:,t))'*teta_hat_n(:,t);
end
   % % estimate of teta and error:
   TETA_HAT=teta_hat(:,L) 
   error=norm(TETA_HAT-teta)/norm(teta)*100