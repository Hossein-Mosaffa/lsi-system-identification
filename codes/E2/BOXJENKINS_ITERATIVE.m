clc
clear 
close all

%======================================
%======================================
%/box jenkins_iterative
%======================================
%======================================


% %data length and noise variance:
L=2000;
sigma_2=1;

% error tolerance:
e=0.001;

%----------------------------------
%----------------------------------
%n_a=2
%n_b=2
%n_c=1
%n_d=1
%----------------------------------
%----------------------------------

p0=10^6;

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

Y=zeros(L,1);            %stacked output matrix
Y(L)=y(1);
Y(L-1)=y(2);

% % collecting system output:
for a=3:L
    phi_n(:,a)=[-w(a-1);v(a-1)];
    w(a)=(phi_n(:,a))'*teta_n+v(a);
    
    phi_s(:,a)=[-x(a-1);-x(a-2);u(a-1);u(a-2)];
    x(a)=(phi_s(:,a))'*teta_s;
    
    phi(:,a)=[phi_s(:,a);phi_n(:,a)];
    
    %  % outpuy of system:
    y(a)=(phi(:,a))'*teta+v(a);
    Y(L+1-a)=y(a);            
end

% initial values:
x_hat=rand([1,L]);
w_hat=rand([1,L]);
v_hat=rand([1,L]);

for k=1:15       %k=iteration variable
    phi_hat_s(:,1)=[-1/p0;-1/p0;0;0];
    phi_hat_n(:,1)=[-1/p0;1/p0];
    phi_hat(:,1)=[phi_hat_s(:,1);phi_hat_n(:,1)];
    PHI_hat(L,:)=(phi_hat(:,1))';
    
    phi_hat_s(:,2)=[-x_hat(1);-1/p0;u(1);0];
    phi_hat_n(:,2)=[-w_hat(1);v_hat(1)];
    phi_hat(:,2)=[phi_hat_s(:,2);phi_hat_n(:,2)];
    PHI_hat(L-1,:)=(phi_hat(:,2))';
    
    for t=3:L
        phi_hat_s(:,t)=[-x_hat(t-1);-x_hat(t-2);u(t-1);u(t-2)];
        phi_hat_n(:,t)=[-w_hat(t-1);v_hat(t-1)];
        phi_hat(:,t)=[phi_hat_s(:,t);phi_hat_n(:,t)];
        PHI_hat(L+1-t,:)=(phi_hat(:,t))';
    end
    
    teta_hat(:,:,k)=(((PHI_hat)'*PHI_hat)^(-1))*(PHI_hat)'*Y
    z=teta_hat(:,:,k);
    teta_hat_s(:,:,k)=z(1:4,:);
    teta_hat_n(:,:,k)=z(5:6,:);
    
    for g=1:L
        x_hat(g)=(phi_hat_s(:,g))'*teta_hat_s(:,:,k);
        w_hat(g)=y(g)-x_hat(g);
        v_hat(g)=y(g)-(phi_hat(:,g))'*teta_hat(:,:,k);
    end
    
    % % estimation error:
    error(k)=norm(teta_hat(:,:,k)-teta)/norm(teta)*100
    
    % %stp condition:
    %-------------------------------------------
    if (k>=2)   &&   (norm(teta_hat(:,:,k)-teta_hat(:,:,k-1))<=e)
        break;
    else
    end
    %-------------------------------------------
end

plot(1:k,error,'.-k');
legend('estimation error');
xlabel('iteration variable');
ylabel('error%');