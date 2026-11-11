clear;
close all;
clc

L=1000;

x(1)=0;
x(2)=0;
f(1)=0;
f(2)=0;

u=normrnd(0,1,[1,L]);

for i=3:L
    x(i) =0.18*x(i-1) -0.44*x(i-2) +0.58*u(i-1) +0.41*u(i-2);
end

v=normrnd(0,0.1,[1,L]);

for i=3:L
    f(i) = x(i) -0.45*x(i)^2 +0.25*x(i)^3;
end

y =f +v;

%%
na =2;
nb =2;
ng =3;
n0 =na +nb +ng -1;
Y = y(na+1:L)';
Y_train=Y(1:700);
Y_test=Y(701:end);
% x_hat(1,:) = normrnd(0,0.6,[1,L]);
% x2=x_hat;
x_hat(1,:) = zeros(1,L);
theta(:,1) = 1e-6*ones(n0,1);

for k =2:100
    
    for t=na+1:L
        fs = [ -x_hat(k-1,t-1:-1:t-na) u(t-1:-1:t-nb) ];
        
        for g=2:ng
            fg(1,g-1) = x_hat(k-1,t)^g;
        end
        fi(:,t-na,k-1) = [fs fg]';
    end
    fi_train =fi(:,1:700,k-1);
    fi_test =fi(:,701:end,k-1);
    
    theta(:,k) = (fi_train*fi_train' +0.1*eye(n0))^-1 *fi_train *Y_train;
    
    x_hat(k,na+1:end) = (fi(1:na+nb,:,k-1)' * theta(1:na+nb ,k) )';
    
    y_hat =fi_test'*theta(:,k);
    e(1,k-1)=norm(Y_test-y_hat);
    
    if norm(theta(:,k)-theta(:,k-1))<0.001
        break;
    end
    
end

%%
semilogy(1:k-1,e);grid on;

fitting=100*(1-((norm(y_hat-Y_test))/norm(Y_test-mean(Y_test))))