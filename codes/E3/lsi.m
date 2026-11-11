clear;
close all;
clc

L=3000;

x(1)=0;
x(2)=0;
f(1)=0;
f(2)=0;

u=normrnd(0,1,[1,L]);

for i=3:L
    x(i) =0.18*x(i-1) -0.44*x(i-2) +0.58*u(i-1) +0.41*u(i-2);
end

v=normrnd(0,0.2,[1,L]);

for i=3:L
    f(i) = x(i) -0.45*x(i)^2 +0.25*x(i)^3;
end

y =f +v;

%%
na =2;
nb =2;
ng =3;
n0 =na +nb +ng -1;
Y = y(3:L)';
% x_hat(1,:) = normrnd(0,0.6,[1,L]);
% x2=x_hat;
x_hat(1,:) = zeros(1,L);
theta(:,1) = 1e-6*ones(n0,1);

for k =2:50
    
    for t=3:L
        fs = [ -x_hat(k-1,t-1:-1:t-na) u(t-1:-1:t-nb) ];
        
        for g=2:ng
            fg(1,g-1) = x_hat(k-1,t)^g;
        end
        fi(:,t-na,k-1) = [fs fg]';
    end
    
    theta(:,k) = (fi(:,:,k-1)*fi(:,:,k-1)' +0.1*eye(n0))^-1 *fi(:,:,k-1) *Y;
    
    x_hat(k,na+1:end) = (fi(1:na+nb,:,k-1)' * theta(1:na+nb ,k) )';
    
    if norm(theta(:,k)-theta(:,k-1))<0.0001
        break;
    end
    
end

t=[-0.18 0.44 0.58 0.41 -0.45 0.25]';
delta =(sqrt(sum((theta-t).^2))/norm(t));

k
theta(:,k)
t

semilogy(1:k,delta);