clear;
close all;
clc

data =load('exchanger.dat');
%data =load('heating_system.dat');

u =data(:,2)';
y =data(:,3)';
L =size(data,1);

num_train=round(0.75*L);


%%
na =3;
nb =2;
ng =2;
n0 =na +nb +ng -1;
Y = y(na+1:L)';
Y_train=Y(1:num_train);
Y_test=Y(num_train+1:end);
x_hat(1,:) = normrnd(0,0.6,[1,L]);

%x_hat(1,:) = zeros(1,L);
theta(:,1) = 1e-6*ones(n0,1);

for k =2:50
    
    for t=na+1:L
        fs = [ -x_hat(k-1,t-1:-1:t-na) u(t-1:-1:t-nb) ];
        
        for g=2:ng
            fg(1,g-1) = x_hat(k-1,t)^g ;
        end
        fi(:,t-na,k-1) = [fs fg]';
    end
    fi_train =fi(:,1:num_train,k-1);
    fi_test =fi(:,num_train+1:end,k-1);
    
    theta(:,k) = (fi_train*fi_train' +eye(n0))^-1 *fi_train *Y_train;
    
    x_hat(k,na+1:end) = (fi(1:na+nb,:,k-1)' * theta(1:na+nb ,k) )';
    
    y_hat =fi_test'*theta(:,k);
    e(1,k-1)=norm(Y_test-y_hat);
    
    if norm(theta(:,k)-theta(:,k-1))<0.2
        break;
    end
    
end

%%
semilogy(1:k-1,e);grid on;

fitting=100*(1-((norm(y_hat-Y_test))/norm(Y_test-mean(Y_test))))