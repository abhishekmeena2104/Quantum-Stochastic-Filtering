

%  Name - Abhishek Meena
%  Department of Electrical Engineering
%  Indian Iinstitute of Technology Kanpur

clc
clear all
%% intializations of parameters
% Linear case
beta = 0.08;
gamma = 2;
m = 0.10;
tau = 480;
h=1;
delx = 0.1;
delt=0.01;
N = 400; % no of neurons
%%
nu=0;
t=0:delt:50;
[q,epoch]=size(t);

psi=rand(400,epoch);
K = 0.5.*(2*rand(N,1)-1);  % weights
v=zeros(N,epoch); % potential function
y_est=zeros(1,epoch);
%%
signal(1,1:epoch)=2;
y=awgn(signal,6);
psi(:,1) = exp(-((-199:200).^2)./(2*1)).*(1/sqrt(2*pi*1));
psiq(:,1)=psi(:,1)/norm(psi(:,1));
X=-199:1:200;
inte=0;
for k=1:N
    inte=inte+(X(k)*(norm(psiq(k,1))^2)*delx);
end
y_est(1) = inte;
nu=y(1)-y_est(1);
%% training
for i=1:epoch-1
    for ga =1:gamma
        nu=y(i)-y_est(i);
        v(:,i)=-(480*nu*K(:,i));
        for k = 2:N-1
            psi(k,i+1) = psi(k,i)+1j*( -delt*v(k,i)*psi(k,i) + ((delt*(psi(k+1,i)-2*psi(k,i)+psi(k-1,i)))/(2*m*delx*delx) ));
        end
        psi(1,i+1) = psi(1,i)+1j*( -delt*v(k,i)*psi(k,i) + (delt*(psi(2,i)-2*psi(1,i)))/(2*m*delx*delx) );
        psi(N,i+1) = psi(N,i)+1j*( -delt*v(k,i)*psi(N,i) + (delt*(-2*psi(N,i)+psi(N-1,i)))/(2*m*delx*delx) );
        
        psi(:,i+1)=psi(:,i+1)/norm(psi(:,i+1));
        
        X=-199:1:200;
        inte=0;
        for k=1:N
            inte=inte+(X(k)*(norm(psi(k,i+1))^2)*delx);
        end
        y_est(i+1) = inte;
        nu=y(i)-y_est(i);
        for k=1:N
            K(k,i+1) = K(k,i)+beta*nu*(norm(psi(k,i))^2);
        end
    end
        figure(3)
        plot(t,signal,'r',t,y_est,'b');
        ylim([-5 5]);
        legend('DC Signal(6db SNR)','Estimated Signal')
        hold on;
        plot(t,y_est,'b');
           ylim([-5 5]);
        hold off
        pause(0.0001)
    
end
%% RMSE AND PLOT

% figure(3)
% plot(t,signal,'r',t,y_est,'b');
% ylim([-5 5]);
legend('DC Signal(6db SNR)','Estimated Signal')
error=(signal(1,100:10001) - y_est(1,100:10001));    % Errors
error=error.^2;   % Squared Error

RMSE = sqrt(mean(error));  % Root Mean Squared ErroR



figure(4)
plot(-199:200,(smooth(psi(:,1).*conj(psi(:,1)))),'y');
hold on
plot(-199:200,(smooth(psi(:,300).*conj(psi(:,300)))),'r');
hold on
plot(-199:200,(smooth(psi(:,2500).*conj(psi(:,2500)))),'g');
hold on
plot(-199:200,(smooth(psi(:,10000).*conj(psi(:,10000)))),'b');
hold off
title('wave packets at different time instants 0s / 10s/ 50s / 100s/')

        
     
     
     
     
     
     