close all; clear all; clc
rng('shuffle')

x0 = [0.9 0 0 0];
dt = 0.1;
alpha = 1;
alpha = alpha*ones(1,4);
beta = 1;
nu = 1;
eps = 0.001;

neighborhood = 0.8;
neighborhood1 = 0.5;

t1 = [];
t2 = [];
t3 = [];
t4 = [];
x = x0;

Nt = 200;
decision = [2 3]; %which nodes do we want the decision split to be between
rho = decision_rho(alpha,beta,nu,decision); %this works!
rho

figs = 0; %change to 1 if you want plots

plotcolors = [1 0 0; 0 0 0; 0 0 1; 0 1 0];
plotcolors = [1 0 0; 0 0 0; 0 0 1; 1 0 1];


for i = 1:Nt-1
    % noise(i,:) = eps*randn(1,size(rho,1)); %gaussian noise
    noise(i,:) = eps*sqrt(dt)*randn(1,size(rho,1)); %gaussian noise

    dx = x(i,:).*(alpha - x(i,:)*rho)*dt + noise(i,:); %with dt
    %dx needs to be bound? something's up with Matlab's randn
    x(i+1,:) = max(min(x(i,:) + dx, 1), 0.0005);

        if any(x(i+1,1) >= neighborhood)
            t1 = [t1 i];
        end
        if any(x(i+1,2) >= neighborhood)
            t2 = [t2 i];
        end
        if any(x(i+1,3) >= neighborhood)
            t3 = [t3 i];
        end
        if any(x(i+1,4) >= neighborhood)
            t4 = [t4 i];
        end

end


figure(2)
for i = 1:4
    hold on
    plot((1:Nt)*dt,x(:,i),'.')
end
colororder(plotcolors)
if any(t1)
    xline(t1(end)*dt,':k',{'$\tau_t$ start'},'Interpreter','latex',...
        'Fontsize',14,'Linewidth',2,'labelverticalalignment','middle')
end
if any(t2)
    xline(t2(1)*dt,':k',{'$\tau_t$ end'},'Interpreter','latex',...
        'Fontsize',14,'Linewidth',2,'labelverticalalignment','middle')
end
if any(t3)
    xline(t3(1)*dt,':k',{'$\tau_t$ end'},'Interpreter','latex',...
        'Fontsize',14,'Linewidth',2,'labelverticalalignment','middle')
end
yline(neighborhood,'--',{'Decision';'Threshold'},'Color',[.5 .5 .5],...
    'Linewidth',2,'labelverticalalignment','middle','Fontsize',12)
hold off
title('Node Activations & Transition Time','Fontsize',16)
ylabel('State Activation, x', 'Fontsize',14)
xlabel('Time (s)','Fontsize',14)


%% DECISION RHO
function rho = decision_rho(alpha,beta,nu,dec)

% Check dimension of parameters
len = [length(alpha) length(beta) length(nu)];
n = max(len);
if n == 1
    n = 3;
elseif any(len ~= 1 & len ~= n)
    error('decision_rho:ParameterDimensionMismatch',...
         ['The input parameters, Alpha, Beta, and Nu, must be scalars '...
          'or equal length vectors.']);
end

% Expand scalars
z = ones(n,1);
alpha = alpha(:).*z;
beta = beta(:).*z;
nu = nu(:).*z;

%Check decision location, and identify row prior
if dec(1) == 1
    dec_pre = n+1;
else
    dec_pre = dec(1)-1;
end

loc = mod(dec,n);

%Rho construction from Horchler 2015
rho = (alpha./beta)*z.';
rho([n+1:n+1:end n]) = -alpha./(beta.*nu);
rho(1:n+1:end) = 0;
rho = rho+(z./beta)*alpha.';

%Adjust for decision
insert1 = rho(1,2);
insert2 = rho(1,3);
rho(dec_pre,dec) = insert1;
for i = dec
    dec_loc = dec(dec~=i);
    rho(i,dec_loc) = insert2;
end
rho(dec,loc(end)+1) = insert1;

end
