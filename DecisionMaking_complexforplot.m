% October 2023
%Rewriting the SMP decision-making file

% January 2024
%Small changes for figures.
%x0: 0.1->0.9
%run: 30->1
%tiled layout: 2 subfigs->just first subfigure

%New updates for scalable decision-making

close all; clear all; clc
rng('shuffle')

x0 = [0.9 0 0 0 0];
dt = 0.15;
alpha = ones(1,5);
beta = 1;
nu = 1;

Nt = 200;
% rho = [.5 0 0 2; 2 .5 2 0; 2 2 .5 0; 0 2 2 .5];
decision = [4,5];
% decision = [3,4];
rho = decision_rho(alpha,beta,nu,decision)
x = x0;

black = 0;
blue = 0;
neither = 0;
other = 0;

for run = 1:1

for i = 1:Nt-1
    noise = 0.001*randn(1,length(alpha));
    % noise(3) = noise(3) + 0.0001;
    dx = x(i,:).*(alpha - x(i,:)*rho)*dt + noise; %with dt
    for j = 1:4
        comp(i,:,j) = x(j)*rho(j,:);
    end
    %dx needs to be bound? something's up with Matlab's randn
    x(i+1,:) = max(min(x(i,:) + dx, 1), 0.0005);
    % x(i+1,:) = x(i,:) + dx;
end

black_diff = diff(x(:,2));
blue_diff = diff(x(:,3));

plotcolors = [1 0 0; 0 0 0; 0 0 1; 0 1 0; 1 0 1];

figure(1)
t = tiledlayout(length(alpha),2,'TileSpacing','tight');
nexttile(1,[length(alpha),1])
for i = 1:length(alpha)
    hold on
    plot((1:Nt)*dt,x(:,i),'.')
end
hold off
colororder(plotcolors)
title('State Space Activations')

% nexttile
% hold on
% plot(black_diff,'Color',[0 0 0])
% plot(blue_diff,'Color',[0 0 1])
% title('Black & Blue Derivatives')


if all(black_diff<0.011) && all(blue_diff<0.011)
    'neither wins'
    neither = neither+1;
    neither_x = x;
    pause
    % continue    
elseif any(black_diff>=0.01)
    'black wins'
    black = black + 1;
elseif any(blue_diff>=0.01)
    'blue wins'
    blue = blue + 1;
else
    'other'
    other = other + 1;
    pause
    % continue
end
% if any(black_diff>0.01) && any(blue_diff>0.01)
%     'neither wins'
%     neither = neither + 1;
%     neither_x(:,:,neither) = x;
% end





barcolors = [.5 .5 .5; plotcolors];
timestep = [10:30]*dt;

% figure(2)
% tiledlayout('flow')
% for i = 10:30
%     compshape = reshape(comp(i,:,:),[4,4])';
%     nexttile
%     bar([alpha*ones(4,1) compshape])
% 
%     title(['t = ', num2str(timestep(i-9))])
% end
% colororder(barcolors)
% legend('alpha','kernel1','kernel2','kernel3','kernel4')




% nexttile
% for i = size(neither_x,3)
%     for j = 1:4
%         hold on
%         plot(plot((1:Nt)*dt,neither_x(:,j,i),'.'))
%     end
% end
% hold off
% colororder(plotcolors)
% title('Neither Wins')

% figure(2)
% title('State Space')
% tiledlayout(4,1)
for j = 1:size(rho,1)
    nexttile
    hold on
    plot((1:Nt)*dt, x(:,j), '-k','Color',plotcolors(j,:),'Linewidth',2)
    plot((1:Nt)*dt, x(:,j), '-.k','Color',plotcolors(j,:),'Linewidth',2)
    axis([0 30 0 1])
    title(['Node ' num2str(j)])

end
xlabel(t,'Time (s)')
ylabel(t,'State Activation, x')
set(gcf,'Position',[100 100 800 400])

end

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





