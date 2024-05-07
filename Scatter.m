%November 2023
%Using Horchler neighbourhood to define a decision (blue or black) or
%indecision.

%Changing the passage time condition to i(1) = "when x >= neighborhood2"
%and i(end) = "when x < neighborhood2" - 1

%Checked that decision can be moved to nodes 3 & 4.

%Generalize construction of new decision-making rho matrix. (Check for 4
%nodes).

%NEXT! CHECK FOR 5 NODES. (Possible decision between 3 nodes?)

%January 2024
%Adjusted for scalable node figures


close all; clear all; clc
% rng('shuffle')

x0 = [0.9 0 0 0];
dt = 0.1;
alpha = ones(1,4);
beta = 1;
nu = 1;
eps = 1e-3;
neighborhood = 0.8;

Nt = 200;

decision = [2 3]; %which nodes do we want the decision split to be between

rho = decision_rho(alpha,beta,nu,decision);

x = x0;
figs = 1; %change to 1 if you want plots

black = 0;
blue = 0;
other = 0;

noiseincrement = [0 1e-5 5e-5 1e-4 5e-4 1e-3];

noisediff = nan(1,Nt-1);

transitiontime_blue = nan(length(noiseincrement),10);
transitiontime_black = nan(length(noiseincrement),10);
transitiontime_neither = nan(length(noiseincrement),10);
avgtransitiontime_blue = nan(1,length(noiseincrement));
avgtransitiontime_black = nan(1,length(noiseincrement));
choice = 1;

newplotcolor = [0 0 0; 0 0 1; .5 .5 .5];

for j = 1:length(noiseincrement)
    for k = 1:25
        % k = 0;
        for blue_black_switch = 2:3 %black bias then blue bias
        
            rng('shuffle')
            tp2 = []; %passage time
            tp3 = [];
            tp1 = [];
            tp4 = [];
        
            for i = 1:Nt-1
                noise(i,:) = eps*sqrt(dt)*randn(1,4);
                noise(i,blue_black_switch) = noise(i,blue_black_switch) + noiseincrement(j);
                dx = x(i,:).*(alpha - x(i,:)*rho)*dt + noise(i,:); %with dt
                %Euler-Maruyama integration
                x(i+1,:) = max(min(x(i,:) + dx, 1), 0.0005);
                noisediff(i) = abs(noise(i,2)-noise(i,3));

                %Check nodes above decision threshold
                if x(i+1,1)>=neighborhood
                    tp1 = [tp1 i];
                end
                if x(i+1,2)>=neighborhood
                    tp2 = [tp2 i];
                end
                if x(i+1,3)>=neighborhood
                    tp3 = [tp3 i];
                end
                if x(i+1,4)>=neighborhood
                    tp4 = [tp4 i];
                end
                
            end

            sumnoisediff = sum(noisediff) / Nt;
            choice = 1;
        
            black_diff = diff(x(:,2));
            black_vel = black_diff/dt;
            blue_diff = diff(x(:,3));
            blue_vel = blue_diff/dt;
            blueblack_diff = abs(black_diff)-abs(blue_diff);
            
            plotcolors = [1 0 0; 0 0 0; 0 0 1; 0 1 0];

            try
            if any(tp2)
                transitiontime_black(j,k) = (tp2(1)-tp1(end))*dt;
                transitiontime = transitiontime_black(j,k);
                scattercolor = [0 0 0];
                passagetime_black(j,k) = (tp2(end)-tp2(1))*dt;
                if noiseincrement(j)~=0 && blue_black_switch == 3
                    choice = 0;
                end
            elseif any(tp3)
                transitiontime_blue(j,k) = (tp3(1)-tp1(end))*dt;
                transitiontime = transitiontime_blue(j,k);
                scattercolor = [0 0 1];
                passagetime_blue(j,k) = (tp3(end)-tp3(1))*dt;
                if noiseincrement(j)~=0 && blue_black_switch == 2
                    scattermkr = "x";
                    choice = 0;
                end
            elseif ~any(tp2) && ~any(tp3)
                transitiontime_neither(j,k) = (tp4(1)-tp1(end))*dt;
                transitiontime = transitiontime_neither(j,k);
                scattercolor = [.5 .5 .5];
                choice = .5;
            end
            catch ME
                if (strcmp(ME.identifier,'MATLAB:badsubscript'))
                    warning('Red magnitude not large enough.')
                    figs = 1;
                    % continue
                end
            end
            k = k+1;

            % figure(1)
            % clf(1)
            % figure(1)
            % for i = 1:size(rho,1)
            %     hold on
            %     plot((1:Nt)*dt,x(:,i),'.-')
            % end
            % colororder(plotcolors)
            % hold off
            snr = (noiseincrement(j)^2)/(0.1*eps^2);


            % figure(3)
            % pause(0.1)
            % hold on
            % if choice == 1
            %     scatter(sumnoisediff,transitiontime,...
            %         'markeredgecolor',scattercolor,...
            %         'markerfacecolor',scattercolor)
            % elseif choice == 0
            %     scatter(sumnoisediff,transitiontime,'x',...
            %         'markeredgecolor',scattercolor)
            % elseif choice == 0.5
            %     scatter(sumnoisediff,transitiontime,...
            %         'markeredgecolor',scattercolor)
            % end
            % ylabel('time to decision')
            % xlabel('sum noise diff')

            figure(4)
            pause(0.1)
            hold on
            if choice == 1
                scatter(snr,transitiontime,...
                    'markeredgecolor',scattercolor,...
                    'markerfacecolor',scattercolor)
                set(gca,'xscale','log')
            elseif choice == 0
                scatter(snr,transitiontime,'x',...
                    'markeredgecolor',scattercolor)
                set(gca,'xscale','log')
            elseif choice == 0.5
                scatter(snr,transitiontime,...
                    'markeredgecolor',scattercolor)
                set(gca,'xscale','log')
            end
            ylabel('time to decision')
            xlabel('snr')
        end
    end

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

