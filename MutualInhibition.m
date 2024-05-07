%April 2024
%Increasing "other inhibition" in rho matrix construction to see if we get
%0 in the heat map!

%May 2024
%Looking for indecisions, correct decisions, and transition times

close all; clear all; clc
rng('shuffle')


x0 = [0.1 0 0 0];
dt = 0.1;
alpha = ones(1,4);
beta = 1;
nu = 1;

eps = 1e-3;
biasval = [0 1e-5 1e-4 1e-3 1e-2];

eps_plot = [1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2];
eps_plot = 1e-2;

neighborhood = 0.8;

Nt = 200;
decision = [2 3]; %which nodes do we want the decision split to be between
trials = 1;

rhoa = 1;
rhob = 0;
rhoc = 1:5;
rhoc = 5;
rhod = 2;

figs = 1; %change to 1 if you want plots

plotcolors = [1 0 0; 0 0 0; 0 0 1; 0 1 0];


for a = 1:numel(rhoc)

    rho = [rhoa rhob rhob rhod; rhod rhoa rhoc(a) rhob;...
        rhod rhoc(a) rhoa rhob; rhob rhod rhod rhoa];

    for b = 1:numel(eps_plot)
        neither = 0;
        k = 0;

        while k < trials
            t1 = [];
            t2 = [];
            t3 = [];
            t4 = [];
            x = x0;
            rng('shuffle')

            for i = 1:Nt-1
                noise(i,:) = eps_plot(b)*sqrt(dt)*randn(1,size(rho,1)); %gaussian noise
                dx = x(i,:).*(alpha - x(i,:)*rho)*dt + noise(i,:); %with dt
                %Euler-Maruyama udpate
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

            if isempty(t1) || t1(1)>100  %red did not activate enough!
                continue
            end
            if isempty(t2) && isempty(t3) %no decision made
                if isempty(t4)
                    continue
                end
            end
            k = k + 1; %Count the number of usable trials
    
            if figs == 1
                figure(1)
                clf(1)
                figure(1)
                for i = 1:size(rho,1)
                    hold on
                    plot((1:Nt)*dt,x(:,i),'.-')
                end
                colororder(plotcolors)
                hold off
                % title({
                %     ['State Space Activations'], ...
                %     ['Nt = ' num2str(Nt)], ...
                %     ['c = ' num2str(rhoc(a))...
                %     ', eps = ' num2str(eps(noisemag))] ...
                %     })
            end

            % mrkr1 = find(diff(t1)~=1);
            % t1 = t1(1:mrkr1(1));

            if any(diff(t1)~=1)
                mrkr = find(diff(t1)~=1);
                t1 = t1(1:mrkr(1));
            end

            check1 = any(t2>=t1(1)); %black
            check2 = any(t3>=t1(1)); %blue
            if ~check1 && ~check2
                neither = neither + 1;
                transtime(k) = nan;
            elseif isempty(t2) && isempty(t3)
                neither = neither + 1;
                transtime(k) = nan;
            elseif any(check1)
                transtime(k) = t2(1)-t1(end);
            elseif any(check2)
                transtime(k) = t3(1)-t1(end);
            end

        end
        neithercount(b,a) = neither;
        neithercount
        timecount(b,a) = mean(transtime,'omitmissing');
        timecount
    end
    pause(0.1)
end

% eps_plot = eps+biasval;
timecount = timecount*dt;

figure(2)
h = heatmap(rhoc,flip(eps_plot),flip(neithercount),...
    'Fontsize',14,'ColorLimits',[0 100]);
title('Rate of Indecision')
colormap(flipud(autumn))

figure(3)
heatmap(rhoc,flip(eps_plot),flip(timecount),'Fontsize',14)
title('Time to Decision (s)')







%% SHC_LV_CREATECYCLE
function rho=shc_lv_createcycle(alpha,bet,nu,direction)
%SHC_LV_CREATECYCLE  Create connection matrix RHO for Lotka-Volterra SHC cycle.
%   RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU) returns an N-dimensional
%   non-negative connection matrix RHO for the growth rates ALPHA, state
%   magnitudes BETA, and saddle values NU. ALPHA, BETA, and NU, must be symbolic
%   or floating-point scalars or length N vectors. If all three parameters are
%   scalar, a 3-by-3 connection matrix will be created. To form an SHC cycle,
%   the elements of RHO are:
%   
%                  { ALPHA(i)/BETA(i),                   if i == j
%       RHO(i,j) = { (ALPHA(i)-ALPHA(j)/NU(j))/BETA(j),  if i == mod(j,N)+1
%                  { (ALPHA(i)+ALPHA(j))/BETA(j),        otherwise
%   
%   where i,j are in {1,..., N}, N >= 3.
%   
%   RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU,-1) is the same as above except the
%   resultant connection matrix RHO is transposed, reversing the direction of
%   the cycle. RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU,1) is equivalent to
%   RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU).
%   
%   See also:
%       SHC_LV_PARAMS, SHC_LV_EIGS, SHC_LV_JACOBIAN, SHC_LV_PASSAGETIME,
%       SHC_LV_MINTRANSITIONTIME

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.

%   Andrew D. Horchler, horchler @ gmail . com, Created 2-27-14
%   Revision: 1.2, 6-2-16


% Check ALPHA
if ~(isfloat(alpha) || isa(alpha,'sym'))
    error('SHCTools:shc_lv_createcycle:InvalidAlpha',...
          'The input Alpha must be a symbolic or floating-point vector.');
end
if ~isreal(alpha) || ~all(isfinitesym(alpha(:)))
    error('SHCTools:shc_lv_createcycle:AlphaNonFiniteReal',...
         ['The input Alpha must be a finite real symbolic or floating-point '...
          'vector.']);
end
if isempty(alpha) || ~isvector(alpha)
    error('SHCTools:shc_lv_createcycle:AlphaDimensionMismatch',...
          'The input Alpha must be a non-empty vector.');
end

% Check BETA
if ~(isfloat(bet) || isa(bet,'sym'))
    error('SHCTools:shc_lv_createcycle:InvalidBeta',...
          'The input Beta must be a symbolic or floating-point vector.');
end
if ~isreal(bet) || ~all(isfinitesym(bet(:)))
    error('SHCTools:shc_lv_createcycle:BetaNonFiniteReal',...
         ['The input Beta must be a finite real symbolic or floating-point '...
          'vector.']);
end
if isempty(bet) || ~isvector(bet)
    error('SHCTools:shc_lv_createcycle:BetaDimensionMismatch',...
          'The input Beta must be a non-empty vector.');
end

% Check NU
if ~(isfloat(nu) || isa(nu,'sym'))
    error('SHCTools:shc_lv_createcycle:InvalidNu',...
          'The input Nu must be a symbolic or floating-point vector.');
end
if ~isreal(nu) || ~all(isfinitesym(nu(:)))
    error('SHCTools:shc_lv_createcycle:NuNonFiniteReal',...
         ['The input Nu must be a finite real symbolic or floating-point '...
          'vector.']);
end
if isempty(nu) || ~isvector(nu)
    error('SHCTools:shc_lv_createcycle:NuDimensionMismatch',...
          'The input Nu must be a non-empty vector.');
end

% Check dimension of parameters
len = [length(alpha) length(bet) length(nu)];
n = max(len);
if n == 1
    n = 3;
elseif any(len ~= 1 & len ~= n)
    error('SHCTools:shc_lv_createcycle:ParameterDimensionMismatch',...
         ['The input parameters, Alpha, Beta, and Nu, must scalars or '...
          'equal length vectors.']);
end

%Check DIRECTION
if nargin > 3
    if ~isscalar(direction) || ~isnumeric(direction) || ~any(direction == [1 -1])
        error('SHCTools:shc_lv_createcycle:InvalidFlag',...
              'The option Direction argument must a scalar 1 or -1.');
    end
else
    direction = 1;
end

% Expand scalars
z = ones(n,1);
alpha = alpha(:).*z;
bet = bet(:).*z;
nu = nu(:).*z;

% c >= 1, Adjust size of non-dominant stable eigenvalues, c == 1 for all equal
c = 1;

% Calculate Rho based on direction
if direction == 1
    rho = c*z*(alpha./bet).';
    rho([2:n+1:end n*(n-1)+1]) = -alpha./(bet.*nu);
    rho(1:n+1:end) = 0;
    rho = rho+alpha*(z./bet).';
else
    rho = c*(alpha./bet)*z.';
    rho([n+1:n+1:end n]) = -alpha./(bet.*nu);
    rho(1:n+1:end) = 0;
    rho = rho+(z./bet)*alpha.';
end

% Make non-negative
if ~isa(rho,'sym')
    rho(:) = max(rho(:),0);
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




