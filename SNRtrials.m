%April 2024
%Remaking Horchler Figure (part B)

close all; clear all; clc

x0 = [0.9 0 0 0];
dt = 0.1;
alpha = ones(1,4);
beta = 2;
nu = 1;
self = x0;

eps = 1e-3;

neighborhood = 0.8;

Nt = 80;
decision = [2 3]; %which nodes do we want the decision split to be between
rho = decision_rho(alpha,beta,nu,decision);


figs = 1; %change to 1 if you want plots

plotcolors = [1 0 0; 0 0 0; 0 0 1; 0 1 0];
trials = 10;
biasval = [1e-5 5e-4 1e-3];



for j = 1:length(biasval)
    for k = 1:trials
        x = x0;
        rng('shuffle')
        % rng(1)
        for i = 1:Nt-1
            noise(i,:) = eps*sqrt(dt)*randn(1,size(rho,1));%gaussian noise
            noise(i,2) = noise(i,2) + biasval(j);
            dx = x(i,:).*(alpha - x(i,:)*rho)*dt + noise(i,:);
            %Euler-Maruyama integration
            x(i+1,:) = max(min(x(i,:) + dx, 1), 0.0005);
        end
    
        figure(1)
        for i = 2
            hold on
            plot((1:Nt)*dt,x(:,i),'-','Color',plotcolors(i,:))
        end
    end
    pause(1)

end
yline(neighborhood,'--',{'Decision';'Threshold'},'Color',[.5 .5 .5],...
    'Linewidth',2,'labelverticalalignment','middle',...
    'labelhorizontalalignment','left','Fontsize',14)
xlabel('Time (s)','Fontsize',16)
ylabel('State Activation, x','Fontsize',16)
hold off


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




