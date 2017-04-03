%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

    %% loading data
    [A B ] =xlsread('stockreturns.xls') ;                                  % A = daily prices; B = Stock Name
    %% Reading in the variables and defining basic quantities
    names   = transpose(B);
    logret  = transpose(A);
    [T N]   = size(logret);                                                % T = No. Stocks; N = No. return days
    % mean and covariance
    Sigma   = cov(logret)*252;                                             % I use 252 for annualizing the cov matrix
    Vars    = diag(Sigma);                                                 % variances of the stocks
    mu      = mean(logret)'*252;                                           % mean log returns
    

% defining auxiliary variables
e       = ones(size(mu));                                                  % NOT SURE WHAT THESE VARIABLES ARE
a       = mu'/Sigma*mu;
b       = mu'/Sigma*e;
c       = e'/Sigma*e;
d       = a*c - b^2;
k1      = (c*mu - b)./d;
k2      = (a - b*mu)./d;

%% Formulae for the efficient frontier in the (muP,sigmaP)-plane
% No short selling constraint - closed form
nport   = 30;                                                             % number of portfolios to draw the efficient frontier
sigmaP  = (1/c:(max(Vars)-1/c)/nport:max(Vars))'.^(0.5);                   % vector of portfolio variance, 1/c is the MinVar!
muP     = b/c + sqrt(d./c.*(sigmaP.^2 - 1./c));
muP     = real(muP);                                                       % sometimes the MinVaR return causes Matlab to return imaginary number...
Wu      = zeros(N,nport);                                                  % unconstrained weight matrix Wu
for i = 1:length(muP)
    Wu(:,i) = Sigma\(muP(i)*k1 + k2);
end
% Short selling constraint - numerics needed
Aeq     = [mu'; ones(1,N)];
LB      = zeros(1,N);                                                      % one could define other linear constraints
UB      = ones(1,N);                                                       % without much effort...
opts    = optimset('Algorithm', 'interior-point-convex', 'Display','off');
b       = min(mu):(max(mu) - min(mu))/nport:max(mu);
Wc      = zeros(N,nport);                                                  % constrained weight matrix Wc
muC     = zeros(1,N);
sigC    = zeros(1,N);
for i = 1:length(b)
    beq             = [b(i); 1];
    [Wc(:,i) varP]  = quadprog(Sigma,[],[],[],Aeq,beq,LB,UB,UB/N,opts);
    muC(i)          = Wc(:,i)'*mu;
    sigC(i)         = sqrt(2*varP);
end

% selecting the efficient frontier for the constrained case;
muMin   = muC(sigC==min(sigC));
muEff   = muC(muC>=muMin);
sigEff  = sigC(muC>=muMin);


%% Figure plot
figure(1)
%%plot(sigmaP,muP,'r','Linewidth', 1.5)
hold on
plot(sigEff,muEff,'b','Linewidth', 1.5)
plot(sqrt(Vars),mu,'ko')
set(gca, 'Box', 'on', 'Linewidth', 1.5, 'Fontsize', 14)
xlabel('Portfolio standard deviation')
ylabel('Portfolio return')
title('Efficient Frontiers')
legend( 'no constraint', 'no shorts', 'assets' , 'Location', 'NorthWest')
grid on


%% Excel Data
xlswrite('Test1',sigEff, 'section1');

