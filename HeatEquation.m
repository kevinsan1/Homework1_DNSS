%% Heat_equation
% Program to solve the diffusion equation
% using the Backward Euler method
%% Parameters
clear all;close all;clc;
savePath = ['/Users/kevin/SkyDrive/KTH Work/' ...
    'Period 3 2014/DN2255/Homework/1/Heat Equation/Figures'];
N = 50; % Number of grid points
L = 1; % The system extends from (x)=(0) to (x)=(L)
h = L/N;
i = 1:(N); % 1:N
[x,y] = meshgrid(h/2:h:L,h/2:h:L);
w = 0.2;
xs = 0.5;
ys = 0.5;
tfinal = .5;
tau = .1*h;
coeff = tau/h^2;
tsteps = ceil(tfinal/tau);
time = linspace(0,tfinal,tsteps);
%% Initialize Source function
xExponent = (x-xs).^2;
yExponent = (y-ys).^2;
S = exp(-(xExponent)/w^2).*exp(-yExponent/w^2);
deltaFunction = zeros(N);
deltaFunction(round(N/2),round(N/2))=2;
S = deltaFunction;
S = reshape(S,[N^2,1]);
%% Initialize Q-matrix
Q = zeros(N^2,1);
%% Compute matrix A
TN = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
% Boundary conditions
TN(1,1)=1;
TN(end,end)=1;
TNxN = kron(eye(N),TN) + kron(TN,eye(N));
mA = eye(N^2) + coeff*TNxN;
sparseA = sparse(mA);
%% Initialize loop and plot variables
Qplot(:,1) = Q; % initial value
stepNumber=round(.25/tau);
%% Main loops
for iter=1:stepNumber
    Q = sparseA\Q + tau*S;
    Qplot(:,iter+1) = Q(:);
end
% Loop after source is gone
for iter2=(iter+2):tsteps
    Q = sparseA\Q;
    Qplot(:,iter2) = Q(:);
end
%% Reshape Q for plotting
Qresh = reshape(Qplot,[N,N,tsteps]);
%% look at dx*dy*Qij for Conservation
cons(1:tsteps,1) = h^2*sum(sum(Qresh(:,:,1:end)));
figure(1);
plot(tau*(1:length(cons)),cons);
axis([0 tau*length(cons) min(cons) 1.2*max(cons)])
% hTitle, hXLabel, hYLabel
hTitle = title('Conservation of Q_{i,j} over time');
hXLabel = xlabel('time (sec)');
hYLabel = ylabel('\Deltax \Deltay Q_{i,j}');
% Configuration
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set( gca             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );
set(gca, ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...
    'XMinorTick'  , 'on'          , ...
    'YMinorTick'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'ZColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1             );
%% Print Plots
fin = length(Qplot(1,:));
maxZ = max(max(max(Qresh)));
%     n1 = [1 ceil(.25/(2*tau)) ceil(.25/tau)...
%        ceil(.25/tau)+5 ceil(.25/tau)+20 ceil(.25/tau)+50];
n1 = [1 floor(fin*6/24) floor(fin*12/24)...
    floor(fin*13/24) floor(fin*14/24) floor(fin*15/24)];
figure(2)
for i = 1:length(n1)
    s1=subplot(3,2,i);
    mesh(x,y,Qresh(:,:,n1(i)))
    axis([0 1 0 1 0 maxZ]);
	% hTitle, hXLabel, hYLabel, hZLabel
	hTitle = title(sprintf('t = %0.4f',time(n1(i))));
    hXLabel = xlabel('x');
    hYLabel = ylabel('y');
    hZLabel = zlabel('Q');
	% Configuration
    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set([hTitle, hXLabel, hYLabel, hZLabel], ...
        'FontName'   , 'AvantGarde');
    set( gca             , ...
        'FontSize'   , 8           );
    set([hXLabel, hYLabel, hZLabel]  , ...
        'FontSize'   , 10          );
    set( hTitle                    , ...
        'FontSize'   , 12          , ...
        'FontWeight' , 'bold'      );
    set(gca, ...
        'Box'         , 'off'         , ...
        'TickDir'     , 'out'         , ...
        'TickLength'  , [.02 .02]     , ...
        'XMinorTick'  , 'on'          , ...
        'YMinorTick'  , 'on'          , ...
        'XColor'      , [.3 .3 .3]    , ...
        'YColor'      , [.3 .3 .3]    , ...
        'LineWidth'   , 1             );
end
%% Save Figure
saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work' ...
    '/Period 3 2014/DN2255/Homework/1/Heat Equation/Figures/'];
addpath(saveFigurePath);
printYesNo = 0;
if printYesNo == 1
    set(figure(1), 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('deltaConservationPlot')]);
end
%% Save Figure 2
saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work' ...
    '/Period 3 2014/DN2255/Homework/1/Heat Equation/Figures/'];
addpath(saveFigurePath);
printYesNo = 0;
if printYesNo == 1
    set(figure(2), 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('deltaFunctionPlot')]);
end