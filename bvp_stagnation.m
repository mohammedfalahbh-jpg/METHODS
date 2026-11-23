clc
clear all
close all
% clc
 
eta_inf = 7;       % y value for the infinity (far from wall) loction
n_eta = 100;        % number of grid points in the wall normal direction

% Y0 --> f
% Y1 --> f'
% Y2 --> f''

% Set initial condition of Y0, Y1, Y2 to be all 0, all 1 and all 0,
% respectively:
solinit = bvpinit(linspace(0,eta_inf,n_eta),[0 0 0]);


% Solve system of ODE
sol = bvp4c(@sode, @bc,solinit);     
eta = sol.x;
Y = sol.y;      % assign solution (Y0,Y1,Y2) to Y matrix (3,n_eta)

BB=eta.';


fprintf('\n');
fprintf('Y2(0) computed using infinity = %g is %7.5f. \n',eta_inf,Y(3,1))

% Plot results
% Set figure parameters for plot:
iplots   = 8;
monitor  = 0.6;
scrsz    = get(0,'ScreenSize');
scrsz(3) = 1260;
width    = monitor*scrsz(3);
height	 = width;
fig_loc = [100 , 100];
fig_size = [0.8*width, 0.6*width];

% Genrate plot of Y0=f, Y1=f', Y2=f'' along eta:
h=figure(1);
set(h,'Position',[fig_loc, [1 1].*fig_size],'PaperPositionMode','auto');
set(h,'DefaultAxesFontsize',20);set(h, 'DefaultLineLineWidth', 2);

eta = linspace(0,eta_inf,n_eta);    % new eta range for plot
Y0 = deval(sol,eta,1);       % evaluate Y0
Y1 = deval(sol,eta,2);       % evaluate Y1
Y2 = deval(sol,eta,3);       % evaluate Y2
plot(Y0, eta,'k', Y1,eta,'r',Y2,eta,'b');
% --------% Differential equations,Y0'=Y1,  Y1'=Y2 ,  Y2'=-(1/2)*Y1*Y2 .  
function dYdeta = sode(eta,Y)
% Note: Y(1), Y(2), Y(3) correspond to Y0, Y1, Y2
dYdeta = [
    BB*Y(2)                    %RHS of Y0'=Y1
    Y(3)                    %RHS of Y1'=Y2
    -(1/2)*Y(1)*Y(3)   %RHS of Y2'=-(1/2)*Y0*Y2
    ];
end

% -------%Boundary conditions Y0(0)=0, Y1(0)=0, Y1(infty) = 1 --> Y2(infinity)-1=0
function res = bc(Y0,Yinf)
%Compute the residual at boundary
% Note: Y0(1), Y0(2), Y0(3) correspond to Y0, Y1, Y2 boundary values @ eta = 0.
% Note: Yinf(1), Yinf(2), Yinf(3) correspond to Y0, Y1, Y2 boundary values @ eta = inf.
res = [
    Y0(1) - 0           % to set Y0(0)=0
    Y0(2) - 0           % to set Y1(0)=0
    Yinf(2) - 1         % to set Y2(0) so that Y1(infty)=1
    ];
end  
axis([0 2 0 eta_inf]);
legend({'$Y_1=f(\eta)\sim -v$', '$Y_2=f''(\eta)\sim u$','$Y_3=f''''(\eta)\sim \tau$'},'interpreter','latex',...
    'location','Northwest')
legend('boxoff')
grid on
title('Plane stagnation point flow')
xlabel('Solution')
ylabel('Distance from Wall (\eta)')
hold on
        
% find eta 0.99
index1 = find(Y1 <= 0.99, 1, 'last');
index2 = find(Y1 >= 0.99, 1, 'first');
A1=Y1(index1);
A2=Y1(index2);
B1=eta(index1);
B2=eta(index2);
eta99=B1 + (B2 - B1) * ((0.99 - A1) / (A2 - A1));
fprintf('eta_0.99 = %g',eta99)
