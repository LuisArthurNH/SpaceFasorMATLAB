
%mag_fied.m
clc; clear variables; close all;

bmax = 1;
freq = 50;
w = 2*pi*freq;
t=0:1/6000:1/freq;

n = 1; % number of rotations that the mag field will do

% Allow plot(true/1) or not (false/0):
plot_abc = 1;
plot_alpha_beta = 1;
plot_dq = 1;

%% Important Angles: Rho and ThetaZero

% RHO: Angle between d and alpha
% rho = pi/4;
% rho = 0.2*pi/2;
rho = pi/6;
% rho = -pi/2;
% rho = 0;

% THETAZERO: Initial Space Vector Phase
% thetazero = pi/3.35;
thetazero = pi/6; 
% thetazero = 0; 
% thetazero = -pi/2; 
% thetazero = rho;

%% Definition of time variant functions
% Electrical quantities

% Definition of Magnetic Fields: 
% addition of pi/2 for the space fasor to start at origin
% thetazero is the initial angle between alpha and the space vector

Baa = sin(w*t + pi/2 + thetazero)         .*(cos(0)+1i*sin(0));
Bbb = sin(w*t - 2*pi/3 + pi/2 + thetazero).*(cos(2*pi/3)+1i*sin(2*pi/3));
Bcc = sin(w*t + 2*pi/3 + pi/2 + thetazero).*(cos(-2*pi/3)+1i*sin(-2*pi/3));

% Definition of SPACE FASOR - Invariant Amplitude Factor
% Bnet = 1*(Baa + Bbb + Bcc);
Bnet = (2/3)*(Baa + Bbb + Bcc);




%% Definition of alpha and beta axis and arrows
% Definition of Alpha and Beta Axis
Aa = [-1.7 1.7];
Ab = [-1.7 1.7];

% Definition of Alpha and Beta Arrows
% alpha
Aa1x = [1.6 1.69];
Aa1y = [-0.05 0];
Aa2x = [1.6 1.69];
Aa2y = [0.05 0];
%beta
Ab1y = [1.6 1.69];
Ab1x = [-0.05 0];
Ab2y = [1.6 1.69];
Ab2x = [0.05 0];

% Alpha and Beta Components

alpha = real(Bnet); 
beta  = imag(Bnet);

%% Definition of the Rotating Axis - DQ

% Definition of AXIS d
daxis = 1.7.*(cos(w*t + rho) + 1i*sin(w*t + rho));
% Definition of AXIS q 
qaxis = 1.7.*(cos(w*t + rho  + pi/2) + 1i*sin(w*t + rho + pi/2));

% angle between d and alpha
theta = w*t + rho;

% Definition of COMPONENT d
d = (alpha.*cos(theta) + beta.*sin(theta)).*(cos(w*t + rho) + 1i*sin(w*t + rho));
% Definition of COMPONENT q
q = (-sin(theta).*alpha + cos(theta).*beta).*(cos(w*t + rho + pi/2) + 1i*sin(w*t + rho+pi/2));


%% Reference Circles
circle =  1*(cos(w*t) + 1i*sin(w*t));
% circle1 = 1.5*circle;



%% Plotting

  
if (plot_abc)
    figure(1);   % 3ph Magnetic Fields and Resultant one
    for ii = 1:n 

        for x = 1:length(t)
    %         plot (circle1,'k');      % se essa linha está comentada,comentar hold on
    %         hold on 
            plot (circle,'k', 'LineWidth', 1.0);
            title('3ph syst and Space Fasor')
            hold on

            % Plot axis ABC
            plot([0 1.7*cos(2*pi/3)],  [0 1.7*sin(2*pi/3)],'--k', 'LineWidth',1.0);
            hold on
            plot([0 1.7*cos(-2*pi/3)], [0 1.7*sin(-2*pi/3)],'--k', 'LineWidth',1.0);
            hold on
            plot([0 1.7], [0,0],'--k', 'LineWidth',1.0);
            hold on





            % Plot Magnetic Fields - Ba, Bb, Bc
            plot([0 real(Baa(x))], [0 imag(Baa(x))],'r', 'LineWidth',1.5);
            hold on
            plot([0 real(Bbb(x))], [0 imag(Bbb(x))],'b', 'LineWidth',1.5);
            hold on
            plot([0 real(Bcc(x))], [0 imag(Bcc(x))],'g', 'LineWidth',1.5);
            hold on

            % Plot Resultant Magnetic Field - Bnet
            plot([0 real(Bnet(x))],[0 imag(Bnet(x))],'m', 'LineWidth',1.5);


            axis square;
            axis([-2 2 -2 2]);
            drawnow;
            hold off

        end

    end
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plot_alpha_beta)
    figure(1);   % Rotating Magnetic Field and Alpha Beta  
    for ii = 1:n 

        for x = 1:length(t)
    %         plot (circle1,'k');      % se essa linha está comentada,comentar hold on
    %         hold on 
            plot (circle,'k', 'LineWidth', 1.0);
            title('Alpha Beta Frames')
            hold on

            % Plot axis alpha and beta
            plot(Aa, [0 0],'--k', 'LineWidth',1.0);
            hold on
            plot([0 0], Ab,'--k', 'LineWidth',1.0);
            hold on
            % plot arrows alpha axis
            plot(Aa1x, Aa1y ,'k', 'LineWidth',1.0);
            hold on
            plot(Aa2x, Aa2y ,'k', 'LineWidth',1.0);
            hold on
            % plot arrows beta axis
            plot(Ab1x, Ab1y ,'k', 'LineWidth',1.0);
            hold on
            plot(Ab2x, Ab2y ,'k', 'LineWidth',1.0);
            hold on


            % Plot Resultant Magnetic Field - Bnet
            plot([0 real(Bnet(x))],[0 imag(Bnet(x))],'r', 'LineWidth',1.5);

            % Plot Alpha and Beta components of Bnet
            plot([0 alpha(x)],[0 0]      ,'b', 'LineWidth',1.5);
            plot([0 0]       ,[0 beta(x)],'g', 'LineWidth',1.5);


            axis square;
            axis([-2 2 -2 2]);
            drawnow;
            hold off

        end

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

if (plot_dq)
    figure(1);   % Rotating Magnetic Field and QD axis
    for ii = 1:n 
        for x = 1:length(t)

            % Plot reference Circle
            plot (circle,'k', 'LineWidth', 1.0);
            title('DQ Reference Frame')
            hold on  

            % Plot axis alpha and beta
            plot(Aa, [0 0],'--k', 'LineWidth',1.0);
            hold on
            plot([0 0], Ab,'--k', 'LineWidth',1.0);
            hold on
            % plot arrows alpha axis
            plot(Aa1x, Aa1y ,'k', 'LineWidth',1.0);
            hold on
            plot(Aa2x, Aa2y ,'k', 'LineWidth',1.0);
            hold on
            % plot arrows beta axis
            plot(Ab1x, Ab1y ,'k', 'LineWidth',1.0);
            hold on
            plot(Ab2x, Ab2y ,'k', 'LineWidth',1.0);
            hold on

             % Plot d and q axis
            plot([0 real(qaxis(x))],[0 imag(qaxis(x))],'--k', 'LineWidth',1.5);
            plot([0 real(daxis(x))],[0 imag(daxis(x))],'--k', 'LineWidth',1.5);

            % Plot Resultant Magnetic Field - Bnet
            plot([0 real(Bnet(x))],[0 imag(Bnet(x))],'r', 'LineWidth',1.5);


            % Plot Alpha and Beta components of Bnet
            plot([0 real(d(x))],[0 imag(d(x))]      ,'b', 'LineWidth',1.5);
            plot([0 real(q(x))],[0 imag(q(x))],'g', 'LineWidth',1.5);




            axis square;
            axis([-2 2 -2 2]);
            drawnow;
            hold off
        end

    end

end