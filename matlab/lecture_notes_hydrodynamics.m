clearvars
close all

%% Definition of backwater objects
R(3)=Backwater; % River has three pieces

[R(:).So]=deal(1e-4); % set the slope
[R(:).Chez]=deal(50); % set the roughness
[R(:).b]=deal(200); % set the width
[R(:).Q]=deal(800); % set the discharge

%% Evaluate initial situation, before river widening

% The river is in equilibrium so the boundary condition is at equilibrium
% depth
[R(:).a0]=deal(R(:).a_equilibrium); % set the boundary depth

% The reaches are contiguous, so we can match the start and end of the
% reaches. The x0 is where the boundary condition is imposed to the reach,
% so for this subcritical case this always has to be downstream. So for the
% downstream reach we have:
R(3).x0=0;

%Let's assume that the downstream reach has a length of 5 km, the we have:
R(3).x_end=-5e3;

% The end of the downstream piece is the start of the middle piece:
R(2).x0=R(3).x_end;
R(2).zb0=R(3).bed_level(end);

% This piece has a length of 20 km so:
R(2).x_end=R(2).x0-20e3;

% Now we define the most upstream piece:
R(1).x0=R(2).x_end;
R(1).zb0=R(2).bed_level(end);

% Let's make this piece also 20 km:
R(1).x_end=R(1).x0-20e3;

%% Now we will change the width in reach 2
R(2).b=220;
figure
[x,a]=R(2).solve;
R(1).a0=a(end);
plot(R)




