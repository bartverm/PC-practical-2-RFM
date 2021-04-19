close all
restoredefaultpath
clearvars

 addpath ~/src/backwater/

%% Create the advective case (example in the practical guide)
B=Backwater; % create a backwater object
[x,~]=B.solve; % get the x-coordinate
delta_z=zeros(size(x)); % create the initial condition to a change in the bed level of zero
delta_z(x<-1.5e4 & x > -2e4)=-2; % modify the initial condition with an erosion pit between -15 and 20 km
[x, z_b_sym, dt]=morf_solver(B,delta_z,3*30*24*3600); % Run the morphological simulation for three months
imagesc(z_b_sym) % create a color plot of the morphodynamic prediction
xlabel('time / 5.6 days')
ylabel('x / 909 m')
hc=colorbar;
ylabel(hc,'bed elevation / 1 m')
axis xy
plot_zb_sym(x,z_b_sym, dt, 5) % animation of the simulation


%% Create the diffusive case (by increasing the spatial scale of the bed disturbance)
B.x_end=-4e5; % increase length to 400 km
[x,~]=B.solve; % get the x coordinates of the backwater solution
delta_z=zeros(size(x)); % initialize the initial condition for the bed
delta_z(x<-1.5e5 & x > -3.5e5)=-2; % create a large erosion pit between -150 and -350 km
[x, z_b_sym, dt]=morf_solver(B,delta_z,20*365*24*3600); % run simulation for 20 years
figure
plot_zb_sym(x,z_b_sym, dt, 5, B.So) % animate the result removing the bed slope (otherwise it's hard to see)
