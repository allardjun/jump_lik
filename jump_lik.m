% estimate the likelihood for a single step that jumps across the boundary
% of different diffusion coefficients and elastic forces
% Sohyeon Park, Jun Allard, allardlab.com


dt_obs = 0.05; % time between observations

dz_numerical = 0.001; %numerical spatial step -- This controls accuracy.

z_max = 2; % domain size


xi = 1; % viscosity, pNs/um
kBT = 4.2e-3; % pNum
UTZ = 1e-1; % height of triangle energy barrier
z_TZ = 1; % center of tz
width_TZ = 0.1; % width of tz

% -- consequential parameters

D = kBT/xi;
v = UTZ/(width_TZ*xi);

Peclet = (v*dz_numerical)/D;

dt_numerical_estimate = min( [dz_numerical^2/D, dz_numerical/v] );
dt_numerical = dt_numerical_estimate;

DZR = v*sqrt(dt_numerical/2*D);

display(DZR)
display(Peclet)

ntmax = round(dt_obs/dt_numerical);

% -- build matrices for Backward Euler
nz = ceil(z_max/dz_numerical)+1;

% Diffusion matrix

D_matrix = D/(dz_numerical^2)*( -2*diag(ones(nz,1)) + diag(ones(nz-1,1),-1)+ diag(ones(nz-1,1),+1) );

% Advection matrix

v_matrix = v/dz_numerical*(  -1*diag(ones(nz,1)) + diag(ones(nz-1,1),+1) );

% invert!

Back_matrix = inv(eye(nz) - (D_matrix + v_matrix)*dt_numerical);

% direct power

tic;
pz_final = mpower(Back_matrix,ntmax);
toc;

pz_final= pz_final/dz_numerical;

% check for unity
%sum(pz_final*dz_numerical)

% Solver loop
if 0
    tic;
    pz_final = Back_matrix;
    for nt = 1:ntmax
        pz_final = Back_matrix*pz_final;
    end % finished time loop
    toc;
end

%% analyze

z_array = 0:dz_numerical:z_max;


figure(1); clf; hold on; box on;

plot(z_array,pz_final(1:100:end,:)')
%set(gca,'ylim', [0,0.05])

sigma_t = sqrt(2*D*dt_obs);
gaussian = 1/sqrt(2*pi*sigma_t^2)*exp(-(z_array-v*dt_obs-1).^2/(2*sigma_t^2));

plot(z_array,gaussian, '--k')
