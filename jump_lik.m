% estimate the likelihood for a single step that jumps across the boundary
% of different diffusion coefficients and elastic forces
% Sohyeon Park, Jun Allard, allardlab.com

verbose = 1; % set to 1 for debugging

% incoming parameters
dt_obs = 0.05; % time between observations

z_max = 2; % domain size

xi = 1; % viscosity, pNs/um
kBT = 4.2e-3; % pNum
UTZ = 0*1e-1; % height of triangle energy barrier
z_TZ = 0.5; % center of tz
width_TZ = 0.1; % width of tz

xi_TZ=2*xi;

% internal parameters
dz_numerical = 0.001; %numerical spatial step -- This controls accuracy.

% -- consequential parameters

D = kBT/xi; % diffusion coefficient in bulk
v = UTZ/(width_TZ*xi_TZ); % magnitude of velocity (not its sign/direction) in TZ

D_TZ = kBT/xi_TZ; % diffusion coefficient in TZ

dt_numerical_estimate = min( [dz_numerical^2/D, dz_numerical/v] );
dt_numerical = 16*dt_numerical_estimate; % this can be quite large since we're using Backward Euler

if (verbose) % tests
    Peclet = (v*dz_numerical)/D;
    DZR = v*sqrt(dt_numerical/2*D);
    
    display(DZR)
    display(Peclet)
end

ntmax = round(dt_obs/dt_numerical);

% -- build matrices for Backward Euler
nz = ceil(z_max/dz_numerical)+1;

z_array = 0:dz_numerical:z_max;
n_TZ    = floor( (z_TZ-(width_TZ/2))/dz_numerical ):ceil( (z_TZ+(width_TZ/2))/dz_numerical );
n_down  = floor( (z_TZ-(width_TZ/2))/dz_numerical ):floor( (z_TZ)/dz_numerical-1 );
n_up    = ceil( (z_TZ)/dz_numerical+2 ):ceil( (z_TZ+(width_TZ/2))/dz_numerical );


% Diffusion matrix
% NOTE Is this the Ito or Stratanovich diffusion matrix?
D_array = D*ones(nz,1);
D_array(n_TZ) = D_TZ;

D_matrix = 1/(dz_numerical^2)*( -2*diag(D_array) + diag(D_array(2:end),+1)+ diag(D_array(1:end-1),-1) );

% periodic BCs for convenience
D_matrix(1,end) = -sum(D_matrix(:,end));
D_matrix(end, 1) = -sum(D_matrix(:,1));

if verbose
    figure(61); clf;
    spy(D_matrix);
    any(sum(D_matrix))   
end

% Advection matrix
v_matrix = zeros(nz);
v_matrix(n_down,n_down) = v/dz_numerical*( -1*diag(ones(numel(n_down),1),0) + diag(ones(numel(n_down)-1,1),+1) );
v_matrix(n_up,n_up)     = v/dz_numerical*( -1*diag(ones(numel(n_up),1)  ,0) + diag(ones(numel(n_up)-1,1),  -1) );

v_matrix(min(n_down)-1,min(n_down)) = + v/dz_numerical;
v_matrix(max(n_up)+1,max(n_up)) = + v/dz_numerical;

if verbose
    figure(62); clf;
    spy(v_matrix);
    any(sum(v_matrix))   
end

% invert!

Back_matrix = inv(eye(nz) - (D_matrix + v_matrix)*dt_numerical);

% SOLUTION METHOD
% we wish to solve dpdt=Ap for IC delta functions (p0=1/dz_numerical at one bin and zero elsewhere).
% Backward Euler method: p(t+1) = A p(t+1) dt + p(t) so p(t+1)=inv(I-A*dt)p(t), and p(t)= (inv(I-A*dt))^nt p0

% direct power

tic;
pz_final = mpower(Back_matrix,ntmax);
toc;

% Solver loop - I tried this and it was slower....
if 0
    tic;
    pz_final = Back_matrix;
    for nt = 1:ntmax
        pz_final = Back_matrix*pz_final;
    end % finished time loop
    toc;
end

pz_final= pz_final/dz_numerical; % impose delta function initial conditions.
% the columns represent different ICs

if verbose
    % check for unity
    sum(pz_final*dz_numerical);
    normalization_error = sum(pz_final*dz_numerical) - 1;
    display(max(abs(normalization_error)));
end

%% analyze

figure(1); clf; hold on; box on;

plot(z_array,pz_final(:,1:100:end))
%set(gca,'ylim', [0,0.05])

sigma_t = sqrt(2*D*dt_obs);
gaussian = 1/sqrt(2*pi*sigma_t^2)*exp(-(z_array-v*dt_obs-1).^2/(2*sigma_t^2));

plot(z_array,gaussian, '--k')
