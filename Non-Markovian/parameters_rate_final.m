%Number of units
num_columns = 6 ;
N = num_columns*2; %number of excitatory units, modular net
n = N; %number of inhibitory units, modular net
npp = 1; %neurons per population, modular net

%Time parameters
dt = 1; %time step
T = 100/dt; %time stimulated
delta = 400;
t_total = 5001; %time of trial
time_t = 1:dt:t_total; %time step array

%Time constants
tau_u = 40; %excitatory time constant
tau_v = 10; %inhibitory time constant
D = 10/dt; %intrinsic delay
delay_time = 50/dt; %added reward delay
T_d = 0/dt; %Trace delay

%Trial parameters
num_trials = 75;
l_array = 1:num_trials; %array of trials to be plotted

%Eligibility trace initalization
T_ijpt = zeros(N,N,numel(time_t)); %LTP eligibility trace at each time
T_ijdt = zeros(N,N,numel(time_t)); %LTD eligibility trace at each time
LTP0 = zeros(numel(time_t),numel(l_array)); %vector for plotting average LTP in population
LTD0 = zeros(numel(time_t),numel(l_array)); %vector for plotting average LTD in population

%Activations and firing rate initialization
u_it = zeros(N,numel(time_t)); %synaptic activation of each excitatory unit i
rate_uit = zeros(N,numel(time_t)); %firing rate of each excitatory unit i
v_kt = zeros(n,numel(time_t)); %synaptic activation of each inhibitory unit k
rate_vkt = zeros(N,numel(time_t)); %firing rate of each inhibitory unit k
plot_rate_uit = zeros(N,numel(time_t),numel(l_array));
plot_rate_vkt = zeros(N,numel(time_t),numel(l_array));

%Activation function parameters
u_c = 2; %critical synpatic activation
v = 2;
theta = 0;

%Learning parameters
tau_p1 = 2000; %LTP time constant
tau_d1 = 850; %LTD time constant
T_max_p1 = .322; %maximum eligibility trace for LTP
T_max_d1 = .334; %maximum eligibility trace for LTD
eta_p1 = 8; %
eta_d1 = 3.6; %
tau_p2 = 300; %LTP time constant
tau_d2 = 800; %LTD time constant
T_max_p2 = .322; %maximum eligibility trace for LTP
T_max_d2 = .33; %maximum eligibility trace for LTD
eta_p2 = 4; %
eta_d2 = 2.2; %
eta_rec = .2; %learning rate
eta_ff = 6; %learning rate
% eta_ff = 0;
t_refractory = zeros(N,N,numel(time_t));
T_intp = zeros(N,N,N/2 + 1);
T_intd = zeros(N,N,N/2 + 1);


%Scaling parameters
scale_v = 1;
I_strength = 3/dt; %strength of stimulus


%% lsm paramters

M = 336;%280;%270;
tau_net = 100; %time constant for excitatory synaptic activation
g = 1.5; %LSM gain factor
W_ji2 = normrnd(0,g/sqrt(M),M,M); %LSM recurrent weights
r_0 = 1;
theta_c2 = 1;
theta_c = 0.0095;
O_ij = normrnd(0,g/M,4*M,M); %LSM to sparse weights
O_ij = O_ij.*round(.525*rand(4*M,M)); %sparse selection of weights
% O_ij = O_ij.*round(.502*rand(4*M,M));
O_tot = O_ij~=0;
O_tot2 = sum(O_tot);
Q_ij = zeros(4*M,N); %sparse net -> modular timer weights
del_Q_ij = 0;
tau_q = 100000;