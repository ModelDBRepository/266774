%% structure parameters
num_columns = 4;
npp = 100;
ppc = 2*npp;
N = ppc*num_columns;
pop = num_columns*2;
n = N;
unit_num = 1:2:(2*num_columns);
num_trials = 100;

%% time parameters
dt = 1; %time step
T = 50/dt; %time stimulated
delta = 600/dt; %time in between stimulation
t_total = 5001; %time of trial
D = 10/dt; %intrinsic delay
tau_w = 40; %estimate time window

%% input parameters
eG = .1; %input gain excitatory
iG = .1; %input gain inhibitory
W_in = input_weights(num_columns,npp,N,n,eG);
W_inin = (iG/eG)*W_in;
unit_stim = unit_num;
t_stim = [1 501 1501 2201 4001];
p_r = .03*dt; %poisson rate 40 Hz

%% weight matrix parameters
l5_rec = .00012;%.00012;%.00015;%.645; DO NOT SET TO ZERO OR RECURRENT IDENTITY WILL BE WRONG
l23_rec = .00000;
l5_l23 = .0002; %.0002
l23_l5 = .0000002;%.4; %DO NOT SET TO ZERO OR FF IDENTITY WILL BE WRONG
i_l5_rec = .0001;
i_l23_rec = .000;
i_l5_l23 = .07;%.05
i_l23_l5 = 0;
i_l23_l23 = .1;
i_l5_l5 = .1;
p_scale = -.1;
p_l5_rec = .0002;
p_l23_rec = .001;
p_l5_l23 = 0;
p_l23_l5 = 0;
W_ji = Sparse_L_ij(num_columns,npp,N,l5_rec,l23_rec,l5_l23,l23_l5,0,0,.000);
M_ki = Sparse_L_ij(num_columns,npp,N,i_l5_rec,i_l23_rec,i_l5_l23,i_l23_l5,i_l23_l23,i_l5_l5,.0001);
P_ik = Sparse_L_ij(num_columns,npp,N,p_l5_rec,p_l23_rec,p_l5_l23,p_l23_l5,0,0,.0001);
rec_identity = W_ji.*L_ij_no_rand(num_columns,npp,N,1,0,0,0,0)>0;
ff_identity = W_ji.*L_ij_no_rand(num_columns,npp,N,0,0,0,1,0)>0;
total_identity = rec_identity + ff_identity;

%% membrane dynamics parameters
rho = 1/7; %percentage change of synaptic activation with input spikes
tau_se = 80; %time constant for excitatory synaptic activation
tau_si = 10; %time constant for inhibitory synaptic activation
norm_noise = 1e1;
C_m = .2; %membrane capacitance
g_L = .01; %leak conductance
E_i = -70;
E_l = -60; %leak reversal potential
E_e = -5; %excitatory reversal potential
v_th = -55; %threshold potential
v_th_i = -50; 
v_rest = -60; %resting potential
v_hold = -61; %return potential
t_refractory = 2/dt;

%% Learning parameters
tau_p = 2000; %LTP time constant
tau_d = 1000; %LTD time constant
T_max_p = .0033; %maximum eligibility trace for LTP
T_max_d = .00345; %maximum eligibility trace for LTD
eta_p = 45*3500; %
eta_d = 25*3500; %
eta_rec = .002; %learning rate
eta_ff = .25; %learning rate
tau_p1 = 200;% tau_p1 = 200; %LTP time constant
tau_d1 = 800;% tau_d1 = 800; %LTD time constant
T_max_p1 = .0034;% T_max_p1 = .003; %maximum eligibility trace for LTP
T_max_d1 = .00345;% T_max_d1 = .0034; %maximum eligibility trace for LTD
eta_p1 = 20*3500;% eta_p1 = 25*3000; %
eta_d1 = 15*3500;% eta_d1 = 15*3000; %
trace_refractory = zeros(N,N);

eta_p_rec = 1;
K = 1;

%% reward parameters
delay_time = 25/dt; %added reward delay
t_reward = t_stim(1:end)+(delay_time);
rew_vect = zeros(t_total,1);
del_time = round(delay_time/2);
for m = t_reward
    rew_vect(m-del_time:m+del_time) = 1;
    rew_vect(m+del_time+1) = 2;
end
