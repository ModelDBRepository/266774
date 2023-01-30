%% initialization of e -> e weight matrix
W_ji = zeros(N);
W_jil = zeros(N,N,num_trials+1);
tim_to_puls = 9.72/1.2;
puls_to_tim = 2.5;%.6;
rec_eet = .75;
rec_eep = 0;
for i = 1:N
    if mod(i,2) == 1
        W_ji(i,i+1) = tim_to_puls;
        W_ji(i,i) = rec_eet;
    elseif mod(i,2) == 0
        W_ji(i,i) = rec_eep;
    end
end
W_jil_initial = W_ji;
%% initialization of i -> e weight matrix
M_ki = zeros(N);
M_kil = zeros(N,N,num_trials+1);
inh_tim_to_puls = 8.1/.5;
rec_iet = 2*.062;
rec_iep = 2*5;
for i = 1:N
    if mod(i,2) == 1
        M_ki(i,i+1) = inh_tim_to_puls;
        M_ki(i,1:2:N) = .2;
        M_ki(i,i) = rec_iet;
    elseif mod(i,2) == 0
        %         M_ki(:,i) = 10;
        M_ki(i,2:2:N) = 0;
        M_ki(i,i) = rec_iep;
    end
end
M_kil_initial = M_ki;
%% initialization of e -> i weight matrix
P_ik = .48*eye(N); %magnitude of synaptic weight from excitatory to inhibitory units
puls_to_puls_inh = 15;
% P_ik(12,1) = 15;
for i = 1:N
    if mod(i,2) == 1
        P_ik(i,i+1) = 0;
        P_ik(i,i) = .6;
    elseif mod(i,2) == 0
        P_ik(i,2:2:N) = puls_to_puls_inh;
        P_ik(i,i) = 0;
    end
end
% P_ik(10,12) = 1.2;
% P_ik(1,3) = .3;