clear all
close all
%% Parameters and Initialization
run('parameters_rate_final.m')

%% initialization of weight matrices
run('weight_initialization_final.m')

%% Stimulation Matrix
run('stim_init_final.m')

u_it0 = normrnd(0,1.4,M,t_total/dt); %LSM inital activation
u_it0copy = u_it0;
v_t0 = normrnd(0,0.012,4*M,t_total/dt); %sparse net initial activation

%% Main Program 
for l = 1:num_trials + 1
    trace_refractory = zeros(N,N,numel(time_t));
    if (l == num_trials + 1) || (l == 1)
        I_ext_it = one_stim;
        I_ext_it_inh = one_stim_inh;
    else
        I_ext_it = all_stim;
        I_ext_it_inh = all_stim_inh;
    end
    v_t1 =  v_t0;
    v_t2 =  zeros(4*M,t_total/dt);
    rate_vt1 =  zeros(4*M,t_total/dt);
    v_it_01 = zeros(num_columns,t_total/dt); %input to LSM from messenger cells
    u_it1 = u_it0;
    rate_it1 = zeros(M,t_total/dt);
    W_jil(:,:,l) = W_ji;
    M_kil(:,:,l) = M_ki;
    [u_it, v_kt, rate_uit, rate_vkt] = deal(zeros(N,t_total));
    [del_ui, del_vk] = deal(zeros(N,1));
    
    for t = 2:t_total/dt - (1/dt)
        for y = 1:num_columns
            u = (y-1)*2 + 2;
            v_it_01(y,t) = rate_uit(u,t-1); %input to LSM from messenger cells, pop level
        end
        v_it_2 = repelem(v_it_01,M/num_columns,1); %input to LSM from messenger cells, over all neurons
        v_it_1 = theta2(v_it_2,.75); %input to LSM from messenger cells, over all neurons, thresholded
        sum_wj = W_ji2(:,:)*phi(u_it1(:,t-1),r_0); %recurrent activation 
        del_ui1 = (2*v_it_1(:,t) -u_it1(:,t-1) + sum_wj)*(dt/tau_net); %change in activation for LSM units
        u_it1(:,t) = u_it1(:,t-1) + del_ui1;
        rate_it1(:,t) = phi(u_it1(:,t),r_0); %LSM activation converted to rates
        sum_oj = O_ij(:,:)*phi(u_it1(:,t-1),r_0); %LSM to sparse net activation
        del_v = (-v_t1(:,t-1) + sum_oj)*(dt/(2*tau_net)); %change in sparse net activation
        v_t1(:,t) = v_t1(:,t-1) + del_v;
        v_toty = repelem(v_it_1(:,t),4,1);
        v_t2(:,t) = v_t1(:,t-1).*theta2(v_toty,0.01); %sparse net activation which overlaps with messenger activation
        rate_vt1(:,t) = theta2(v_toty,0.01).*theta_ish(v_t1(:,t),theta_c,theta_c2); %sparse net rates (binary)
        
        for i = 1:N
            sum_uj = 0;
            sum_vk = 0;
            sum_ui = 0;
            sum_pj = 0;
            sum_qj = 0;
            for j = 1:N
                sum_uj = sum_uj + W_ji(j,i)*(phi_2(u_it(j,t-1),v,theta,u_c)); % e to e activation, modular net
                sum_vk = sum_vk + M_ki(j,i)*(phi_2(v_kt(j,t-1),v,theta,u_c)); %i to e activation, modular net
                sum_ui = sum_ui + P_ik(j,i)*(phi_2(u_it(j,t-1),v,theta,u_c)); %e to i activation, modular net
            end
            for k = 1:4*M
                sum_qj = sum_qj + Q_ij(k,i)*rate_vt1(k,t-1); %sparse net to Timer activation
            end
            del_ui(i) = (I_ext_it(i,t) - u_it(i,t-1) + sum_uj-sum_vk + sum_qj)*(dt/tau_u); % change in activation, exciatory in modular net
            del_vk(i) = ((I_ext_it_inh(i,t)*scale_v) - v_kt(i,t-1) + sum_ui)*(dt/tau_v); % change in activation, inhibitory in modular net
        end
        u_it(:,t) = u_it(:,t-1) + del_ui;
        v_kt(:,t) = v_kt(:,t-1) + del_vk;
        for i = 1:N
            rate_uit(i,t) = phi_2(u_it(i,t),v,theta,u_c);
            rate_vkt(i,t) = phi_2(v_kt(i,t),v,theta,u_c);
            for j = 1:N
                if (i==j && mod(j,2) == 1) && (t > D) && (l < num_trials + 1) && (t< t_total/dt) %only recurrent learning in modular net
                        T_max_p = T_max_p1;
                        T_max_d = T_max_d1;
                        tau_p = tau_p1;
                        tau_d = tau_d1;
                        eta_d = eta_d1;
                        eta_p = eta_p1;
                    if trace_refractory(i,j,t) == 0
                        if rate_uit(i,t-D) > .01 && rate_uit(j,t-1) > .01
                            H_d = eta_d*rate_uit(i,t-D)*rate_uit(j,t-1); %Hebbian learning term for depression
                            H_p = eta_p*rate_uit(i,t-D)*rate_uit(j,t-1); %Hebbian learning term for potentiation
                        else
                            H_d = 0; %Hebbian learning term for depression
                            H_p = 0; %Hebbian learning term for potentiation
                        end
                        del_T_ijp = (-T_ijpt(i,j,t-1) + (H_p*(T_max_p - T_ijpt(i,j,t-1))))*(dt/tau_p); %change in LTP eligibility trace
                        del_T_ijd = (-T_ijdt(i,j,t-1) + (H_d*(T_max_d - T_ijdt(i,j,t-1))))*(dt/tau_d); %change in LTD eligibility trace
                        T_ijpt(i,j,t) = T_ijpt(i,j,t-1) + del_T_ijp; %update LTP eligibility trace
                        T_ijdt(i,j,t) = T_ijdt(i,j,t-1) + del_T_ijd; %update LTD eligibility trace
                    else
                        T_ijpt(i,j,t) = 0;%T_ijpt(i,j,t)/2;
                        T_ijdt(i,j,t) = 0;%T_ijdt(i,j,t)/2;
                    end
                    for m = t_reward
                        del_time = round(delay_time/4);
                        if t/(m+ 5*del_time) == 1 %at time of reward
                            trace_refractory(i,j,t+1:t+75/dt) = 1;
                        end
                        if t > (m + 3*del_time) && t < (m+ 5*del_time)
                            if l>1
                                del_W_ji = eta_rec*(T_ijpt(i,j,t)-T_ijdt(i,j,t))*(2*dt/delay_time); %change in synaptic weights at reward time
                            else
                                del_W_ji = 0;
                            end
                            if (W_ji(i,j) + del_W_ji) >= 0 && trace_refractory(i,j,t) == 0 
                                W_ji(i,j) = W_ji(i,j) + del_W_ji; %update weights
                            end
                        end 
                    end
                end
            end
        end
        if t > 150/dt
            H_ij1 = rate_uit(:,t)*rate_vt1(:,t-150)';
            %the following is a diagonal mapping from modular timers to
            %sparse net
            a = 1120/5;
            H_ij1(2:2:N,:) = 0; %messengers do not participate
            H_ij1(1,1:a) = 0; %first timer connects to first "a" sparse net cells
            H_ij1(3,a+1:2*a) = 0; %second timer connects to a+1 to 2a sparse net cells
            H_ij1(5,2*a+1:3*a) = 0; %etc.
            H_ij1(7,3*a+1:4*a) = 0;
            H_ij1(9,4*a+1:5*a) = 0;
            H_ij1 = H_ij1';
            del_Q_ij = del_Q_ij + H_ij1; %change in sparse net -> modular Timer weights
        end
    end
    del_Q_ij = del_Q_ij*(dt/tau_q);
    Q_ij = Q_ij + del_Q_ij;
    Q_ij(Q_ij<0) = 0;
    Q_ij(Q_ij>0.15) = .15; %upper limit on sparse -> timer weights
    Q_ijl(:,:,l) = Q_ij;
    del_Q_ij = 0;
    if l == 1
        old_rate_uit = rate_uit;
        old_rate_vkt = rate_vkt;
    end
    if ismember(l,l_array)
        plot_rate_uit(:,:,l) = rate_uit;
        plot_rate_vkt(:,:,l) = rate_vkt;
        LTP0(:,l) = 40*T_ijpt(1,1,:);
        LTD0(:,l) = 40*T_ijdt(1,1,:);
        plot_rate_vt1(:,:,l) = rate_vt1;
    end
    sprintf('Trial %d complete',l)
    
end
%end main program

slider_graph_microcircuit_final(rate_uit,rate_vkt,plot_rate_uit, plot_rate_vkt, LTP0, LTD0,...
    num_trials, time_t,l_array,t_total,W_jil,one_stim,all_stim,old_rate_uit,old_rate_vkt,W_jil_initial,T,delta,dt,plot_rate_vt1,Q_ijl,t_stim);



function y = theta_ish(u,theta_c,theta_c2)
y = (theta_c2>u) & (u>theta_c);
end

function y = theta2(u,theta_c)
y = u>theta_c;
end

function y = phi(u,r_0)
y = r_0*tanh(u/r_0).*(u<=0) + (2-r_0)*tanh(u/(2-r_0)).*(u>0);
end