

%% main program
for l = num_trials+1:num_trials+1
    %% initialization
    [v_yt,v_it,v_kt] = deal(zeros(N,t_total/dt)+v_rest); %membrane potentials
    [R_yt,R_it,R_kt,sc_R] = deal(zeros(N,t_total/dt)); %rates (sc_R is spikes)
    [s_yt,s_it,s_kt] = deal(zeros(N,t_total/dt)); %activations
    [g_Ey,g_Ei,g_Ii,g_Ek] = deal(zeros(N,t_total/dt)); %conductances
    [t_ref_y,t_ref_i,t_ref_k] = deal(zeros(N,t_total/dt)); %refractory periods
    [T_ijp,T_ijd,del_W_ji] = deal(zeros(N,N)); %synapse specific traces, weight update
    [T_pt,T_dt] = deal(zeros(1,t_total/dt + dt)); %mean trace for population at time t
    
    all_stim1 = t_mik1(n,dt,T,t_total,p_r,npp,unit_stim,t_stim); %neuron by neuron stimulation, all pops stimulated
    one_stim1 = t_mik1(n,dt,T,t_total,p_r,npp,1,t_stim); %neuron by neuron stimulation, first pop stimulated
    for i = 1:pop
        temp = (i-1)*10;
        all_stim(i,:,l) = mean(all_stim1(temp+1:temp+10,:),1); %population stimulation, all pops stimulated
        one_stim(i,:,l) = mean(one_stim1(temp+1:temp+10,:),1); %population stimulation, first pops stimulated
    end
 
    t_miy = one_stim1;

    %% time step loop
    for t = 2:((t_total)/dt) %at each time step
        %% pre-synaptic loop
        for i = 1:N %over each pre-synaptic neuron
            %% input layer dynamics
            if t_ref_y(i,t) == 1 %refractory period
                v_yt(i,t) = v_rest; %set voltage of neuron i to ressting potential if in refractory period
                t_miy(i,t) = 0; %no spike for neuron i at this time
            elseif (v_it(i,t-1) < v_th) && t_miy(i,t) == 0 %subthreshold update
                del_v_y = (g_L*(E_l-v_yt(i,t-1)))*(dt/C_m); %change in the membrane potential at each time step
                v_yt(i,t) = v_yt(i,t-1) + del_v_y; %update membrane potential
            elseif (v_yt(i,t-1) >= v_th) || t_miy(i,t) == 1
                v_yt(i,t) = v_hold; %voltage resets, neuron enter refractory phase
                t_miy(i,t) = 1; %spike for neuron i at time k
            end
            if v_yt(i,t) == v_hold %if neuron spikes
                del_R_yt = (1/dt-R_yt(i,t-1))*(dt/tau_w); 
                R_yt(i,t) = R_yt(i,t-1) + del_R_yt; %update firing rate
                del_s_y = -(s_yt(i,t-1)*dt/tau_si) + rho*(1-s_yt(i,t-1));
                s_yt(i,t) = s_yt(i,t-1) + del_s_y;
                t_ref_y(i,t:t+t_refractory) = 1;
            else %if neuron does not spike
                del_R_yt = -R_yt(i,t-1)*(dt/(tau_w));
                R_yt(i,t) = R_yt(i,t-1) + del_R_yt;
                del_s_y = -s_yt(i,t-1)*(dt/tau_si);
                s_yt(i,t) = s_yt(i,t-1) + del_s_y;
            end
            
            %% excitatory dynamics
            if t_ref_i(i,t) == 1 %refractory period
                v_it(i,t) = v_rest; %set voltage of neuron i to ressting potential if in refractory period
            elseif (v_it(i,t-1) < v_th) %subthreshold update
                del_v_i = ((randn/norm_noise)+g_L*(E_l-v_it(i,t-1)) + (g_Ei(i,t-1) + g_Ey(i,t-1))*(E_e - v_it(i,t-1)) + g_Ii(i,t-1)*(E_i - v_it(i,t-1)))*(dt/C_m);
                v_it(i,t) = v_it(i,t-1) + del_v_i; %update membrane potential
            elseif (v_it(i,t-1) >= v_th)
                v_it(i,t) = v_hold; %voltage resets, neuron enter refractory phase
            end
            if v_it(i,t) == v_hold %if neuron spikes
                sc_R(i,t) = 1;
                del_R_it = (1/dt-R_it(i,t-1))*(dt/tau_w); 
                R_it(i,t) = R_it(i,t-1) + del_R_it; %update firing rate
                del_s_j = -(s_it(i,t-1)*dt/tau_se) + rho*(1-s_it(i,t-1));
                s_it(i,t) = s_it(i,t-1) + del_s_j;
                t_ref_i(i,t:t+t_refractory) = 1;
            else %if neuron does not spike
                del_R_it = -R_it(i,t-1)*(dt/tau_w);
                R_it(i,t) = R_it(i,t-1) + del_R_it;
                del_s_j = -s_it(i,t-1)*(dt/tau_se);
                s_it(i,t) = s_it(i,t-1) + del_s_j;
            end
            
            %% inhibitory dynamics
            if t_ref_k(i,t) == 1 %refractory period
                v_kt(i,t) = v_rest; %set voltage of neuron i to ressting potential if in refractory period
            elseif (v_kt(i,t-1) < v_th_i) %subthreshold update
                del_v_k = ((randn/norm_noise)+ g_L*(E_l-v_kt(i,t-1)) + (g_Ek(i,t-1) + (iG/eG)*g_Ey(i,t-1))*(E_e - v_kt(i,t-1)))*(dt/C_m); 
                v_kt(i,t) = v_kt(i,t-1) + del_v_k; %update membrane potential
            elseif (v_kt(i,t-1) >= v_th_i) 
                v_kt(i,t) = v_hold; %voltage resets, neuron enter refractory phase
            end
            if v_kt(i,t) == v_hold %if neuron spikes
                del_R_kt = (1/dt-R_kt(i,t-1))*(dt/tau_w); 
                R_kt(i,t) = R_kt(i,t-1) + del_R_kt; %update firing rate
                del_s_k = -(s_kt(i,t-1)*dt/tau_si) + rho*(1-s_kt(i,t-1));
                s_kt(i,t) = s_kt(i,t-1) + del_s_k;
                t_ref_k(i,t:t+t_refractory) = 1;
            else %if neuron does not spike
                del_R_kt = -R_kt(i,t-1)*(dt/tau_w);
                R_kt(i,t) = R_kt(i,t-1) + del_R_kt;
                del_s_k = -s_kt(i,t-1)*(dt/tau_si);
                s_kt(i,t) = s_kt(i,t-1) + del_s_k;
            end
            
        end  
        %% updating conductances
        g_Ey(:,t) = W_in(:,:)*s_yt(:,t); %input conductance
        g_Ei(:,t) = W_ji(:,:)*s_it(:,t); %recurrent excitatory conductance
        g_Ii(:,t) = M_ki(:,:)*s_kt(:,t); %I to E conductance
        g_Ek(:,t) = P_ik(:,:)*s_it(:,t); %E to I conductance
        

        T_pt(t) = mean(T_ijp(2*npp+1:3*npp,npp+1:2*npp),'all')*100000; %ff traces for graphing
%         T_pt(t) = mean(T_ijp(1:npp,1:npp),'all')*100000; %recurrent
%         traces for graphing
        T_dt(t) = mean(T_ijd(2*npp+1:3*npp,npp+1:2*npp),'all')*100000; %ff traces for graphing
%         T_dt(t) = mean(T_ijd(1:npp,1:npp),'all')*100000; %recurrent
%         traces for graphing

    end
    
    
    %rescaling firing rates
    R_it = R_it*1000;
    R_kt = R_kt*1000;
    R_yt = R_yt*1000;


    
    
    sprintf('Trial %d complete',l)
    %plotting
    if l == 1
        old_R_it = plot_R_it(:,:,l);
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,R_it,0*R_kt,0*T_pt,0*T_dt,W_ji,'Before Learning',t_stim)
    elseif l == num_trials
        tit = sprintf('Final Learning Trial %d (Full Stim)',num_trials);
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,R_it,0*R_kt,0*T_pt,0*T_dt,W_ji,tit,t_stim)
    elseif l == num_trials+1
        new_R_it = plot_R_it(:,:,l);
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,R_it,0*R_kt,T_pt,T_dt,W_ji,'After Learning (one stim)',t_stim)
    else
        tit = sprintf('During Learning Trial %d',l);
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,R_it,0*R_kt,T_pt,T_dt,W_ji,tit,t_stim)
    end
drawnow    
end
