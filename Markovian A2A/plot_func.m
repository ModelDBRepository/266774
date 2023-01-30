function plot_func(T,delta,l,num_columns,t_j,t_total,dt,npp,N,R_it,R_kt,T_pt,T_dt,W_ji,Title,t_stim)
%info
new_rate = zeros(2*num_columns,numel(t_j));
new_rate_ele = 0;
new_rate2 = zeros(2*num_columns,numel(t_j));
new_rate_ele2 = 0;
pop_vect = 0:npp:N;
set(0,'DefaultAxesColorOrder',brewermap(12,'Paired')) 
for i = 1:2*num_columns
    for t = 1:((t_total)/dt)
        for k = pop_vect(i)+1:pop_vect(i+1)
            new_rate_ele = new_rate_ele + R_it(k,t);
            new_rate_ele2 = new_rate_ele2 + R_kt(k,t);
        end

        new_rate(i,t) = new_rate_ele/(npp);
        new_rate_ele = 0;
        new_rate2(i,t) = new_rate_ele2/(npp);
        new_rate_ele2 = 0;
    end
end


Plot_array = cell(4*num_columns,1);
for i = 1:(4*num_columns)
    if mod(i,2) == 0
        Plot_array{i,1} = new_rate(i/2,:);
    else
        Plot_array{i,1} = t_j;
    end
end

Plot_array2 = cell(4*num_columns,1);
for i = 1:(4*num_columns)
    if mod(i,2) == 0
        Plot_array2{i,1} = new_rate2(i/2,:);
    else
        Plot_array2{i,1} = t_j;
    end
end

if l == 1
    f = figure('rend','painters','pos',[100 100 600 600]);
    p = uipanel('Parent',f,'BorderType','none');
    p.BackgroundColor = [1 1 1];
    subplot(2,1,1,'Parent',p)
    plot(Plot_array{:,1},Plot_array2{:,1},'linewidth',4);
    axis([0.01 t_total 1 140]);
    title(Title);
    ylabel('Firing Rate (Hz)');
    xticks(t_stim-ones(1,length(t_stim)));
 


    subplot(2,1,2,'Parent',p)
    plot(t_j,T_dt, 'b--',t_j,T_pt, 'r','linewidth',4);
    axis([0.01 t_total 1 140]);
    xticks(t_stim-ones(1,length(t_stim)));
    xlabel('Time(ms)');
    ylabel('Eligibility Trace (AU)');


    p.FontSize = 34;
else
    subplot(2,1,1)
    plot(Plot_array{:,1},Plot_array2{:,1},'linewidth',4);
    axis([0.01 t_total 1 140]);
    title(Title);
    ylabel('Firing Rate (Hz)');
    xticks(t_stim-ones(1,length(t_stim)));
    
    subplot(2,1,2)
    plot(t_j,T_dt, 'b--',t_j,T_pt, 'r','linewidth',4);
    axis([0.01 t_total 1 140]);
    xticks(t_stim-ones(1,length(t_stim)));
    xlabel('Time(ms)');
    ylabel('Eligibility Trace (AU)');
end
drawnow
end


