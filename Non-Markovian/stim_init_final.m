I_ext_it = zeros(N,t_total); %initiation of external input matrix
unit_stim = [1 3 1 5 1 7];
% t_stim = 1:(T+delta):t_total;
t_stim = [1 501 1201 1701 3001 3501 4201];
t_reward = t_stim;
for i = 1:numel(unit_stim)
    s = randn(1)/2;
    I_ext_it(unit_stim(i),t_stim(i):t_stim(i) + T) = I_strength;
end


all_stim = I_ext_it;
one_stim = zeros(N,t_total);
one_stim(unit_stim(1),:) =  all_stim(unit_stim(1),:);
one_stim(:,500:t_total) = 0;
one_stim_inh = one_stim;
all_stim_inh = all_stim;     