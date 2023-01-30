function y = t_mik1(n,dt,T,t_total,p_r,npp,unit_stim,t_stim)
t_mik11 = zeros(n,t_total/dt);
ndd = npp;
for k = 1:length(unit_stim)
    j = (k-1)*ndd;
    t_mik11(j+1:j+ndd,(t_stim(k)+dt)/dt:(t_stim(k)+T)/dt) = (rand(ndd,T/dt) < p_r);
%     t_mik1(j+npp+1:j+2*npp,1:t_stim/dt) = (rand(npp,t_stim/dt) < .5*p_r);
end
y = t_mik11;
end