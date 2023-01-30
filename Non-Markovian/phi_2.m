function y = phi_2(u,v,theta,u_c)
if isnan(u) == 1
    y = 0;
elseif theta >= u
    y = 0;
elseif (theta < u) && (u < u_c)
    y = v*(((u-theta)/(u_c-theta))^2);
elseif u_c <= u 
    y = 2*v*sqrt(((u-theta)/(u_c - theta))-.75);
end
        
end

%%
%activation transfer function phi

% function y = phi(u,v,theta,u_c)
% if theta > u
%     y = 0;
% elseif (theta <= u) && (u <= u_c)
%     y = v*(u-theta);
% elseif u_c < u 
%     y = v*(u_c - theta);
% end
%         
% end