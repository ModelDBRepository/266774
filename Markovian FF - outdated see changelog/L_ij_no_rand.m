function y = L_ij_no_rand(num_columns,npp,N,l5_rec,l23_rec,l5_l23,l23_l5,l23_l23)
%instantiation of weight matrix without randomness
L_ij1 = zeros(N,N);
for i = 1:num_columns
    j = (i-1)*npp*2;
    if i < num_columns
        L_ij1(j+npp+1:j+(2*npp),j+(2*npp)+1:j+npp*3) =  l23_l5;
        for k = 0:2*npp:(num_columns - 1)*2*npp
            L_ij1(j+npp+1:j+2*npp,k+npp+1:k+2*npp) =  l23_l23;
        end
    end
    L_ij1(j+1:j+npp,j+1:j+npp) =  l5_rec;
    L_ij1(j+npp+1:j+2*npp,j+npp+1:j+2*npp) =  l23_rec;
    L_ij1(j+1:j+npp,j+npp+1:j+2*npp) =  l5_l23;
end
y = abs(transpose(L_ij1(1:N,1:N)));
end
