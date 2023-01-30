function y = input_weights(num_columns,npp,N,n,gain)
%instantiation of weight matrix
L_ij1 = zeros(N,N);
ndd = npp;
for i = 1:num_columns
    j = (i-1)*npp*2;
    l = (i-1)*ndd;
    L_ij1(j+1:j+npp,l+1:l+ndd) =  gain*eye(npp);
end
y = abs(L_ij1);
end
