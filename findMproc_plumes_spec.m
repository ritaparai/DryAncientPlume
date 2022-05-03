function output = findMproc_plumes_spec(late_t,startidx,Nres,Mres)

tx = late_t(startidx-1:end) - late_t(startidx-1);
T = tx(end);
Qp = 2.33 * 2.9e15 * 10; 

Talpha = fzero(@(x) findTalpha(x,Mres,Nres,T,Qp),0.5);


%% MAKE MPROC BLOCK ALONG TIME AXIS
tp = T-tx;
Mproc = zeros(size(Talpha,1), size(Talpha,2), numel(tx)-1);
for n = 2:numel(tx)
    tnow = tp(n);
    tlast = tp(n-1);
    Mproc(:,:,n-1) = (Qp.*T./Talpha) .* (exp(Talpha.*tlast./T) - exp(Talpha.*tnow./T));
end

output = Mproc;
%%MPROC is a cheeseblock that is Nres on dim1, Ms on dim2, time on dim3
