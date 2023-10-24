function rxx = autocorrelation(x)

    N = length(x);

    rxx = zeros(1,2*N-1);

    for i = 0:N-1
        for m = 1:N-i
            rxx(i+N) = rxx(i+N) + x(m)'*x(m+i);
        end
    end
    
    for i = 1:N-1
        rxx(i) = rxx(2*N - i)';
    end

    rxx = 1/N*rxx;

end