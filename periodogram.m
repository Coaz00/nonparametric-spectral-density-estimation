function Pper = periodogram(x, f)

    N = length(x);
    Nf = length(f);
    Pper = zeros(1,Nf);

    for i = 1:Nf
        for j = 1:N
            Pper(i) = Pper(i) + x(j)*exp(-1i*2*pi*f(i)*(j-1));
        end
    end
    
    Pper = abs(Pper).^2/N;
end