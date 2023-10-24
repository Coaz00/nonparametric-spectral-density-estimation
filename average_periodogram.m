function Pavper = average_periodogram(x,f,K)
    
    N = length(x);
    Nf = length(f);
    L = floor(N/K);

    Pavper = zeros(1,Nf);
    for i = 1:K-1
        Pavper = Pavper + periodogram(x((i-1)*L+1:i*L), f);
    end

    % ostatak
    Pavper = Pavper + periodogram(x((K-1)*L + 1:end), f);

    Pavper = 1/K*Pavper;

end