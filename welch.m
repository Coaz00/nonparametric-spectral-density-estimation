function Pwelch = welch(x, f, K, p)

    N = length(x);
    Nf = length(f);
    Pwelch = zeros(1,Nf);

    L = floor(N/(1 + (1-p)*(K-1)));

    % prvih L pa shift pocetak i kraj svaki put
    shift = floor(L*(1-p));
    
    for i = 0:K-2
        Pwelch = Pwelch + periodogram(x(1+i*shift:L+i*shift).*hamming(L)', f);
    end

    %ostatak
    Pwelch = Pwelch + periodogram(x(1+(K-1)*shift:end).*hamming(length(x(1+(K-1)*shift:end)))',f);

    Pwelch = 1/K*Pwelch;
    
end