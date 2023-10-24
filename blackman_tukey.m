function Pbt = blackman_tukey(x, f, window, M)
    % Prozorske funkcije mogu biti Bartletova ili Parzenova
    % jer one imaju strogo pozitivne spektre

    N = length(x);
    Nf = length(f);
    Pbt = zeros(1, Nf);

    if window == "bartlet"
        w = bartlett(M)';
    elseif window == "parzen"
        w = parzenwin(M)';
    end
    
    w = [zeros(1,(2*N-1-M)/2) w zeros(1,(2*N-1-M)/2)];
   

    rxx = autocorrelation(x);
    rxx_win = rxx.*w;

    for i = 1:Nf
        for j = 1:2*N-1
            Pbt(i) = Pbt(i) + rxx_win(j)*exp(-1i*2*pi*f(i)*(j-N));
        end
    end
    
    Pbt = real(Pbt);

end