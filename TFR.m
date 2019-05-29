function [Tk] = TFR(x,Rk,Wk,nk)
    m = 5;  % Orden Polinomio
    X = reshape(x,[m,m]);
    Tk = zeros(length(Rk),length(Wk));
    T = zeros(m,m);
    for rk = Rk
        for wk = 1:length(Wk)
            for i = 1:m
                for j = 1:m
                    T(i,j) = X(i,j)*nk*Wk(wk)^i*Rk(rk)^j;
                end 
            end
            Tk(rk,wk) = sum(sum(T));
        end
    end 
end

