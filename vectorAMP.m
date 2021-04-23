function [xhat,f_tau] = vectorAMP(K,N,M,L,y,a,pathloss,iter)

    r = y;
    x = zeros(N,M);
    tau = zeros(iter,1); % 1 -> iter
    
    for t = 1:iter
     
        xhat = a'*r + x;
        avg_rnorm = 0;
        for m = 1:M
            avg_rnorm = avg_rnorm + norm(r(:,m))^2;
        end
        tau(t) = sqrt(avg_rnorm/M/L);
        [x, nprime] = MMSE(pathloss,tau(t),K/N,M,N,xhat);
        %[x, nprime] = threshPrimeThreshComplexGaussian(xhat,N,M,tau(t),K/N,pathloss);
        r = y - a*x + N/L.*r*nprime;
        
    end
    avg_rnorm = 0;
    for m = 1:M
        avg_rnorm = avg_rnorm + norm(r(:,m))^2;
    end
    f_tau = sqrt(avg_rnorm/M/L);
    %f_tau = f_tau(end);
    
end

