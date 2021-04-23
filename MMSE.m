function [n_op, nprime_op] = MMSE(pathloss,tau,lambda,M,N,xhat)

    n_op = zeros(N,M);
    nprime_op = zeros(M,M);

    for n = 1:N

        a = pathloss(n)/(pathloss(n) + tau^2);
        b = (1 - lambda)/lambda*((pathloss(n) + tau^2)/tau^2)^M;
        c = pathloss(n)/(tau^2*(pathloss(n) + tau^2));

        d = 1 + b*exp(-c*(norm(xhat(n,:))^2));

        n_op(n,:) = a/d*xhat(n,:);
 
        nprime = a/d*eye(M) + a*b*c*(xhat(n,:)'*xhat(n,:))*exp(-c*norm(xhat(n,:))^2)/d^2;
        nprime_op = nprime_op + nprime;
    end

    nprime_op = nprime_op./N;
    
end

