function xhat = bomp(K,N,M,L,y,a,sigma_w,sigma_x,iter) 

%Moore-Penrose, how it links to standard SR
%BBOMP, group sparsity, line 21 use average 

    r = y;
    xhat = zeros(N,M);
    x = zeros(N,M);
    shat = zeros(N,1);
    s = zeros(N,1);
    temp = zeros(M,1);
    f = zeros(N,1);
    lambda = K/N;
    j = 1;
    
    for t = 1:iter
        
        for n = 1:N
            
            T = 2*sigma_w^2*(sigma_w^2 + sigma_x^2)/sigma_x^2*log((1-lambda)/lambda);
            
            %for m = 1:M
            %    temp(m) = (r(:,m) + a(:,n)*xhat(n,m))'*a(:,n);
            %end
            
            temp = (r + a(:,n)*xhat(n,:))'*a(:,n);
            
            if ((sum(temp)/M)^2 > T)
                s(n) = 1;
            else
                s(n) = 0;
            end
           
            %for m = 1:M
            %    x(n,m) = s(n)*sigma_x^2/(sigma_x^2 + sigma_w^2)*(xhat(n,m)+a(:,n)'*r(:,m));
            %end
            
            x(n,:) = s(j)*sigma_x^2/(sigma_x^2 + sigma_w^2)*(xhat(n,:)+a(:,n)'*r);
            
            t1 = sigma_w^2/sigma_x^2*norm(x(n,:))^2;
            t2 = sigma_w^2*log10((1 - lambda)/lambda)*s(n);
            t3 = r + a(:,n)*(shat(n).*xhat(n,:) - s(n).*x(n,:));
            t3 = abs(t3).^2;
            t3 = sum(t3,[1 2]);
            
            f(n) = t1 - t2 - t3;
        end
        
        [~,j] = max(f);
        shat(j) = s(j);
        
        idx = find(shat);
        a_s = a(:,idx);
        
        xhat_temp = (inv(a_s'*a_s + sigma_w^2/sigma_x^2*eye(length(idx)))*a_s')*y;
        
        i = 1;
        for n = 1:N
            if i <= length(idx) && n == idx(i)
                xhat(n,:) = xhat_temp(i,:);
                i = i+1;
            else
                xhat(n,:) = 0;
            end
        end
        
        r = y - a*xhat;
    end
    
end

