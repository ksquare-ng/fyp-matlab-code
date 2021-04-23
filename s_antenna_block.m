clear;

%initialize variables
K = 100; %Number of active UE (P_active = 0.05 of N)
p = 0.05;
N = 2000; %Total number of UE
M_set = [1 2 4 8 16 32 64]; %Number of antennas %%
L = 100; %Pilot sequence length %%
f = 10; %Bandwidth (MHz)
d_min = 100; %Cell radius (m)
d_max = 750; %Cell radius (m)
iter = 10; %Number of AMP iterations
b = 10; %b-block sparsity

N_md = zeros(length(M_set),1); 
N_fa = zeros(length(M_set),1); 
P_md = zeros(length(M_set),1); 
P_fa = zeros(length(M_set),1); 

figure 

for i = 1:length(M_set)
    M = M_set(i);
    display(strcat("M = ", num2str(M)));
    
    MonteCarlo = 5000;
    
    active = false(N,MonteCarlo); 
    signal = zeros(N,MonteCarlo);
    
    for j = 1:MonteCarlo
        display(strcat(num2str(j), "-th iteration"));
        
        %Pilot sequence
        a = sqrt(0.5)*(1/L)^0.5.*randn(L,N) + sqrt(0.5)*(1/L)^0.5*1i.*randn(L,N);
    
        %Block fading model
        h = zeros(M,N);
        d = d_min + (d_max-d_min).*rand(N,1);
        pathloss = zeros(N,1);
        for n = 1:N
            pathloss(n) = 10^((-128.1 - 37.6*log10(d(n)*1e-3))/10);
            h(:,n) = sqrt(pathloss(n))*sqrt(1/2).*(randn(M,1) + 1i*randn(M,1)); %ignore shadowing for simplicity
        end

        %Active/Inactive
        block = N/b;
        index = randperm(block-1);
        for bidx = 1:(K/b)
            start = index(bidx)*b+1;
            last = start+b-1;
            active(start:last,j) = true;
        end
        
        %Noise and power
        tx = 10^(23/10)*1e-3;
        noise_psd = 10^(-174/10)*1e-3;
        noise = noise_psd*f*1e6;
        noise_var = noise/(tx*L);
        noise_sd = sqrt(noise_var);
        
        w = noise_sd*(sqrt(0.5).*randn(L,M) + sqrt(0.5)*1i.*randn(L,M));
        
        x = zeros(N,M);
        x_sd = sqrt(N*p*(1 - p));
        for n = 1:N
            if active(n,j)
                x(n,:) = h(:,n).';
            end
        end
        x = sqrt(L*tx).*x;
        
        y = a*x + w;
        
        %Choose the algorithm
        [xhat, tau] = vectorAMP(K,N,M,L,y,a,pathloss,iter);
        %xhat = bomp(K,N,M,L,y,a,noise_sd,x_sd,iter);
        
        for n = 1:N
            signal(n,j) = norm(xhat(n,:))^2;
        end
        
    end
    
    %Calculation for 50/50 FA/MD (VAMP)
    %{
    t_min = 0;
    t_max = 1;
    tolerance = 1e-15;
    while t_max - t_min >= tolerance
        t = (t_max + t_min)/2;
        N_md(i) = 0;
        N_fa(i) = 0;

        for j = 1:MonteCarlo
            for n=1:N
                if signal(n,j) < t && active(n,j) == 1
                    N_md(i) = N_md(i) + 1;
                end
                if signal(n,j) > t && active(n,j) == 0
                    N_fa(i) = N_fa(i) + 1;
                end
            end
        end

        P_md(i) = N_md(i)/(N*MonteCarlo);
        P_fa(i) = N_fa(i)/(N*MonteCarlo);

        if P_md(i) > P_fa(i)
            t_max = t;
        else
            t_min = t;
        end
    end
    %}
    
    %Calculation for FA/MD (BaGOMP)
    for j = 1:MonteCarlo
        for n=1:N
            if signal(n,j) == 0 && active(n,j) == 1
                N_md(i) = N_md(i) + 1;
            end
            if signal(n,j) > 0 && active(n,j) == 0
                N_fa(i) = N_fa(i) + 1;
            end
        end
    end
    
    P_md(i) = N_md(i)/(N*MonteCarlo);
    P_fa(i) = N_fa(i)/(N*MonteCarlo);
    
end

M=M_set;
semilogy(M,P_md,'ko-',M,P_fa,'k+-');
grid on
hold
xlabel('Antenna');
ylabel('Probability');
legend show;
