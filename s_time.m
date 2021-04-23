%initialize variables
K = 100; %Number of active UE (P_active = 0.05 of N)
p = 0.05;
N = 2000; %Total number of UE
M = 16; %Number of antennas %%
L = 100; %Pilot sequence length %%
f = 10; %Bandwidth (MHz)
d_min = 100; %Cell radius (m)
d_max = 750; %Cell radius (m)
iter = 10; %Number of iterations

%figure 

MonteCarlo = 1000; %10 iter = 10000

active = false(N,MonteCarlo); 
signal = zeros(N,MonteCarlo);
time_set = zeros(3,1);
time = 0;

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
    index = randperm(N);
    active(index(1:K),j) = true;
    %active = false(N,1);
    %active(index(1:K)) = true;

    %Noise and power
    tx = 10^(23/10)*1e-3;
    noise_psd = 10^(-174/10)*1e-3;
    noise = noise_psd*f*1e6;
    noise_var = noise/(tx*L);
    noise_sd = sqrt(noise_var);

    w = noise_sd*(sqrt(0.5).*randn(L,M) + sqrt(0.5)*1i.*randn(L,M));

    x = zeros(N,M);
    x_sd = sqrt(N*p*(1 - p));
    for m = 1:M
        x(index(1:K),:) = h(:,index(1:K)).';
    end
    x = sqrt(L*tx).*x;

    y = a*x + w;
    
    %Choose the algorithm
    tic
    [xhat1, tau] = vectorAMP(K,N,M,L,y,a,pathloss,iter);
    %xhat2 = bomp(K,N,M,L,y,a,noise_sd,x_sd,20);
    time = time + toc;
end

avr_time = time/MonteCarlo/iter;

%algo = {'VAMP', 'BOMP'};
%num = [1 2];
%plot(num, time_set);
%xlabel('Algorithms');
%ylabel('Time per iteration');
%set(gca,'xtick',num,'xticklabel',algo);
%legend show;
