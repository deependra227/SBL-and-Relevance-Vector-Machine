% EE698Z: Machine Learning for Wireless Communications
% MATLAB ASSIGNMENT : SBL and Relevance Vector Machine
%180227-Deependra Chansoliya
%180589-Ramandeep
%180621-Rohan Mahnot

M=40;
N=20;
Do=7;
density = Do/M;

noise = [-20, -15, -10, -5, 0];% dB matrix
vars = 10.^(noise/10);% noise variance matrix
NMSE = zeros(5,1);

for itr = 1:2000
    nmse = zeros(5);
    
    for j = 1:5
        
        phi = randn(N,M);%Qno 2 part-1
        w = sprandn(M,1, density);%Qno 2 part-2
        mu = ones(M,1);
        A = diag(100* ones(M,1)); 
        t = phi*w + sqrt(vars(j)).*randn(N,1);%Qno 2 part-3
        %updating parameters
        while 1
            
            Sig = inv((1/vars(j)) * phi' * phi + A);% eqn-12
            mu_new = (1/vars(j)) * Sig * phi' * t;% eqn-13
            
            if (norm(mu-mu_new)^2)/(norm(mu)^2) <= 1e-3
                mu = mu_new;
                break   
            else
                mu = mu_new;
            end
            G = 1 - diag(A) .* diag(Sig);% gamma
            A = diag(G./(mu.^2 +eps));% alpha
            
        end
        
        Wmp = mu;
        nmse(j) = (norm(Wmp-w)^2/norm(w)^2);
    end
    
    NMSE= NMSE + nmse; 
end

NMSE = NMSE ./ itr;% Average value after itr iterations

semilogy(noise,NMSE,'-o');%Qno 5 plot of NMSE
xlabel('Noise Variance (in dB)');
ylabel('Average NMSE');
title('NMSE vs Noise Variance');