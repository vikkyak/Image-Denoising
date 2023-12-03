
%% TypeTwoFuzzyMembershipGeneration

function [T_min,T_max,T_min_max,PI,H,sigma,average_mu]=T2FM(vector)
p=numel(vector);
H=floor((p+1)/2);
mu=zeros(1,H);
for q=1:H;
    mu(1,q)=KMM(q,vector);
end
average_mu=((sum(mu))/H)*ones(1,p);
lembda= 4.5*abs(vector-average_mu);
mu_Mat=repmat(mu,p,1)';
lembda_vector=repmat(vector,H,1);
sigma=KMM(H,lembda);
if sigma < 0.0001
    sigma=0.0001;
end
PI= exp(-0.5*((lembda_vector-mu_Mat)/sigma).^2);
T_max=max(max(PI));
T_min_max=min(max(PI));
T_min=max(min(PI)); 
end










