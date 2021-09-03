function A=create_A(EE,EI,IE,II,Ne,Ni)
% rng('shuffle')

% Probability of connection to each neuron is between 0 and 1
% based on number of Ne and Ni each neuron recieve specefic number of
% connections.
% 1st input:    excitatory->excitatory
% 2nd input:    excitatory->inhibitory
% 3rd input:    inhibitory->excitatory
% 4th input:    inhibitory->inhibitory
% 5th input:    Number of excitatory neuron(s)
% 6th input:    Number of inhibitory neuron(s)

probab_exci_exci=round(EE*Ne);  
probab_exci_inhi=round(EI*Ne);
probab_inhi_exci=round(IE*Ni);
probab_inhi_inhi=round(II*Ni);
N=Ne+Ni;
A=zeros(Ne+Ni);

for ii=1:Ne % ii is Post
    X=randperm(Ne,probab_exci_exci); % X is Pre
    A(ii,X)=1;
    Y=randperm(Ni,probab_inhi_exci); %Y is Pre
    A(ii,Ne+Y)=1;
end
for ii=Ne+1:Ni+Ne %ii is Post
    X=randperm(Ne,probab_exci_inhi); %X is Pre
    A(ii,X)=1;
    Y=randperm(Ni,probab_inhi_inhi); %Y is Pre
    A(ii,Ne+Y)=1;
end
C=eye(N,N);
A=A.*(~C);
end