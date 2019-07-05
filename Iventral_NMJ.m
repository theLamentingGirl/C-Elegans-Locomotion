function I=Iventral_NMJ(n,iter,M)

s=S(iter,n,I_dorsal,I_ventral);
%N_out=4;
I=zeros(n,M);
w_NMJ=1;
w_NMJ_bar=-1;

for m=1:M
    for t=2:iter
        
    %w_NMJ->Excitatory neuromuscular jn weights
    %w_NMJ'->GABAergic neuromuscular weights(Inhibitory)
    I(t,m)=w_NMJ*s(t,n)+w_NMJ_bar*s(t,n);
    end
end