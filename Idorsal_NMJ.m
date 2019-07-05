function I=Idorsal_NMJ(M,n)

I_dorsal=Idorsal(iter,R1,L_dorsal_l,L0l,M,P);
s=S(iter,M,R1,L_dorsal_l,L0l,P);


%N_out=4;
I=zeros(n,M);
w_NMJ=1;
w_NMJ_bar=-1;

%at t=0,I->for all 48 seg  = 0
for m=1:M
    for t=2:iter
    
    %w_NMJ->Excitatory neuromuscular jn weights
    %w_NMJ'->GABAergic neuromuscular weights(Inhibitory)
    I(t,m)=w_NMJ*s(t,n)+w_NMJ_bar*s(t,n);
    end
end++