function A=activation_ventral(h,tau,iter,M)

I_ven_NMJ=I_ventral_NMJ(iter,M);

A=zeros(iter,M);

for m=1:M
    for t=2:iter
    
    dAdt=(1/tau)*(I_ven_NMJ(t-1,m)-dA(t-1,m));
    A(t,m)=A(t-1,m)+h*dAdt;

    end
end
end