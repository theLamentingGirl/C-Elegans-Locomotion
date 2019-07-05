function dAdt=activation_dorsal(t,A,tau);
    I_dors_NMJ=Idorsal_NMJ(M);
%M different ode
    dAdt=zeros(M,1);
    for m=1:M
    
        dAdt(m)=(1/tau)*(I_dors_NMJ(m)-A(m));

    end

end