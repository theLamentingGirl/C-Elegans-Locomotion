%?Is this a recursive func.? Since, S calls itself in the if statement. 
%is also ocalled by the total current input I_dorsal/ventral

%%---------------Function is for Neuronal state variable for n neural segments--------
function s=Sdors(iter,M,R1,L_dorsal_l,L0l,P)
N=12;
s=zeros(iter,N);

I_dorsal=Idorsal(iter,R1,L_dorsal_l,L0l,M,P);
I_ventral=Iventral(iter,R1,L_dorsal_l,L0l,M,P);
eps_hys=0.5;%hysteresis behaviour introducing state dep. activating/deactivatingthreshold
N_out=4;
for t=1:iter
    for m=1:M
        n=ceil(m/N_out);
        if I_dorsal(t,n) >0.5+eps_hys*(0.5-s(t-1,n))
            s(t,n)=1;
        elseif I_ventral(t,n) >0.5+eps_hys*(0.5-s(t-1,n))
            s(t,n)=1;
        else
            s=0;
        end
    end
end
end