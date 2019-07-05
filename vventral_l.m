
%----------lateral passive muscles, Ventral------------------
function v=vventral_l(h,iter,M,L_ventral_l)
    v=zeros(iter,M);
    %how length changes with time
    %at t=0 v=0
for i=2:P/2
    for t=2:iter
    
        dLdt=L_ventral_l(t-1,i-1);
        v(t,i-1)=v(t-1,i-1)+h*dLdt;
    end
end
end