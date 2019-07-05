 
%-----------diagonal passive muscles, velocity----------------------
function v=vventral_d(h,iter,M,L_ventral_d)
    v=zeros(iter,M);
    %how length changes with time
    %at t=0 v=0
for i=2:P/2
    
    for t=2:iter
    
        dLdt=L_ventral_d(t-1,i-1);
        v(t,i-1)=v(t-1,i-1)+h*dLdt;
    end
    
end
end