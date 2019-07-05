function L=yCoM(iter,P,C_para,f_para_even,f_dorsal_perp,spacing,C_perp,phi,h,f_ventral_perp)
L=zeros(iter,P/2);
% y=0 at CoM
for t=2:iter
    for i=1:P/2
        dLdt=(-1/C_para)*(2*f_para_even(t-1,i))*sin((pi/2)-phi(t-1,i))+(1/C_perp)*(f_dorsal_perp(t-1,i)+f_ventral_perp(t-1,i))*cos((pi/2)-phi(t-1,i));
        
        L(t,i)=L(t-1)+h*dLdt;
    end
end
end