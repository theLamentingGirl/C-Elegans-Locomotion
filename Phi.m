function omega=Phi(iter,P,C_para,f_para_odd,R1,h)
    omega=zeros(iter,P/2);
    omega(1,:)=90;
for t=2:iter
    for i=1:P/2
        dphidt=(1/R1(i)*C_para)*(2*f_para_odd(t-1,i));
        omega(t,i)=omega(t-1,i)+h*dphidt;
    end
end
end