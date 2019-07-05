function I=Iventral(iter,R1,L_ventral_l,L0l,M,P)
N=12;%-> no. of neural units
N_out=4;%->M/N
N_SR=M/2;

I_ventral_AVB=zeros(iter,N);%const. input current in forward locomotion 


G_SR=zeros(iter,N);%initialsing vector for conductance parameter;
A=zeros(iter,N);%initialising vector for weighing prefactor 
lambda=zeros(iter,M);
h_ventral=zeros(iter,M);
h_ventral_sum=zeros(iter,N);
for t=1:iter
    for i=2:P/2
        %stretch receptor current
        %n is a function of m
        n=ceil(i-1/N_out);
        
        
        lambda(t,i)=2*R1(i-1)+R1(i);
        
        for o=1+(n-1)*N_out:min([M,N_SR+(n-1)*N_out])
            %2nd weighting term
            if L_ventral_l(t,i-1)>L0l(t,i-1)
                gamma_ventral=0.8;
            else
                gamma_ventral=1.2;
            end
            %effective mechanosensory activation fn
            h_ventral(t,o)=lambda(t,o)*gamma_ventral*((L_ventral_l(t,o)-L0l(t,o))/L0l(t,o));
            
            h_ventral_sum(iter,n)=h_ventral(t,o)+h_ventral(t,o+1);
            
        end
        
    %conductance parameter; inc. linearly H->t t ocompensate for decresing
    %curvature of undulatns 
        G_SR(t,n)=(0.224+0.056*n)/N_out;
     
    %weighing prefactor A
        if (n-1)*N_out<=M-N_SR
            A(t,n)=1;
        else
            A(t,n)=sqrt(N_SR/(M-(n-1)*N_out));
        end
        
    I_ventral_SR=A(t,n)*G_SR(t,n)+h_ventral_sum;
    end
end
w_ventral=-1;
%Total I_ventral output
%Here S is neuronal state variable defined by another fn
I=I_ventral_AVB(iter,n)+I_ventral_SR(iter,n)+w_ventral*S(iter,n);

end