%------------dimensions of the worm------------------
R=40*10^-6;%Minor axis radius of the worm
M=48;%no. of segments
P=2*(M+1);%-->no. of points (pk & pk_bar)=98
L=(1*10^-3);%-->Length of the worm = 1mm
ti=0;
tf=10;
iter=100;%->Time points
h=tf-ti/iter;
%equally spaced 48 segments
spacing=linspace(0,L,P/2);

%-----------varying the minor radius-------------------
R1=zeros(1,P/2);
for i=1:P/2
    R1(i)=R*sin(acos((i-M/2)/((M+1)/2+0.2)));
end

%-----------plotting the shape of the worm as described by the R1 eqn->prolate ellipse

plot(spacing,R1,'k')
hold on
plot(spacing,-R1,'k')
xlim([0 2*10^-3])
ylim([-20*10^-5 20*10^-5])
%count=0;%->to check if 48 segments present

%-------Representing the lateral and diagonal muscles through joining pd and
%pv points-----------------------------
for i=1:P/2
    if (i+1)==P/2
        break
    else
    plot([spacing(i) spacing(i)],[-R1(i) R1(i)],'b')
    plot([spacing(i) spacing(i+1)],[-R1(i) R1(i+1)])
    
    plot([spacing(i) spacing(i+1)],[R1(i) -R1(i+1)])
    
    end
    %count=count+1;
end
%%
%-------------------------properties of the environment--------------------
N=12;
C_para=3.2*10^-3;
C_perp=128*10^-3;


%========================Equations of motion========================================
%initialising vectors

f_dorsal=zeros(iter,P/2);
f_ventral=zeros(iter,P/2);
f=zeros(iter,P/2);
f_dorsal_para=zeros(iter,P/2);
f_dorsal_perp=zeros(iter,P/2);
f_ventral_para=zeros(iter,P/2);
f_ventral_perp=zeros(iter,P/2);
f_para_even=zeros(iter,P/2);
f_para_odd=zeros(iter,P/2);


%%%-=-=-=-=-=-=-=-=-=-We need to describe the postion x,y, phi of 49 rods
%%%over time period (tf-ti)
%-------------renaming/storing pts to x & y coordinates-----------
%initial conditions for the ode of eqns of motion

phi=Phi(iter,P,C_para,f_para_odd,R1,h);
x_pd=zeros(iter,P/2);
x_pv=zeros(iter,P/2);
y_pd=zeros(iter,P/2);
y_pv=zeros(iter,P/2);
x_CoM=xCoM(iter,P,C_para,f_para_even,f_dorsal_perp,spacing,C_perp,phi,h,f_ventral_perp);
y_CoM=yCoM(iter,P,C_para,f_para_even,f_dorsal_perp,spacing,C_perp,phi,h,f_ventral_perp);

L_CoM=zeros(iter,P/2);

%%========================Body forces===========================
%--------------------------parameters---------------------------
L_seg=L/M;
kl=(M/24)*0.01;
kd=kl*350;
beta_l=kl*0.025;
beta_d=kd*0.01;
beta_0M=kd*100;
k0M=kl*20;
del_m=0.75;%-->parameter part of neural circuitry

%%d->diagonal;l->lateral elements
%initialising vectors to hold varying lengths b/w pi & pi+1
L_dorsal_l=zeros(iter,M);
L_ventral_l=zeros(iter,M);
L_dorsal_d=zeros(iter,M);
L_ventral_d=zeros(iter,M);
Lmin=zeros(iter,M);
L_dorsal_0M=zeros(iter,M);
L_ventral_0M=zeros(iter,M);

%rest lengths
L0l=zeros(iter,M);
L0d=zeros(iter,M);

%initialising vectors to hold passive forces 
f_dorsal_l=zeros(iter,M);
f_ventral_l=zeros(iter,M);
f_dorsal_d=zeros(iter,M);
f_ventral_d=zeros(iter,M);

%initialising vectors to hold active forces
Fmax=zeros(iter,M);
f_dorsal_M=zeros(iter,M);
f_ventral_M=zeros(iter,M);

%params for active muscle forces
k_dorsal_M=zeros(iter,M);
k_ventral_M=zeros(iter,M);
beta_ventral_M=zeros(iter,M);
beta_dorsal_M=zeros(iter,M);

%---------Accessing the matrix of functions---------------------
I_dorsal=Idorsal(iter,R1,L_dorsal_l,L0l,M,P);
I_ventral=Iventral(iter,R1,L_dorsal_l,L0l,M,P);
S=S(iter,M,I_dorsal,I_ventral);
v_dorsal_d=vdorsal_d(h,iter,M,L_dorsal_d);
v_dorsal_l=vdorsal_l(h,iter,M,L_dorsal_l);
v_ventral_l=v_ventral_l(h,iter,M,L_ventral_l);
v_ventral_d=v_ventral_d(h,iter,M,L_ventral_d);

%------------Neural Activation term
Av=activation_ventral(h,tau,iter,M);
Ad=activation_dorsal(h,tau,iter,M);


%%--------------------Check the logic of how t & i change%%
for t=1:iter
for i=1:P/2
%%------------------passive body forces-----------------------------------
    %defining varying length b/w pi & pi+1 pts

    x_pd(t,i)=x_CoM(t,i)-spacing(t,i)/2;
    x_pv(t,i)=x_CoM(t,i)+spacing(t,i)/2;

    y_pd(t,i)=y_CoM(t,i)+R1(t,i);
    y_pv(t,i)=y_CoM(t,i)-R1(t,i);
    
    if i>1
    L_dorsal_l(t,i-1)=sqrt((x_pd(t,i-1)-x_pd(t,i))^2+(y_pd(t,i-1)-y_pd(t,i))^2);
    L_ventral_l(t,i-1)=sqrt((x_pv(t,i-1)-x_pv(t,i))^2+(y_pv(t,i-1)-y_pv(t,i))^2);
    
    L0l(t,i-1)=sqrt(L_seg^2+(R1(t,i-1)-R1(t,i))^2);
    L0d(t,i-1)=sqrt(L_seg^2+(R1(t,i-1)+R1(t,i))^2);
    
    %----------------------------------------------------------------------
    %%-------------lateral
    %lateral dorsal forces
    if L_dorsal_l(t,i-1)<L0l(t,i-1)
        %to define v_dorsal_l as a function
        f_dorsal_l(t,i-1)=kl*(L0l(t,i-1)-L_dorsal_l(t,i-1))+beta_l*v_dorsal_l(t,i-1);
    else
        f_dorsal_l(t,i-1)=kl*(L0l(t,i-1)-L_dorsal_l(t,i-1))+(2*(L0l(t,i-1)-L_dorsal_l(t,i-1))^4)+beta_l*v_dorsal_l(t,i-1); 
    end

    %lateral ventral forces
    if L_ventral_l(t,i-1)<L0l(t,i-1)
        %to define v_ventral_l as a function
        f_ventral_l(t,i-1)=kl*(L0l(t,i-1)-L_ventral_l(t,i-1))+beta_l*v_ventral_l(t,i-1);
    else
        f_ventral_l(t,i-1)=kl*(L0l(t,i-1)-L_ventral_l(t,i-1)) +(2*(L0l(t,i-1)-L_ventral_l(t,i-1))^4)+beta_l*v_ventral_l(t,i-1); 
    end
    
    %%------------diagonal
    %diagonal dorsal forces
    L_dorsal_d(t,i-1)=sqrt((x_pd(t,i-1)-x_pv(t,i))^2+(y_pd(t,i-1)-y_pv(t,i))^2);
    
    f_dorsal_d(t,i-1)=kd*(L0d(t,i-1)-L_dorsal_d(t,i-1))+beta_d*v_dorsal_d(t,i-1);
    
    %diagonal ventral forces
    L_ventral_d(t,i-1)=sqrt((x_pv(t,i-1)-xpd(t,i))^2+(y_pv(t,i-1)-y_pd(t,i))^2);
    
    f_ventral_d(t,i-1)=kd*(L0d(t,i-1)-L_ventral_d(t,i-1))+beta_d*v_ventral_d(t,i-1);

%%-------------------Active muscle forces-----------
    if i>2
        Fmax(t,i-1)=0.70-0.42*((i-1)-1)/M;
    end
    
    Lmin(t,i-1)=L0l(t,i-1)*(1-del_m*(R1(i-1)+R1(i)/2*R));
    %sigma is the neural activation fn;input to be given->2nd part of code
    
    %dorsal
    k_dorsal_M(t,i-1)=k0M*Fmax(t,i-1)*sigma(Ad(t,i-1));
    L_dorsal_0M(t,i-1)=L0l(t,i-1)-Fmax(t,i-1)*sigma(Ad(t,i-1))*(L0l(t,i-1)-Lmin(t,i-1));
    beta_dorsal_M(t,i-1)=beta_0M*Fmax(t,i-1)*sigma(Ad(t,i-1));
    
    f_dorsal_M(t,i-1)=k_dorsal_M(t,i-1)*(L_dorsal_0M(t,i-1)-L_dorsal_l(t,i-1))+beta_dorsal_M(t,i-1)*v_dorsal_l(t,i-1);
    
    %ventral
    
    k_ventral_M(t,i-1)=k0M*Fmax(t,i-1)*sigma(Av(t,i-1));
    L_ventral_0M(t,i-1)=L0l(t,i-1)-Fmax(t,i-1)*sigma(Av(t,i-1))*(L0l(t,i-1)-Lmin(t,i-1));
    beta_ventral_M(t,i-1)=beta_0M*Fmax(t,i-1)*sigma(Av(t,i-1));
   
    f_ventral_M(t,i-1)=k_ventral_M(t,i-1)*(L_ventral_0M(t,i-1)-L_ventral_l(t,i-1))+beta_ventral_M(t,i-1)*v_ventral_l(t,i-1);
    
    end
    %--------------------Applying forces at each pt pi-------------------------
    if i==1
        
        
        f_dorsal(t,1)=f_dorsal_l(t,1)+f_dorsal_d(t,1)+f_dorsal_M(t,1);
    elseif i==P/2
        
        f_ventral(t,49)=f_ventral_l(t,48)+f_ventral_d(t,48)+f_ventral_M(t,48);
    else
        f_dorsal(t,i)=(f_dorsal_l(t,i-1)-f_dorsal_l(t,i))+(f_dorsal_d(t,i-1)-f_dorsal_d(t,i))+(f_dorsal_M(t,i-1)-f_dorsal_M(t,i));
        f_ventral(t,i)=(f_dorsal_l(t,i-1)-f_dorsal_l(t,i))+(f_ventral_d(t,i-1)-f_ventral_d(t,i))+(f_ventral_M(t,i-1)-f_ventral_M(t,i));
    
        
    %Converting the dorsal and ventral components of the force vector 
    %to perp and parallel
        
        f_dorsal_para(t,i)=f_dorsal(t,i)*cos(90-phi(t,i))-f_dorsal(t,i)*sin(90-phi(t,i));
        f_dorsal_perp(t,i)=f_dorsal(t,i)*sin(90-phi(t,i))+f_dorsal(t,i)*cos(90-phi(t,i));
        
        f_ventral_para(t,i)=f_ventral(t,i)*cos(90-phi(t,i))-f_ventral(t,i)*sin(90-phi(t,i));
        f_ventral_perp(t,i)=f_ventral(t,i)*sin(90-phi(t,i))+f_ventral(t,i)*cos(90-phi(t,i));
        
        f_para_even(t,i)=(f_dorsal_para(t,i)+f_ventral_para(t,i))/2;
        f_para_odd(t,i)=(f_dorsal_para(t,i)-f_ventral_para(t,i))/2;
        
        
    end
end
end

