using Plots
using LaTeXStrings
using Calculus
using LinearAlgebra
tmax=0.17;
t0=0.01;#Time step Δt
c=0.000001;#NEB's parameter c

h=0.0125#Space step
N=3000;# Maximum steps of iteration
M=Int32(1/h);#Numbers of images
i=1; 


E(x)= (1-x[1]^2-x[2]^2)^2+((x[2]^2)./(x[1]^2+x[2]^2))
#Energy function E(x)
GE1(x)=[-4*x[1]*(1-x[1]^2-x[2]^2)-2*x[1]*(x[2]^2)./(x[1]^2+x[2]^2),-4*x[2]*(1-x[1]^2-x[2]^2)+2*x[2]*(x[1]^2)./(x[1]^2+x[2]^2)]
#Gradient ∇E(x)

A1=[-1,0];
A2=[1,0];
#two minima

f1(t)=A1+t*(A2-A1)+[0,t-t.^2]
#initial value
fs(x)=[-cos(pi*x),sin(pi*x)];
# the MEP of E(x)


j=1;
s=zeros(2*(M+1));#Save the discrete MEP
for j=1:1:M+1
 s[2*j-1:2*j]=fs(h*(j-1)); 
end
m=1
tol=10^(-7);# Tolerance



P=zeros(2*(M-1),2*(M-1));
P1=zeros(2*(M-1),2*(M-1));

for j=2:1:M-1
    P1[2*j-1:2*j,2*j-1:2*j]=-I(2);
    P1[2*j-1:2*j,2*(j-1)-1:2*(j-1)]=I(2);
end
P1[2*1-1:2*1,2*1-1:2*1]=-I(2);
P=-inv(P1)
#Preconditioner P

e1=ones(1,N+1);#∞norm of F_{h}(ϕ)
e2=ones(1,N+1);#∞norm of the normal gradient force -∇E(ϕ)^⟂
e3=ones(1,N+1);#∞norm of tangent force
e4=ones(1,N+1);#∞norm of ϕ minus MEP
 
i=1;
U=zeros(2*(M+1),N+1);#Save the evolution of the strings 
T=zeros(2*(M-1),N+1);#Save the tangent force 
PT=zeros(2*(M-1),N+1);#Save the tangent force after precondition
GEF=zeros(2*(M-1),N+1);#Save the normal gradient force -∇E(x)^⟂
K=zeros(2*(M+1),N+1);#Save the unit tangent vector φ dot 
K2=zeros(2*(M+1),N+1);#Save |ϕ_{i+1}-ϕ_{i}|  
    
for j=1:1:M+1
U[2*j-1:2*j,1]=f1(h*(j-1))
end 
#input the initial value    

while((i<N+1)&&(e1[1,i]>tol))
    
    
   for j = 2:1:M
            if E(U[2*(j-1)-1:2*(j-1),i])<E(U[2*j-1:2*j,i])<E(U[2*(j+1)-1:2*(j+1),i])
                K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*j-1:2*j,i])./norm((U[2*(j+1)-1:2*(j+1),i] - U[2*j-1:2*j,i]))
            else
                if E(U[2*(j-1)-1:2*(j-1),i])>E(U[2*j-1:2*j,i])>E(U[2*(j+1)-1:2*(j+1),i])
                    K[2*j-1:2*j,i] = (U[2*j-1:2*j,i] - U[2*(j-1)-1:2*(j-1),i])./norm((U[2*j-1:2*j,i] - U[2*(j-1)-1:2*(j-1),i]))
                else
                    K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*(j-1)-1:2*(j-1),i])./norm((U[2*(j+1)-1:2*(j+1),i] - U[2*(j-1)-1:2*(j-1),i]))     
            end
        end
    end   
     K[2*1-1:2*1,i]=(U[2*2-1:2*2,i]-U[2*(2-1)-1:2*(2-1),i])./norm(U[2*2-1:2*2,i]-U[2*(2-1)-1:2*(2-1),i]);
    K[2*(M+1)-1:2*(M+1),i]=(U[2*(M+1)-1:2*(M+1),i]-U[2*(M)-1:2*(M),i])./norm(U[2*(M+1)-1:2*(M+1),i]-U[2*(M)-1:2*(M),i]);
    #Calculate the unit tangent vector φ dot
    
    for j = 2:1:M
    K2[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*(j)-1:2*(j),i])
    end   
     K2[2*1-1:2*1,i]=U[2*2-1:2*2,i]-U[2*(2-1)-1:2*(2-1),i];
    K2[2*(M+1)-1:2*(M+1),i]=U[2*(M+1)-1:2*(M+1),i]-U[2*(M)-1:2*(M),i];
    #Calculate |ϕ_{i+1}-ϕ_{i}|  
    
      
    
    
        
    for j=2:1:M
    T[2*(j-1)-1:2*(j-1),i]=-c*(norm(K2[2*j-1:2*j,i])-norm(K2[2*(j+1)-1:2*(j+1),i]))*(K[2*j-1:2*j,i]);
    end
    #Calculate the tangent force
    PT[:,i]=P*T[:,i];  
    #Precondition
        
    for j=2:1:M
    GEF[2*(j-1)-1:2*(j-1),i]=-GE1(U[2*j-1:2*j,i])+(dot(GE1(U[2*j-1:2*j,i]),K[2*j-1:2*j,i]))*(K[2*j-1:2*j,i]);
    end
    #Calculate the normal gradient force -∇E(x)^⟂
    
    
    U[2*2-1:2*M,i+1]=t0*(GEF[:,i]+T[:,i])+U[2*2-1:2*M,i]
    U[2(M+1)-1:2(M+1),i+1]=A2;
    U[2*1-1:2,i+1]=A1;  
    #Calculate ϕ^{n} 
    
    
    
    e1[1,i+1]=maximum(abs.(GEF[:,i]+T[:,i]));#Calculate ∞norm of F_{h}(ϕ)
    e2[1,i+1]=maximum(abs.(GEF[:,i]));#Calculate ∞norm of the normal gradient force -∇E(ϕ)^⟂
    e3[1,i+1]=maximum(abs.(T[:,i]));#Calculate ∞norm of tangent force
    e4[1,i+1]=maximum(abs.(U[:,i]-s));#Calculate ∞norm of ϕ minus MEP
    i=i+1;  
end


e2[1,1:i]

e5=zeros(M-1)
for j=2:1:M
e5[j-1]=U[2j-1,i]^2+U[2j,i]^2-1
end
e6=maximum(abs.(e5))

i


U[:,i]

P

using Plots
using LaTeXStrings
using Calculus
using LinearAlgebra
tmax=0.2;
t0=0.05
N2=50;
TS=zeros(1,N2);#time step
IS=zeros(1,N2);
c=1;

h=0.025#Space step
N=1000;# Maximum steps of iteration
M=Int32(1/h);
i=1;


E(x)= (1-x[1]^2-x[2]^2)^2+((x[2]^2)./(x[1]^2+x[2]^2))
GE1(x)=[-4*x[1]*(1-x[1]^2-x[2]^2)-2*x[1]*(x[2]^2)./(x[1]^2+x[2]^2),-4*x[2]*(1-x[1]^2-x[2]^2)+2*x[2]*(x[1]^2)./(x[1]^2+x[2]^2)]
A1=[-1,0];
A2=[1,0];
f1(t)=[-1+2t,t*(1-t)]
fs(x)=[-cos(pi*x),sin(pi*x)];# the function expected to convergent
j=1;
s=zeros(2*(M+1));
Ks=zeros(M+1,2);
e=zeros(1,N+1);
for j=1:1:M+1
 s[2j-1:2*j]=fs(h*(j-1)); 
end
m=1
tol=10^(-7);# tolerance
P=zeros(2*(M-1),2*(M-1));
P1=zeros(2*(M-1),2*(M-1));

for j=2:1:M-1
    P1[2*j-1:2*j,2*j-1:2*j]=-I(2);
    P1[2*j-1:2*j,2*(j-1)-1:2*(j-1)]=I(2);
end
P1[2*1-1:2*1,2*1-1:2*1]=-I(2);

P=-inv(P1)



e1=ones(1,N+1);
e2=ones(1,N+1);
e3=ones(1,N+1);
e4=ones(1,N+1);
 
i=1;
U=zeros(2*(M+1),N+1);
T=zeros(2*(M-1),N+1);
PT=zeros(2*(M-1),N+1);
GEF=zeros(2*(M-1),N+1);
K=zeros(2*(M+1),N+1);
 K2=zeros(2*(M+1),N+1);  
l=zeros(M+1,N+1);
    
for j=1:1:M+1
U[2*j-1:2*j,1]=f1(h*(j-1))
end 
    

while((i<N+1)&&(e1[1,i]>tol))
    
    
     for j = 2:1:M
            if E(U[2*(j-1)-1:2*(j-1),i])<E(U[2*j-1:2*j,i])<E(U[2*(j+1)-1:2*(j+1),i])
                K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*j-1:2*j,i])
            else
                if E(U[2*(j-1)-1:2*(j-1),i])>E(U[2*j-1:2*j,i])>E(U[2*(j+1)-1:2*(j+1),i])
                    K[2*j-1:2*j,i] = (U[2*j-1:2*j,i] - U[2*(j-1)-1:2*(j-1),i])
                else
                    K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*(j-1)-1:2*(j-1),i])
                
            end
        end
    end
     
     K[2*1-1:2*1,i]=U[2*2-1:2*2,i]-U[2*(2-1)-1:2*(2-1),i];
    K[2*(M+1)-1:2*(M+1),i]=U[2*(M+1)-1:2*(M+1),i]-U[2*(M)-1:2*(M),i];
    
    
    
    for j=2:1:M+1
        l[j,i]=l[j-1,i]+norm(K[2*j-1:2*j,i]);
    end
    
        
    for j=2:1:M
    T[2*(j-1)-1:2*(j-1),i]=-(norm(K[2*j-1:2*j,i])./h-l[M+1,i])*(K[2*j-1:2*j,i]);
    end
    PT[:,i]=P*T[:,i];  
        
    for j=2:1:M
    GEF[2*(j-1)-1:2*(j-1),i]=-GE1(U[2*j-1:2*j,i])+((dot(GE1(U[2*j-1:2*j,i]),K[2*j-1:2*j,i])./(norm(K[2*j-1:2*j,i],2)^2))*(K[2*j-1:2*j,i]));
    end
    
    U[2*2-1:2*M,i+1]=t0*(GEF[:,i]+PT[:,i])+U[2*2-1:2*M,i]
    U[2(M+1)-1:2(M+1),i+1]=A2;
    U[2*1-1:2,i+1]=A1;  
    
    
    
    
    e1[1,i+1]=maximum(abs.(GEF[:,i]+T[:,i]));
    e2[1,i+1]=maximum(abs.(GEF[:,i]));
    e3[1,i+1]=maximum(abs.(T[:,i]));
    e4[1,i+1]=maximum(abs.(U[:,i]-s));
    i=i+1; 
        
    end




e1[1,i]

 for j = 2:1:M
            if E(U[2*(j-1)-1:2*(j-1),i])<E(U[2*j-1:2*j,i])<E(U[2*(j+1)-1:2*(j+1),i])
                K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*j-1:2*j,i])
            else
                if E(U[2*(j-1)-1:2*(j-1),i])>E(U[2*j-1:2*j,i])>E(U[2*(j+1)-1:2*(j+1),i])
                    K[2*j-1:2*j,i] = (U[2*j-1:2*j,i] - U[2*(j-1)-1:2*(j-1),i])
                else
                    K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*(j-1)-1:2*(j-1),i])
                end
            end
    

e5=zeros(M-1)
for j=2:1:M
e5[j-1]=U[2j-1,i]^2+U[2j,i]^2-1
end
e6=maximum(abs.(e5))

e2[1,i-1]

l[M+1,1]

using Plots
using LaTeXStrings
using Calculus
using LinearAlgebra
tmax=0.7;
t0=0.02;
N2=50;
TS=zeros(1,N2);#time step
IS=zeros(1,N2);
c=1;

h=0.01#Space step
N=1000;# Maximum steps of iteration
M=Int32(1/h);
i=1;


E(x)= (1-x[1]^2-x[2]^2)^2+((x[2]^2)./(x[1]^2+x[2]^2))
GE1(x)=[-4*x[1]*(1-x[1]^2-x[2]^2)-2*x[1]*(x[2]^2)./(x[1]^2+x[2]^2),-4*x[2]*(1-x[1]^2-x[2]^2)+2*x[2]*(x[1]^2)./(x[1]^2+x[2]^2)]
A1=[-1,0];
A2=[1,0];
f1(t)=[-1+2t,t*(1-t)]
fs(x)=[-cos(pi*x),sin(pi*x)];# the function expected to convergent
j=1;
s=zeros(2*(M+1));
Ks=zeros(M+1,2);
e=zeros(1,N+1);
for j=1:1:M+1
 s[2j-1:2*j]=fs(h*(j-1)); 
end
m=1
tol=10^(-7);# tolerance
#P=I(2*(M-1));
P=zeros(2*(M-1),2*(M-1));
P1=zeros(2*(M-1),2*(M-1));

for j=2:1:M-1
    P1[2*j-1:2*j,2*j-1:2*j]=-I(2);
    P1[2*j-1:2*j,2*(j-1)-1:2*(j-1)]=I(2);
end
P1[2*1-1:2*1,2*1-1:2*1]=-I(2);

P=-inv(P1)



e1=ones(1,N+1);
e2=ones(1,N+1);
e3=ones(1,N+1);
e4=ones(1,N+1);
 
i=1;
U=zeros(2*(M+1),N+1);
T=zeros(2*(M-1),N+1);
PT=zeros(2*(M-1),N+1);
GEF=zeros(2*(M-1),N+1);
K=zeros(2*(M+1),N+1);
 K2=zeros(2*(M+1),N+1);  
l=zeros(M+1,N+1);
    

U[:,1]=s2
 
    

while((i<N+1)&&(e1[1,i]>tol))
    
    
    for j = 2:1:M
            if E(U[2*(j-1)-1:2*(j-1),i])<E(U[2*j-1:2*j,i])<E(U[2*(j+1)-1:2*(j+1),i])
                K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*j-1:2*j,i])
            else
                if E(U[2*(j-1)-1:2*(j-1),i])>E(U[2*j-1:2*j,i])>E(U[2*(j+1)-1:2*(j+1),i])
                    K[2*j-1:2*j,i] = (U[2*j-1:2*j,i] - U[2*(j-1)-1:2*(j-1),i])
                else
                    K[2*j-1:2*j,i] = (U[2*(j+1)-1:2*(j+1),i] - U[2*(j-1)-1:2*(j-1),i])
                
            end
        end
    end
     
     K[2*1-1:2*1,i]=U[2*2-1:2*2,i]-U[2*(2-1)-1:2*(2-1),i];
    K[2*(M+1)-1:2*(M+1),i]=U[2*(M+1)-1:2*(M+1),i]-U[2*(M)-1:2*(M),i];
    
      
    for j=2:1:M+1
        l[j,i]=l[j-1,i]+norm(K[2*j-1:2*j,i]);
    end
    
        
    for j=2:1:M
    T[2*(j-1)-1:2*(j-1),i]=-c*(norm(K[2*j-1:2*j,i])./h-l[M+1,i])*(K[2*j-1:2*j,i]);
    end
    PT[:,i]=P*T[:,i];  
        
    for j=2:1:M
    GEF[2*(j-1)-1:2*(j-1),i]=-GE1(U[2*j-1:2*j,i])+((dot(GE1(U[2*j-1:2*j,i]),K[2*j-1:2*j,i])./(norm(K[2*j-1:2*j,i],2)^2))*(K[2*j-1:2*j,i]));
    end
    
    U[2*2-1:2*M,i+1]=t0*(T[:,i])+U[2*2-1:2*M,i]
    U[2(M+1)-1:2(M+1),i+1]=A2;
    U[2*1-1:2,i+1]=A1;  
    
    
    
    
    e1[1,i+1]=maximum(abs.(GEF[:,i]+T[:,i]));
    e2[1,i+1]=maximum(abs.(GEF[:,i]));
    e3[1,i+1]=maximum(abs.(T[:,i]));
    e4[1,i+1]=maximum(abs.(U[:,i]-s));
    i=i+1;  
end




e2


