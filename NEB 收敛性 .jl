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












