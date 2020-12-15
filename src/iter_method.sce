//format("e",16)
//--------------------------- methode iterative ------------------------------//

//créer une matrice tridiagonal de taille n
function [res]=Create_Mdiag(a,b,c,n)compl
     B=diag(ones(n-1,1)*a,1) +  diag(ones(n,1)*b) +diag(ones(n-1,1)*c,-1)
    res=B
endfunction
// methode de Jacobi //

function [res]= jacobi(A,b,epsilon,n)
    D=diag(diag(A))
    x=zeros(n,1)
    r=[]
    cpt=1
    
    while (norm(b-A*x) > epsilon)
        x=x+inv(D)*(b-A*x)
        r(cpt)=norm(b-A*x) //residu
        cpt=cpt+1
    end
    
    res=r
endfunction

// Methode de Gauss Sidel //

//norme residu / norme b 

function [res]= Gauss_Seidel(A,b,epsilon,n)
    D=diag(diag(A))
    E=tril(A)
    r=[]
    x=zeros(n,1)
    cpt=1
    
    while (norm(b-A*x) > epsilon)
        x=x+inv(E)*(b-A*x)
        r(cpt)=norm(b-A*x) //residu
        cpt=cpt+1
        
    end
    
    res=r
endfunction

//-----------------TEST-----------------//
n=50
A=Create_Mdiag(-1,2,-1,n)
//disp(A)
b=rand(n,1)

x_ex=A\b

eps=0.00000001
jacobi=log10(jacobi(A,b,eps,n))
GS=log10(Gauss_Seidel(A,b,eps,n))

//disp("result Jacobi",jacobi)
//disp("result Gauss_sidel",GS)
//disp("result exact",x_ex)


// graphique //

s1=size(GS,1) //Gauss Seidel
r_GS=1:1:s1

s2=size(jacobi,1) // Jacobi
r_jacobi=1:1:s2

plot(r_GS,GS,r_jacobi,jacobi)
hl=legend(['Gauss Seidel';'Jacobi';])
xlabel("nb iteration");
ylabel("log10(erreur residuelle)");
//disp(s1)
