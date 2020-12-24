format("e",16)
//------------------------------ Richardson ----------------------------------//

//créer une matrice tridiagonal de taille n
function [res]=Create_Mdiag(a,b,c,n)
     B=diag(ones(n-1,1)*a,1) +  diag(ones(n,1)*b) +diag(ones(n-1,1)*c,-1)
    res=B
endfunction


function [x,relres,resvec,it]= richardson(A,b,tol,max_it,x0,alpha)
    it=0
    x=0
    res=b-A*x0
    relres=norm(res)/norm(b)
    
    
    while( (relres > tol) & (it< max_it) )
        it=it+1
       
        x= x0+ alpha * res //daxpy
        res=b-A*x0      //gaxpy
        relres=norm(res)/norm(b)
        resvec(it)=relres
        x0=x
        
    end
endfunction

function [res]= richardson(A,b,epsilon,n,alpha)
    D=diag(diag(A))
    x=zeros(n,1)
    r=[]
    cpt=1
    
    while (norm(b-A*x) > epsilon)
        x=x+alpha*(b-A*x)
        r(cpt)=norm(b-A*x) //residu
        cpt=cpt+1
    end
    
    res=r
endfunction


//-------------------------------|| TEST ||-----------------------------------//
n=50
tol=10^(-8)

A=Create_Mdiag(-1,2,-1,n)


b=rand(n,1)

max_it=100000
x_0=zeros(n,1)


for i= 2:1:8
    alpha=1/i
    richard=log10(richardson(A,b,tol,n,alpha))
    
    s1=size(richard,1) // Jacobi
    r_richard=1:1:s1
    plot2d(r_richard,richard,style = i)

end


xlabel("nb iteration");
ylabel("log10(erreur residuelle)");
hl=legend(['i=2';'i=3';'i=4';'i=5';'i=6';'i=7';])
