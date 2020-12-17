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
    relres=norm(b-A*x0)/norm(b)
    
    
    while( (relres>tol) & (it< max_it) )
        it=it+1
        x=x0+ alpha*res //daxpy
        res=b-A*x0      //gaxpy
        relres=norm(b-A*x0)/norm(b)
        resvec(it)=relres
        x0=x
    end
    
    
endfunction

//-------------------------------|| TEST ||-----------------------------------//
n=5
tol=0,00000001
A=Create_Mdiag(-1,2-1,n)
b=rand(n,1)
max_it=10000
x_0=ones(n,1)
alpha=1/2
[x,relres,resvec,it]=richardson(A,b,tol,max_it,x_0,alpha)

for i= 1:1:4
    alpha=A/i
    [x,relres,resvec,it]=richardson(A,b,tol,max_it,x_0,alpha)
    plot(1:1:it,resvec)

end


