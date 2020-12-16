format("e",16)
//------------------- Facto Lu pour matrice tridiagonal -----------------------//

/*
Apres l'etude theorique nous povons calculer les coef des diagonales a l'aide de
suite.
*/

//créer une matrice tridiagonal de taille n
function [res]=Create_Mdiag(n)
    //B=diag(ones(n,1),1) +  diag(ones(n,1),1) + diag(ones(n,1)*-1,-1)
     B=diag(ones(n-1,1)*-1,1) +  diag(ones(n,1)*2) +diag(ones(n-1,1)*-1,-1)
    res=B
endfunction

//calcul facto LU d'une matrice tridiagonal
function [res]=facto_LU(A)
    result=[]
    n=size(A,1)
    c=diag(A,1)  // recuperation surdiagonal
    a=diag(A)    // recuperation diag
    b=diag(A,-1) // recuperation sous diag 
    
    m=[]
    m(1)=b(1)/a(1)
    
    for i= 1:1:n-2 // jusqu'a n+1
        m(i+1)= b(i+1) / ( a(i+1)-( m(i)*c(i) ))
    end
    
    d=[]
    d(1)=a(1)
    for i=2:1:n
        d(i)=a(i) - (m(i-1) * c(i-1) )
    end
    
    A= diag(c,1)+ diag(d) +diag(m,-1)
    res=A
endfunction

function[matrix_L]=generate_L(A)
    n=size(A,1)
    matrix_L=diag(ones(n,1)) + diag(diag(A,-1),-1)
endfunction


function[matrix_U]=generate_U(A)
    n=size(A,1)
    matrix_U=diag(diag(A))+diag(diag(A,1),1)
endfunction


//----------------------------------test------------------------------------- //

A=Create_Mdiag(1000)

AB=facto_LU(A)
matrix_L=generate_L(AB)
matrix_U=generate_U(AB)

relres=norm(matrix_L*matrix_U)/norm(A)
disp(relres)
//disp("difference A-LU",A- (matrix_L*matrix_U))




