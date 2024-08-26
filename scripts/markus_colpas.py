import numpy as np

N=8

if (N==6):
    data=np.loadtxt('Hk_bare_phonons.dat')
if (N==8):
    data=np.loadtxt('Hk_bare_magnons_and_phonons.dat')
    #data=np.loadtxt('Hk_full.dat')
    

for i in range (1,len(data[:,0])):
    Hk=data[i,3::2].reshape(N,N)
    Hk=np.round(Hk,10)
    if (N==6):
        g=np.diag([1,1,1,-1,-1,-1])
    elif (N==8):
        g=np.diag([1,1,1,1,-1,-1,-1,-1])
    
    #COLPAS METHOD
    #Check if Hk is positive definite    
    eigenvalues, eigenvectors =np.linalg.eig(Hk)
    for n in range(0,len(eigenvalues)):
        if np.real(eigenvalues[n]) < 0.0:
            print('Hk is not positive definite')
    if(np.max(np.abs(Hk-Hk.conj().T)) >1e-10):
            print('Hk is not Hermitian')
    #Cholesky decomposition
    K=np.linalg.cholesky(Hk)
    K=K.conj().T
    if(np.max(np.abs(K.conj().T@K-Hk)) > 1e-10):
        print('Cholesky decomposition has failed')
    #Solve for L matrix
    eigenvalues, eigenvectors =np.linalg.eig(K@g@K.conj().T)
    if(np.max(np.abs(np.linalg.inv(eigenvectors)@K@g@K.conj().T@eigenvectors-np.diag(eigenvalues)))>1e-10):
        print('Solving EV problem has failed')
    #sort by descending eigenvalues
    ib = np.argsort(-eigenvalues)
    L=np.diag(eigenvalues[ib])
    U=eigenvectors[:,ib]
    U=np.linalg.qr(U)[0]
    if(np.max(np.abs(np.linalg.inv(U)@K@g@K.conj().T@U-L))>1e-10):
        print('Reordering U and L has failed')    
    #check if U is unitary
    if ( (np.max(np.abs(U.conj().T-np.linalg.inv(U))))> 1e-10):
        print('U is not unitary!')
    #get diagonalized matrix of original diagonalization problem
    E=np.linalg.inv(g)@L
    #Check if E is diagonal
    if( (np.max(np.abs(E[~np.eye(E.shape[0],dtype=bool)].reshape(E.shape[0],-1)))) > 1e-10):
        print('E is not diagonal')
    
    #Check if eigenvalues conincide with what can be obtained from the White method
    eigenvaluesWhite, eigenvectorsWhite =np.linalg.eig(g@Hk)
    ibWhite = np.argsort(-eigenvaluesWhite)
    eigenvaluesWhite[:]=eigenvaluesWhite[ibWhite]
    ibnew = np.argsort(-eigenvalues)
    eigenvalues = eigenvalues[ibnew]
    if ( np.max(np.abs(eigenvaluesWhite[0:3]-eigenvalues[0:3]))>1e-10):
        print('EV of Colpas and Whites method do not conincide')

    #get Transformation matrices
    Esqrt=np.zeros(N*N).reshape(N,N)
    for n in range (0,N):    
        Esqrt[n,n]=np.sqrt(E[n,n])
    Q=np.linalg.inv(K)@U@Esqrt
    
    #Check if bosonic commutation relation is conserved
    if(np.max(np.abs(Q.conj().T@g@Q.conj()-g)) >1e-10 ):
        print('Bosonic commuation relation is not conserved for line', i)
    #Check if Q^\dagger H Q= E is fulfilled
    if (np.max(np.abs(Q.conj().T@Hk@Q-E)) > 1e-10):
        print('Diagonalization of Hamiltonian failed for line', i)
        print('\t',*data[i,0:3],E[0,0],E[1,1],E[2,2])
        print('\t Error',np.max(np.abs(Q.conj().T@Hk@Q-E)))
            
    if (np.abs(data[i,0]) + np.abs(data[i,1]) < 0.001):    
        print(i,*data[i,0:3]/2./np.pi,E[0,0],E[1,1],E[2,2],E[3,3], Hk[:,:])