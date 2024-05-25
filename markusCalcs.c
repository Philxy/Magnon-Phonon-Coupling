#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>    /* Standard Library of Complex Numbers */
#define LAT_CONST (1.0)
#define PI (3.1415926)
#define mRy_to_eV (13.605693122994e-3)
#define HBAR (6.582119569e-16) //eV*s
#define KB (8.617e-5) //eV per kelvin
#define J_sd  (0.75/27.211) // in Hartree; extracted by comparison of sp- bands with up/down at the fermi level for 20^3 kpoints



#define GAMMA (1.76085963023e11) // per second per tesla
#define MUS (2.2*5.788381806e-5) //eV per tesla

#define S_QN (1.0/2.2)

#ifndef SWIDTH
	#define SWIDTH (0.1*0.01*27.211) // Gaussian smearing in eV (0.01 Hartree is the standard smearing in elk)
#endif

#define N_q_per_line 20

#define N_QPTSNR (N_q_per_line*N_q_per_line*N_q_per_line)
#if N_q_per_line==10
	#define N_QPTS 47
#elif N_q_per_line==16
	#define N_QPTS 145
#elif N_q_per_line==18
	#define N_QPTS 195
#elif N_q_per_line==19
	#define N_QPTS 220
#elif N_q_per_line==20
	#define N_QPTS 256
#elif N_q_per_line==24
	#define N_QPTS 413
#elif N_q_per_line==30
	#define N_QPTS 752
#endif
#define N_SYM 48


#define DZ (6.97e-6) //in eV ;  I. Razdolski et al Nat comm 2017
//shell_index; data from Mraysov et al., (1996); Units in mRy; listed in PHYSICAL REVIEW B 82, 144304 (2010)
//1 1.24
//2 0.646
//3 0.007
//4 -0.108
//5 -0.071
//6 0.035
//7 0.002
//8 0.014

double dot_product( double a[3], double b[3]){
	double sum =0.0;
	sum+=a[0]*b[0];
	sum+=a[1]*b[1];
	sum+=a[2]*b[2];
	return sum;
}


int cross_product(double *a, double *b, double *c){
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = -a[0]*b[2]+a[2]*b[0];
	c[2] = a[0]*b[1]-a[1]*b[0];
	return 0;
}



int initialize_lattice (double *a1,double *a2,double *a3,double *b1,double *b2,double *b3){
	a1[0]=LAT_CONST;	a1[1]=LAT_CONST;	a1[2]=-LAT_CONST;
	a2[0]=LAT_CONST;	a2[1]=-LAT_CONST;	a2[2]=LAT_CONST;
	a3[0]=-LAT_CONST;	a3[1]=LAT_CONST;	a3[2]=LAT_CONST;

	cross_product(a2,a3,b1);
	cross_product(a3,a1,b2);
	cross_product(a1,a2,b3);
	
	double transform_factor;
	transform_factor=2*PI/dot_product(a1,b1);
	int n;
	for (n=0;n<3;n++){
		b1[n]*=transform_factor;
		b2[n]*=transform_factor;
		b3[n]*=transform_factor;
	}

	return 0;
}





double calculate_frequency(double a1[3],double a2[3],double a3[3],double k_pos[3]){
	double d;
	double proj;
	double frequency;
	double pos[3],J_list[9];
	int i,j,k,n,list_index;
	J_list[0]=0.0;
	J_list[1]=1.24;
	J_list[2]=0.646;
	J_list[3]=0.007;
	J_list[4]=-0.108;
	J_list[5]=-0.071;
	J_list[6]=0.035;
	J_list[7]=0.002;
	J_list[8]=0.014;

	frequency=0.0;
	for (int i=-3;i<=3;i++){
	for (int j=-3;j<=3;j++){
	for (int k=-3;k<=3;k++){
		for (int n=0;n<3;n++){ 
			pos[n]=i*a1[n]	+ j*a2[n] +k*a3[n];	
		}
		d=sqrt(dot_product(pos,pos));
		if (d< 0.001){
			d+=1000.0;
		}
		if (d<1.8){ list_index=1;}
		else if (d<2.2) {list_index=2;}
		else if (d<3.0) {list_index=3;}
		else if (d<3.4) {list_index=4;}
		else if (d<3.5) {list_index=5;}
		else if (d<4.1) {list_index=6;}
		else if (d<4.4) {list_index=7;}
		else if (d<4.5) {list_index=8;}
		else {list_index=0;}
		if (d<1000.0){
			proj=dot_product(k_pos,pos);
			frequency+=J_list[list_index]*(1.-cos(proj)); //FACTOR TWO DUE TO DIFFERENT CONVENTIONS IN MY MODEL AND THE MODEL WITHIN WHICH THE VALUES ARE DEFINED
		}
	}
	}
	}
	frequency=frequency*mRy_to_eV;
	return (2.*DZ+frequency)/HBAR*S_QN/(S_QN*S_QN); // 1/S^2 is to get from a spin model with unit vectors to vectors of length S; multiplication with S is due to HP Trafo for spins with length S
}


double d_N_BE_by_dT (double energy,double T){
	double x=energy/(KB*T);
	return exp(x)/((exp(x)-1)*(exp(x)-1))*x/T;
}





int get_equivalent_index(double k_posl[3], double k_pos_index[N_QPTS][3],double b1[3],double b2[3],double b3[3]){
	int equivalent_index=0;
	int A[48][3][3],dummy,n,i,j,k,l,m; // A = symmetry matrices
	
	//READ MATRICES
	FILE *f=fopen("symmetry_matrices.dat","r");
	fscanf(f,"%d",&(dummy));
	if (dummy != N_SYM){
		printf("Wrong number of symmetry matrices!\n");
		return -1;
	}
	for (n=0;n<N_SYM/2;n++){
		fscanf(f,"%d %d %d %d %d %d %d %d %d",&(A[n][0][0]),&(A[n][1][0]),&(A[n][2][0]),&(A[n][0][1]),&(A[n][1][1]),&(A[n][2][1]),&(A[n][0][2]),&(A[n][1][2]),&(A[n][2][2]));
	}
	fclose(f);

	//ADD REST OF MATRICES (THEY ARE JUST THE NEGATIVE OF THE OTHERS)
	for (n=N_SYM/2;n<N_SYM;n++){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				A[n][i][j]=-A[n-N_SYM/2][i][j];
			}
		}	
					

	}
	for (n=0;n<N_SYM;n++){
		//printf("%d %d %d \n",A[n][0][0],A[n][0][1],A[n][0][2]);
		//printf("%d %d %d \n",A[n][1][0],A[n][1][1],A[n][1][2]);
		//printf("%d %d %d \n\n",A[n][2][0],A[n][2][1],A[n][2][2]);	
	}

	

	//SEARCH FOR EQUIVALENT K-point
	double k_comparec[3], k_pos_tf[3],k_pos_tf_shifted[3],diff[3],distance;
	for (n=0;n<N_SYM;n++){
		

		//MULTIPLY WITH MATRIX
		for(i=0;i<3;i++){
			k_pos_tf[i]=k_posl[0]*A[n][i][0] +	k_posl[1]*A[n][i][1] + k_posl[2]*A[n][i][2];
		}
		//COMPARE WITH K-point list
		//MAP BACK TO first BZ
		for(i=-2;i<=2;i++){
		for(j=-2;j<=2;j++){
		for(k=-2;k<=2;k++){

			k_pos_tf_shifted[0]=k_pos_tf[0]+i*1.0;
			k_pos_tf_shifted[1]=k_pos_tf[1]+j*1.0;
			k_pos_tf_shifted[2]=k_pos_tf[2]+k*1.0;
			//CALCULATE DIFFERENCE TO ENTRY IN LIST OF IRREDUCIBLE K-POINTS
			for(l=0;l<N_QPTS;l++){

				k_comparec[0] =k_pos_index[l][0];
				k_comparec[1] =k_pos_index[l][1];
				k_comparec[2] =k_pos_index[l][2];

				distance=0.0;
				for(m=0;m<3;m++){
					distance+=pow(k_pos_tf_shifted[m]-k_comparec[m],2);
				}
				distance = sqrt(distance);
				if (distance<0.001){
					equivalent_index = l;
					return equivalent_index;
				}

			}

		}
		}
		}

	}
	return -1;
}




double delta_distribution(double x){
	return 1.0/(sqrt(PI))*exp(-x*x);
}







int main (){
	double a1[3],a2[3],a3[3];	
	double b1[3],b2[3],b3[3];
	double k_pos_index[N_QPTS][3];
	double k_pos[3];
	double Lambda_k[N_QPTS],Lambda_k_uu[N_QPTS],Lambda_k_dd[N_QPTS];
	double frequency_k[N_QPTS],gamma_k[N_QPTS],alpha_eff_k[N_QPTS];
	double C_k[N_QPTS];
	double C_mag,G_emag_1,G_emag_2;
	int n,dummy,multiplicity_k[N_QPTS],total_nqpt;
	initialize_lattice(a1,a2,a3,b1,b2,b3);


	total_nqpt=0;
	C_mag = 0.0;
	G_emag_1 = 0.0;
	G_emag_2 = 0.0;



	//READ DATA AND CALCULATE k dependent quantities
#if N_q_per_line==10
//	FILE *f=fopen("lifetimes_10_qgrid_30_kgrid_g_ud_g_uu_g_dd.dat","r");
	FILE *f=fopen("lifetimes_10_qgrid_20_kgrid_g_ud_g_uu_g_dd.dat","r");
//	FILE *f=fopen("lifetimes_10_qgrid_10_kgrid_g_ud_g_uu_g_dd.dat","r");
#elif N_q_per_line==16
	FILE *f=fopen("lifetimes_16_qgrid_16_kgrid_g_ud_g_uu_g_dd.dat","r");
#elif N_q_per_line==18
	FILE *f=fopen("lifetimes_18_qgrid_18_kgrid_g_ud_g_uu_g_dd.dat","r");
#elif N_q_per_line==19
	FILE *f=fopen("lifetimes_19_qgrid_19_kgrid_g_ud_g_uu_g_dd.dat","r");
#elif N_q_per_line==20
//	FILE *f=fopen("lifetimes_20_qgrid_20_kgrid_g_ud_g_uu_g_dd.dat","r");
	FILE *f=fopen("lifetimes_20_qgrid_20_kgrid_g_ud_g_uu_g_dd_no_interband.dat","r");
#elif N_q_per_line==24
	FILE *f=fopen("lifetimes_24_qgrid_24_kgrid_g_ud_g_uu_g_dd.dat","r");
#elif N_q_per_line==30
	FILE *f=fopen("lifetimes_30_qgrid_30_kgrid_g_ud_g_uu_g_dd.dat","r");
#endif
	




	
	//printf("#i vql(1) vql(2) vql(3) cart_q_x cart_q_y cart_q_z w(Hz) Lambda_k (1/Ha^2) alpha_eff_k gamma_k C_k\n");
	for(n=0;n<N_QPTS;n++){
		fscanf(f,"%d %lf %lf %lf %lf %lf %lf %d",&(dummy),&(k_pos_index[n][0]),&(k_pos_index[n][1]),&(k_pos_index[n][2]),&(Lambda_k[n]),&(Lambda_k_uu[n]),&(Lambda_k_dd[n]), &(multiplicity_k[n]));

		#if N_q_per_line==18
			multiplicity_k[n]++;
		#else
		if ( multiplicity_k[n] > 1 && multiplicity_k[n]%2 != 0){
			multiplicity_k[n]++;
		}
		#endif
		k_pos[0] = k_pos_index[n][0]*b1[0] + k_pos_index[n][1]*b2[0] + k_pos_index[n][2]*b3[0];
		k_pos[1] = k_pos_index[n][0]*b1[1] + k_pos_index[n][1]*b2[1] + k_pos_index[n][2]*b3[1];
		k_pos[2] = k_pos_index[n][0]*b1[2] + k_pos_index[n][1]*b2[2] + k_pos_index[n][2]*b3[2];
		
		frequency_k[n]=calculate_frequency(a1,a2,a3,k_pos);
		alpha_eff_k[n]=Lambda_k[n]*J_sd*J_sd*4.*PI/S_QN;
		gamma_k[n]=frequency_k[n]*alpha_eff_k[n];
		
		C_k[n]=HBAR*frequency_k[n]*d_N_BE_by_dT(HBAR*frequency_k[n],300.);

		
		total_nqpt+=multiplicity_k[n];
		
		C_mag+=C_k[n]*multiplicity_k[n]*1.0;
		G_emag_1+=C_k[n]*gamma_k[n]*multiplicity_k[n]*1.0;

		//print
		printf("%d %e %e %e %e %e %d \n",n, k_pos_index[n][0],k_pos_index[n][1],k_pos_index[n][2],frequency_k[n],gamma_k[n],multiplicity_k[n]);
		
	}
	fclose(f);
	if (total_nqpt != N_QPTSNR){	
		printf("Wrong number of q-points : %d instead of %d\n",total_nqpt,N_QPTSNR);
	}

	



	//CALCULATE 2nd order contribution to e-mag scattering
	int n1,n2,n_equiv,total_nqpt_sqr;
	total_nqpt_sqr=0;
	double t1;
	double k_diff[3];
	double gamma_k_k_prime;
	for(n1=0;n1<N_QPTS;n1++){
		for(n2=0;n2<N_QPTS;n2++){

			t1 = -2.*PI*J_sd*J_sd/(S_QN*S_QN*N_QPTSNR);
			t1 = t1*(frequency_k[n1]-frequency_k[n2]);
						
			k_diff[0] = k_pos_index[n1][0] - k_pos_index[n2][0];
			k_diff[1] = k_pos_index[n1][1] - k_pos_index[n2][1];
			k_diff[2] = k_pos_index[n1][2] - k_pos_index[n2][2];

			n_equiv=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);

			t1 = t1*(Lambda_k_uu[n_equiv]+Lambda_k_dd[n_equiv]);
			gamma_k_k_prime = t1*multiplicity_k[n2];
			

			G_emag_2 +=multiplicity_k[n1]* gamma_k_k_prime*C_k[n1]*KB*300.0/(HBAR*frequency_k[n1]);

			total_nqpt_sqr+=multiplicity_k[n1]*multiplicity_k[n2];

			printf("%d %d %e \n",n1,n2,gamma_k_k_prime);
		}
	}
	
	
	printf("%e\n",C_mag/N_QPTSNR);
	printf("%e\n",G_emag_1/N_QPTSNR);
	printf("%e\n",G_emag_2/N_QPTSNR);
	if (total_nqpt_sqr != N_QPTSNR*N_QPTSNR){	
		printf("Wrong number of q-points : %d instead of %d\n",total_nqpt_sqr,N_QPTSNR*N_QPTSNR);
	}

	
	G_emag_2=0.0;
	total_nqpt_sqr=0;
	int i1,i2,j1,j2,k1,k2;
	double k_pos_1[3],k_pos_2[3],k_diff_1[3],k_diff_2[3];
	double gamma_k_k_prime_matrix[N_QPTS][N_QPTS];
	for(n1=0;n1<N_QPTS;n1++){
	for(n2=0;n2<N_QPTS;n2++){
		gamma_k_k_prime_matrix[n1][n2]=0.0;
	}
	}
	//LOOP OVER (FIRST) IR Q-POINT
	for(n1=0;n1<N_QPTS;n1++){
		for(i1=0;i1<N_q_per_line;i1++){
		for(j1=0;j1<N_q_per_line;j1++){
		for(k1=0;k1<N_q_per_line;k1++){	
			//FIND IR Q-POINT OF SECOND Q-POINT
			k_pos_1[0] = i1*1.0/(1.0*N_q_per_line);		k_pos_1[1] = j1*1.0/(1.0*N_q_per_line);		k_pos_1[2] = k1*1.0/(1.0*N_q_per_line);
			n2=get_equivalent_index(k_pos_1,k_pos_index,b1,b2,b3);
			
			//FIND IR Q-POINT FOR THE Q1-Q2
			k_diff[0] = k_pos_index[n1][0] - k_pos_1[0];
			k_diff[1] = k_pos_index[n1][1] - k_pos_1[1];
			k_diff[2] = k_pos_index[n1][2] - k_pos_1[2];
			n_equiv=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);


			gamma_k_k_prime = -2.*PI*J_sd*J_sd/(S_QN*S_QN*N_QPTSNR)*(frequency_k[n1]-frequency_k[n2])*(Lambda_k_uu[n_equiv]+Lambda_k_dd[n_equiv]);
			gamma_k_k_prime_matrix[n1][n2]+=gamma_k_k_prime;

			G_emag_2 +=multiplicity_k[n1]* gamma_k_k_prime*C_k[n1]*KB*300.0/(HBAR*frequency_k[n1]);

			total_nqpt_sqr+=multiplicity_k[n1];

			//printf("%d %d %d %e %e %e %e %e %e %e %e\n",n1+1,n2+1,n_equiv+1,k_diff[0],k_diff[1],k_diff[2],C_k[n1]*KB*300.0/(HBAR*frequency_k[n1]),(frequency_k[n1]-frequency_k[n2]), Lambda_k_uu[n_equiv]-Lambda_k_dd[n_equiv],gamma_k_k_prime,G_emag_2);
		}
		}
		}
	}
	

	for(n1=0;n1<N_QPTS;n1++){
	for(n2=0;n2<N_QPTS;n2++){
		printf("%d %d %e\n",n1,n2,gamma_k_k_prime_matrix[n1][n2]);
	}
	}
	
	printf("%e\n",C_mag/N_QPTSNR);
	printf("%e\n",G_emag_1/N_QPTSNR);
	printf("%e\n",G_emag_2/N_QPTSNR);
	if (total_nqpt_sqr != N_QPTSNR*N_QPTSNR){	
		printf("Wrong number of q-points : %d instead of %d\n",total_nqpt_sqr,N_QPTSNR*N_QPTSNR);
	}
	
	
	//int n3,n4,n5,n6,n7,n8;
	//double prefactor= 4.0*PI*HBAR/(4.0*N_QPTSNR*4.0*N_QPTSNR);
	//double delta_omega;
	//int	total_nqpt_cubed=0;

	/*for(n1=0;n1<N_QPTS;n1++){ //loop over k
		for(n2=0;n2<N_QPTS;n2++){ // loop over k'
			for(n3=0;n3<N_QPTS;n3++){ // loop over q

				//find q' = k + k' - q 
				k_diff[0] = k_pos_index[n1][0] + k_pos_index[n2][0] - k_pos_index[n3][0];
				k_diff[1] = k_pos_index[n1][1] + k_pos_index[n2][1] - k_pos_index[n3][1];
				k_diff[2] = k_pos_index[n1][2] + k_pos_index[n2][2] - k_pos_index[n3][2];

				n4=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);
				delta_omega = frequency_k[n1] + frequency_k[n2] - frequency_k[n3] - frequency_k[n4];
				
				if( (n1==n3 && n2==n4) || (n2==n3 && n1==n4) || fabs(HBAR*delta_omega)>5.0*SWIDTH ) {
					//ONLY LIST COMBINATIONS WHERE NOT THE SAME TWO MAGNONS ARE CREATED, THAT HAVE BEEN ANNIHILATED 
					//ONLY LIST COMBINATIONS WITH REASONABLE CONTRIBUTION DUE TO DELTA DISTRIBUTION TO FACILIATE 	SUBSEQUENT SIMULATIONS
				} else {
						//CALCULATE |GAMMA_kk'q|^2
						//find index for the k'-q s
						k_diff[0] = k_pos_index[n1][0] - k_pos_index[n3][0];
						k_diff[1] = k_pos_index[n1][1] - k_pos_index[n3][1];
						k_diff[2] = k_pos_index[n1][2] - k_pos_index[n3][2];
						n5=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);

						k_diff[0] = k_pos_index[n2][0] - k_pos_index[n4][0];
						k_diff[1] = k_pos_index[n2][1] - k_pos_index[n4][1];
						k_diff[2] = k_pos_index[n2][2] - k_pos_index[n4][2];
						n6=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);	

						k_diff[0] = k_pos_index[n1][0] - k_pos_index[n4][0];
						k_diff[1] = k_pos_index[n1][1] - k_pos_index[n4][1];
						k_diff[2] = k_pos_index[n1][2] - k_pos_index[n4][2];
						n7=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);

						k_diff[0] = k_pos_index[n2][0] - k_pos_index[n3][0];
						k_diff[1] = k_pos_index[n2][1] - k_pos_index[n3][1];
						k_diff[2] = k_pos_index[n2][2] - k_pos_index[n3][2];
						n8=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);	

						t1=prefactor;
						t1=t1*pow(0.5*(frequency_k[n1] + frequency_k[n2] + frequency_k[n3] + frequency_k[n4] - frequency_k[n5]- frequency_k[n6]- frequency_k[n7]- frequency_k[n8]),2.0);
						//Multiply delta distribution
						t1=t1*delta_distribution(HBAR*delta_omega/SWIDTH)/SWIDTH;
						//Multiply multiplicities of n2 and n3 
						t1=t1*multiplicity_k[n2]*multiplicity_k[n3];			
						if (fabs(t1)> 1e11){

						//printf("%d %d %d %d %e %e %e %e %e %e %e %d %e\n", n1+1,n2+1,n3+1,n4+1,frequency_k[n1],frequency_k[n2],frequency_k[n3],frequency_k[n4],HBAR*delta_omega,delta_distribution(HBAR*delta_omega/SWIDTH)/SWIDTH,t1,multiplicity_k[n1],frequency_k[n1]);
						printf("%d %d %d %d %e %d %e\n", n1+1,n2+1,n3+1,n4+1,t1,multiplicity_k[n1],frequency_k[n1]);
						}
					}
					total_nqpt_cubed+=multiplicity_k[n1]*multiplicity_k[n2]*multiplicity_k[n3];			

			}
		}

	}

	if (total_nqpt_cubed != N_QPTSNR*N_QPTSNR*N_QPTSNR){	
		printf("Wrong number of q-points : %d instead of %d\n",total_nqpt_cubed, N_QPTSNR*N_QPTSNR*N_QPTSNR);
	}

	*/
	
	//TEST SYMMETRY OF MAG MAG INTERACTION
	/*int n3,n4,n5,n6,n7,n8;
	double prefactor= 4.0*PI*HBAR/(4.0*N_QPTSNR*4.0*N_QPTSNR);
	double delta_omega;
	double Gamma_k_kp_q_qp[N_QPTS][N_QPTS][N_QPTS][N_QPTS];
	
	for(n1=0;n1<N_QPTS;n1++){
	for(n2=0;n2<N_QPTS;n2++){
	for(n3=0;n3<N_QPTS;n3++){
	for(n4=0;n4<N_QPTS;n4++){
		Gamma_k_kp_q_qp[n1][n2][n3][n4]=0.0;
	}
	}
	}
	}

	total_nqpt_sqr=0;
	int i1,i2,j1,j2,k1,k2;
	double k_pos_1[3],k_pos_2[3],k_diff_1[3],k_diff_2[3];
	for(n1=0;n1<N_QPTS;n1++){
		for(i1=0;i1<N_q_per_line;i1++){
		for(j1=0;j1<N_q_per_line;j1++){
		for(k1=0;k1<N_q_per_line;k1++){	
			k_pos_1[0] = i1*1.0/(1.0*N_q_per_line);		k_pos_1[1] = j1*1.0/(1.0*N_q_per_line);		k_pos_1[2] = k1*1.0/(1.0*N_q_per_line);
			n2=get_equivalent_index(k_pos_1,k_pos_index,b1,b2,b3);
			for(i2=0;i2<N_q_per_line;i2++){
			for(j2=0;j2<N_q_per_line;j2++){
			for(k2=0;k2<N_q_per_line;k2++){
				k_pos_2[0] = i2*1.0/(1.0*N_q_per_line);		k_pos_2[1] = j2*1.0/(1.0*N_q_per_line);		k_pos_2[2] = k2*1.0/(1.0*N_q_per_line);
				n3=get_equivalent_index(k_pos_2,k_pos_index,b1,b2,b3);


				for(int i=0;i<3;i++){	
					k_diff_1[i] = k_pos_index[n1][i] + k_pos_1[i] - k_pos_2[i];
				}

				n4=get_equivalent_index(k_diff_1,k_pos_index,b1,b2,b3);
				delta_omega = frequency_k[n1] + frequency_k[n2] - frequency_k[n3] - frequency_k[n4];

		


				if( (n1==n3 && n2==n4) || (n2==n3 && n1==n4) || fabs(HBAR*delta_omega)>5.0*SWIDTH ) {
					//ONLY LIST COMBINATIONS WHERE NOT THE SAME TWO MAGNONS ARE CREATED, THAT HAVE BEEN ANNIHILATED 
					//ONLY LIST COMBINATIONS WITH REASONABLE CONTRIBUTION DUE TO DELTA DISTRIBUTION TO FACILIATE 	SUBSEQUENT SIMULATIONS
				} else {
						//CALCULATE |GAMMA_kk'q|^2
						//find index for the k'-q s
						k_diff[0] = k_pos_index[n1][0] - k_pos_index[n3][0];
						k_diff[1] = k_pos_index[n1][1] - k_pos_index[n3][1];
						k_diff[2] = k_pos_index[n1][2] - k_pos_index[n3][2];
						n5=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);

						k_diff[0] = k_pos_index[n2][0] - k_pos_index[n4][0];
						k_diff[1] = k_pos_index[n2][1] - k_pos_index[n4][1];
						k_diff[2] = k_pos_index[n2][2] - k_pos_index[n4][2];
						n6=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);

						k_diff[0] = k_pos_index[n1][0] - k_pos_index[n4][0];
						k_diff[1] = k_pos_index[n1][1] - k_pos_index[n4][1];
						k_diff[2] = k_pos_index[n1][2] - k_pos_index[n4][2];
						n7=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);

						k_diff[0] = k_pos_index[n2][0] - k_pos_index[n3][0];
						k_diff[1] = k_pos_index[n2][1] - k_pos_index[n3][1];
						k_diff[2] = k_pos_index[n2][2] - k_pos_index[n3][2];
						n8=get_equivalent_index(k_diff,k_pos_index,b1,b2,b3);	

						t1=prefactor;
						t1=t1*pow(0.5*(frequency_k[n1] + frequency_k[n2] + frequency_k[n3] + frequency_k[n4] - frequency_k[n5]- frequency_k[n6]- frequency_k[n7]- frequency_k[n8]),2.0);
						//Multiply delta distribution
						t1=t1*delta_distribution(HBAR*delta_omega/SWIDTH)/SWIDTH;

						if (fabs(t1)> 1e9){

						//printf("%d %d %d %d %e %e %e %e %e %e %e %d %e\n", n1+1,n2+1,n3+1,n4+1,frequency_k[n1],frequency_k[n2],frequency_k[n3],frequency_k[n4],HBAR*delta_omega,delta_distribution(HBAR*delta_omega/SWIDTH)/SWIDTH,t1,multiplicity_k[n1],frequency_k[n1]);
						//printf("%d %d %d %d %e %d %e\n", n1+1,n2+1,n3+1,n4+1,t1,multiplicity_k[n1],frequency_k[n1]);
						}
						Gamma_k_kp_q_qp[n1][n2][n3][n4]+=t1;
					}

						total_nqpt_sqr++;
		
			}
			}
			}
			printf("%d\n",total_nqpt_sqr);

		}
		}
		}

	}
	
	f=fopen("Gamma_k_kp_q_qp.dat","w");
	for(n1=0;n1<N_QPTS;n1++){
	for(n2=0;n2<N_QPTS;n2++){
	for(n3=0;n3<N_QPTS;n3++){
	for(n4=0;n4<N_QPTS;n4++){
		fprintf(f,"%d %d %d %d %e\n",n1,n2,n3,n4,Gamma_k_kp_q_qp[n1][n2][n3][n4]);
	}
	}
	}
	fclose(f);
	}

	*/
	return 0;
}