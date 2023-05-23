/******************************************************************************************************************************
 *  
 * klasa zawiera metody do:
 *
 * 1) liczenia calkowitej energii oddzialywania w klastrze
 *    parametryzacja potencjalu: BRENNER
 *    
 *    przekazujemy tablice 2D [n][3]: n-atomow, 0-x, 1-y, 2-z
 *    
 * 2) minimalizacji energii klastra   
 * 
 *    przekazujemy tablice 2D [n][3]: n-atomow, 0-x, 1-y, 2-z
 *    oraz krok dr do liczenia ilorazu roznicowego (pochodna energii calkowitej)
 * 
 * 
 *------   energie klastra liczymy przy uzyciu porcedury:
 * 
 *       compute_potential_all(const int & n, const vector<vector<double>> & atoms, int itype_potential);
 * 
 * 		n -liczba atomow w klastrze
 * 		atoms[n][3] - polozenia (xi,yi,zi) atomow w klastrze
 *		itype_potential = 1(BRENNER), 2(MORSE)
 * 
 * -----  mnimalizacje energii klastra wykonujemy przy uzyciu procedry
 * 
 * 	local_optimization(const int & n,  vector<vector<double>> & atoms, double dr, int iter_max, double step_size, double tol, int itype_potential);
 * 
 *		n -liczba atomow w klastrze
 * 		atoms[n][3] - polozenia (xi,yi,zi) atomow w klastrze
 * 		dr - krok do liczenia ilorazu roznicowego (pochodne czastkowe)   (np. dr=0.005)
 * 		step_size - wielkosc pierwszego kroku jaki robi procedura  (np. step_size=0.3)
 * 		tol - tolerancja rozwiazania   (np. tol=0.1)
 * 		iter_max - maksymalna liczba iteracji  (np. iter_max=30)
 *		itype_potential = 1(BRENNER), 2(MORSE)
 * 
 ******************************************************************************************************************************/
#include <cmath>
#include <random>
#include <vector>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <functional>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>



using namespace std;

inline double atoms_distance(const vector<vector<double>> & atoms, int i, int j);
inline void   compute_dist_ij(const int & n, const vector<vector<double>> & atoms, vector<vector<double>> & dist_ij);

inline double compute_potential_all(const int & n, const vector<vector<double>> & atoms, int itype_potential);
inline double compute_potential_vi(const int & n, int i, const vector<vector<double>> & atoms, int itype_potential);
inline double compute_potential_vij(const int & n, const vector<vector<double>> & atoms, int i,int j,
							  const vector<vector<double>> & dist_ij,int itype_potential);

inline double compute_bij(const int & n, const vector<vector<double>> & atoms,int i, int j,const vector<vector<double>> & dist_ij);
inline double compute_pot_repulsive(double r);
inline double compute_pot_attractive(double r);
inline double compute_fcut(double r);


void   local_optimization(const int & n,  vector<vector<double>> & atoms, double dr, int iter_max, double step_size, double tol, int itype_potential);
double f_energy(const gsl_vector * xvec, void * params);
void   df_energy(const gsl_vector * xvec, void * params, gsl_vector * df);
void   df_energy_fast(const gsl_vector * xvec, void * params, gsl_vector * df);
void   fdf_energy(const gsl_vector *xvec, void *params, double *f, gsl_vector *df);
void   fdf_energy_fast(const gsl_vector *xvec, void *params, double *f, gsl_vector *df);

class EN_BR
{
	public:
		
	/*
	 *  parametry potencjalu Brennera dla fullerenow
	 */	
	
		static constexpr double R1=1.70;
		static constexpr double R2=2.00;
		static constexpr double De=6.325;
		static constexpr double S=1.29;
		static constexpr double lambda=1.5;
		static constexpr double R0=1.315;
		static constexpr double delta=0.80469;
		static constexpr double a0=0.011304;
		static constexpr double c0=19;
		static constexpr double d0=2.5;
};




class EN_MORSE
{
	public:
		
	/*
	 *  parametry potencjalu MORSE'A
	 */	
	
		static constexpr double De=1.0;
		static constexpr double re=1.0;
		static constexpr double alfa=14.0;
		
};		



/******************************************************************************************************************************
 * ****************************************************************************************************************************
 * ****************************************************************************************************************************
 *  
 * 			  ENERGIA KLASTRA
 * 
 * ****************************************************************************************************************************
 * ****************************************************************************************************************************
 ******************************************************************************************************************************/


/*******************************************************************************
 * liczymy odleglosci miedzy atomami i zapisujemyw tablicy - bedzie szybciej
 * 
 *******************************************************************************/

inline void compute_dist_ij(const int & n, const vector<vector<double>> & atoms, vector<vector<double>> & dist_ij){
	
	double xi,yi,zi, xj,yj,zj,r;
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			r=atoms_distance(atoms,i,j);
			dist_ij[i][j]=r;
			dist_ij[j][i]=r;
		}
	}
	return;
}	
		
inline double atoms_distance(const vector<vector<double>> & atoms, int i, int j){
	return sqrt(pow(atoms[i][0]-atoms[j][0],2)+pow(atoms[i][1]-atoms[j][1],2)+pow(atoms[i][2]-atoms[j][2],2));
}		
		
		
		
		
		
/*******************************************************************************
 * liczymy potencjal calkowity:  1-BRENNER, 2-MORSE (itype_potential)
 * 
 *******************************************************************************/
inline double compute_potential_all(const int & n, const vector<vector<double>> & atoms, int itype_potential)
{
	double res=0.;
	vector<vector<double>> dist_ij(n,vector<double>(n,0.));
	compute_dist_ij(n,atoms,dist_ij);
	
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
				res+= compute_potential_vij(n,atoms,i,j,dist_ij,itype_potential);
		}
	}
	return res;	
}

/*******************************************************************************
 * liczymy potencjal i-tego atomu: 1-BRENNER, 2-MORSE
 *******************************************************************************/
inline double compute_potential_vi(const int & n, int i, const vector<vector<double>> & atoms, int itype_potential)
{
	double res=0.;
	vector<vector<double>> dist_ij(n,vector<double>(n,0.));
	compute_dist_ij(n,atoms,dist_ij);
		for(int j=0;j<n;j++){
			if(i!=j){
				res+=compute_potential_vij(n,atoms,i,j,dist_ij,itype_potential);
			}
		}
	return res;	
}	
	
/*******************************************************************************
 * liczymy potencjal pary atomow:  1-BRENNER, 2-MORSE 
 *******************************************************************************/

inline double compute_potential_vij(const int & n, const vector<vector<double>> & atoms, int i,int j,
							  const vector<vector<double>> & dist_ij, int itype_potential)
{
	double vij=0;
	
	if(itype_potential==1){ //BRENNER
		double xi,yi,zi,xj,yj,zj;
		double ur,ua,bij,fcut;
		double r=dist_ij[i][j];
		if(r<EN_BR::R2){
			fcut=compute_fcut(r);
			ur=compute_pot_repulsive(r);
			ua=compute_pot_attractive(r);
			bij=(compute_bij(n,atoms,i,j,dist_ij)+compute_bij(n,atoms,j,i,dist_ij))/2.;
			vij=fcut*(ur-bij*ua);
		} else{
			vij=0.;
		}
	}else if(itype_potential==2){ //MORSE
		double rij=dist_ij[i][j];
		double wsp=exp(EN_MORSE::alfa*(1-rij/EN_MORSE::re));
		vij=EN_MORSE::De*(wsp*(wsp-2.0));
	}else{
		vij=0.;
	}	

	return vij;
}









/*******************************************************************************
 * liczymy wartosc parametru bij - BRENNER
 * 
 *******************************************************************************/

inline double compute_bij(const int & n, const vector<vector<double>> & atoms,int i, int j,
				const vector<vector<double>> & dist_ij)
{
	double bij=0.;
	double zeta=0.;
	double rik,fik, cos_teta,gteta;
	double xi,yi,zi,xj,yj,zj;
	double xk,yk,zk;
	double rij;
			
	xi=atoms[i][0];
	yi=atoms[i][1];
	zi=atoms[i][2];
	
	xj=atoms[j][0];
	yj=atoms[j][1];
	zj=atoms[j][2];
	
	//rij=sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));
	rij=dist_ij[i][j];
	zeta=0.;
	for(int k=0;k<n;k++){
		if(k!=i && k!=j){
			xk=atoms[k][0];
			yk=atoms[k][1];
			zk=atoms[k][2];
			//rik=sqrt(pow(xi-xk,2)+pow(yi-yk,2)+pow(zi-zk,2));
			rik=dist_ij[i][k];
			if(rik<=EN_BR::R2){
				fik=compute_fcut(rik);
					cos_teta=((xj-xi)*(xk-xi)+(yj-yi)*(yk-yi)+(zj-zi)*(zk-zi))/rij/rik;
					gteta=EN_BR::a0*(1+pow(EN_BR::c0/EN_BR::d0,2)-pow(EN_BR::c0,2)/(pow(EN_BR::d0,2)+pow(1+cos_teta,2)));
					zeta+=fik*gteta;
			}
		}
	}

	bij=pow(1+zeta,-EN_BR::delta);
	return bij;
}



/*******************************************************************************
 * liczymy wartosc potencjalu odpychania UR -BRENNER
 * 
 *******************************************************************************/

inline double compute_pot_repulsive(double r){
		double res=0.;
		res=EN_BR::De/(EN_BR::S-1)*exp(-sqrt(2*EN_BR::S)*EN_BR::lambda*(r-EN_BR::R0));
		return res;
	}

/*******************************************************************************
 * liczymy wartosc potencjalu przyciagania UA -BRENNER
 * 
 *******************************************************************************/

inline double compute_pot_attractive(double r){
		double res=0.;
		res=EN_BR::De*EN_BR::S/(EN_BR::S-1)*exp(-sqrt(2/EN_BR::S)*EN_BR::lambda*(r-EN_BR::R0));
		return res;
	}


/*******************************************************************************
 * liczymy wartosc funkcji odciecia  - BRENNER
 * 
 *******************************************************************************/

inline double compute_fcut(double r){
		double res;
		if(r<=EN_BR::R1){
			res=1.;
		}else if( r> EN_BR::R1 && r<= EN_BR::R2){
			res=(1+cos(M_PI*(r-EN_BR::R1)/(EN_BR::R2-EN_BR::R1)))/2.;
		}else{
			res=0.;
		}
		return res;
	}

	
	
	
	
	
	
	
	
/******************************************************************************************************************************
 * ****************************************************************************************************************************
 * ****************************************************************************************************************************
 *  
 * 			2)   MINIMALIZACJA ENERGII KLASTRA
 * 
 * 
 * procedury do lokalnej minimalizacji energii metoda CG 
 * 
 * przekazujemy:
 * 			tablice 2D [n][3]: n-atomow, 0-x, 1-y, 2-z
 * 			dr(=dx,dy,dz) - krok do liczenia ilorazu roznicowego energii
 * 
 * ****************************************************************************************************************************
 * ****************************************************************************************************************************
 ******************************************************************************************************************************/





		
void local_optimization(const int & n,  vector<vector<double>> & atoms, double dr, int iter_max, double step_size, double tol, 
				int itype_potential)
{
	
	int ndim=n*3; //wymiar problemu: n atomow w klastrze * 3 wymiary przestrzenne
	
	
/* wektor startowy + wektor rozwiazan + wektor gradientu*/
	gsl_vector *xvec=gsl_vector_alloc(ndim);
	
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			gsl_vector_set(xvec,i*3+j,atoms[i][j]);
		}
	}
	

/*  ustawianie procedury minimalizacyjnej */	
	
	double *par;
	par=(double*)malloc(10*sizeof(double));
	par[0]=n;
	par[1]=dr;
	double rpot_max=4.0; //zasieg potencjalu 2*rcut
	par[2]=rpot_max;
	par[3]=itype_potential;
	
	
	gsl_multimin_function_fdf my_func;
		my_func.n=ndim;
		my_func.f= &f_energy;
		//my_func.df= &df_energy_fast;
		my_func.df= &df_energy_fast;
		my_func.fdf= &fdf_energy;
		my_func.params=(void *)par;
	
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;	
	//T = gsl_multimin_fdfminimizer_conjugate_fr; // algorytm
	T=gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc (T, ndim); //stan procedury
	
	gsl_multimin_fdfminimizer_set(s,&my_func,xvec,step_size,tol);

	
	
/* iteracyjne poszukiwanie minimum */	

	int iter = 0;
	int status;	
	
	
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
		if (status)	break;
		status=gsl_multimin_test_gradient(s->gradient,1E-10);
		double energy=gsl_multimin_fdfminimizer_minimum(s);
		//printf("%10d   %15.5E\n",iter,energy);
	}
	while (status == GSL_CONTINUE && iter < iter_max);

	
/* zachowujemy nowe polozenia atomow w klastrze  */	
	
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			atoms[i][j]=gsl_vector_get(s->x,i*3+j);
		}
	}
	
	
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(xvec);

	return ;
	
}		




/***************************************************************************
 * energia klastra  
 * 
 ***************************************************************************/
	
double f_energy(const gsl_vector * xvec, void * params){
	double *p = (double *)params;
	int n=(int)(p[0]+0.1);
	int itype_potential=(int)(p[3]+0.1);
	vector<vector<double>> atoms;
	atoms.resize(n,vector<double>(3,0.));
	
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			atoms[i][j]=gsl_vector_get(xvec,i*3+j);
		}
	}
	 
	double en=compute_potential_all(n,atoms,itype_potential);
	return en;
}
	
	
	
/***************************************************************************
 * wektor pochodnych energii klastra  
 * 
 ***************************************************************************/	
	
void df_energy(const gsl_vector * xvec, void * params, gsl_vector * df){
	double *p = (double *)params;
	int n=(int)(p[0]+0.1);
	int itype_potential=(int)(p[3]+0.1);
	double dr=p[1];
	vector<vector<double>> atoms(n,vector<double>(3,0.));
	
	/*
	 * wczytujemy dane do wektora 2D
	 * 
	 */
	
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			atoms[i][j]=gsl_vector_get(xvec,i*3+j);
		}
	}
	
	
	/*
	 * pochodna liczona z ilorazu roznicowego
	 * 
	 */
	
	
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			double old_value=atoms[i][j];
			double rp=old_value+dr;	
			double rm=old_value-dr;
			
			atoms[i][j]=rp;
			double en_p=compute_potential_all(n,atoms,itype_potential);
			
			atoms[i][j]=rm;
			double en_m=compute_potential_all(n,atoms,itype_potential);
			
			double den_dr=(en_p-en_m)/2./dr;
			gsl_vector_set(df,i*3+j,den_dr);
			
			atoms[i][j]=old_value;
		}
	}
	return;
}
		
		

		
		
		
/***************************************************************************
 * 		POCHODNE LICZNE SZYBKO   - FAST
 * 
 * wektor pochodnych energii klastra  
 * 
 ***************************************************************************/	
	
void df_energy_fast(const gsl_vector * xvec, void * params, gsl_vector * df){
	double *p = (double *)params;
	int n=(int)(p[0]+0.1);
	int itype_potential=(int)(p[3]+0.1);
	double dr=p[1];
	double rpot_max=p[2];
	vector<vector<double>> atoms(n,vector<double>(3,0.));
	
	vector<vector<double>> dist_ij(n,vector<double>(n,0.));
	
	
	
	/*
	 * wczytujemy dane do wektora 2D
	 * 
	 */
	
	for(int i=0;i<n;i++){
		for(int j=0;j<3;j++){
			atoms[i][j]=gsl_vector_get(xvec,i*3+j);
		}
	}
	
	
	
	
	// szukamy najblizszych sasiednich atomow rij<rpot_max i liczmy oddzialywanie tylko z nimi
	
	int ile;
	vector<vector<int>> neigh(n,vector<int>(n+1,0));
	vector<vector<double>> atom_energy(n,vector<double>(3,0.0));
	vector<double> ri(3,0.0);
	
	// ---- liczymy stare odleglosci ---------
	compute_dist_ij(n,atoms,dist_ij);
	
	for(int i=0;i<n;i++){
		ile=0;
		for(int j=0;j<n;j++){
			if(i!=j){
				if(dist_ij[i][j]<rpot_max){
					double eij=compute_potential_vij(n,atoms,i,j,dist_ij,itype_potential);
					atom_energy[i][0]+=eij;
					ile++;
					neigh[i][0]=ile; //liczba sasiadow
					neigh[i][ile]=j; //numer sasiada
				}
			}		
		}
	}
	
	
	//energia calkowita z tablicy energii atomow
	double etot=0.;
	for(int i=0;i<n;i++){
		etot+=atom_energy[i][0]/2; //energia liczona podwojnie (vij+vji) - dzielimy przez 2
	}
	
	
	// liczymy pochodne uwzgledniajac tylko lokalna zmiane polozenia atomu
	
	for(int i=0;i<n;i++){
		ri[0]=atoms[i][0];
		ri[1]=atoms[i][1];
		ri[2]=atoms[i][2];
		
		double eip,eim,den_dr;
		
		//kolejnosc pochodnych: dx,dy,dz --> 0,1,2
		for(int ktory_dr=0; ktory_dr<3; ktory_dr++){
					
			//krok do przodu
			atoms[i][ktory_dr]=ri[ktory_dr]+dr;
			
			// ---- liczymy nowe odleglosci ---------
			compute_dist_ij(n,atoms,dist_ij);
			
			
			eip=etot-atom_energy[i][0]/2; //usuwamy wklad od atomu i-tego
			for(int j=1;j<=neigh[i][0];j++){
				int jj=neigh[i][j]; //indeks sasiada
				eip+=compute_potential_vij(n,atoms,i,jj,dist_ij,itype_potential)/2; //liczymy nowy wklad do atomu i-tego
				eip=eip-atom_energy[jj][0]/2;//usuwamy STARY wklad od atomu j-tego
				
				for(int m=1;m<=neigh[jj][0];m++){
					int ll=neigh[jj][m];// indeks sasiada-sasiada jj-tego
					eip+=compute_potential_vij(n,atoms,jj,ll,dist_ij,itype_potential)/2; //liczymy nowy wklad od atomu j-tego
				}
			}
			atoms[i][ktory_dr]=ri[ktory_dr]; //przywracamy stare polozenie
			
			
			//krok do tylu
			atoms[i][ktory_dr]=ri[ktory_dr]-dr;
			
			// ---- liczymy nowe odleglosci ---------
			compute_dist_ij(n,atoms,dist_ij);
			
			eim=etot-atom_energy[i][0]/2; //usuwamy wklad od atomu i-tego
			for(int j=1;j<=neigh[i][0];j++){
				int jj=neigh[i][j]; //indeks sasiada
				eim+=compute_potential_vij(n,atoms,i,jj,dist_ij,itype_potential)/2; //liczymy nowy wklad do atomu i-tego
				eim=eim-atom_energy[jj][0]/2;//usuwamy wklad od atomu j-tego
				
				for(int m=1;m<=neigh[jj][0];m++){
					int ll=neigh[jj][m];// indeks sasiada-sasiada
					eim+=compute_potential_vij(n,atoms,jj,ll,dist_ij,itype_potential)/2; //liczymy nowy wklad od atomu j-tego
				}
			}
			atoms[i][ktory_dr]=ri[ktory_dr]; //przywracamy stare polozenie
			
			den_dr=(eip-eim)/2./dr;
			gsl_vector_set(df,i*3+ktory_dr,den_dr);
		}
	}	
	return;
}
		
		


		
		
		
		
/***************************************************************************
 * energia i wektor pochodnych energii klastra  
 * 
 ***************************************************************************/

void fdf_energy(const gsl_vector *xvec, void *params, double *f, gsl_vector *df){
	*f=f_energy(xvec,params);
	df_energy(xvec,params,df);
	return;
}		


/***************************************************************************
 * energia i wektor pochodnych energii klastra   - FAST
 * 
 ***************************************************************************/

void fdf_energy_fast(const gsl_vector *xvec, void *params, double *f, gsl_vector *df){
	*f=f_energy(xvec,params);
	df_energy_fast(xvec,params,df);
	return;
}

