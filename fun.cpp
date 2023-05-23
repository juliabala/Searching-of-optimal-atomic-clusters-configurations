#include<vector>
#include<iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include<fstream> 

#include "vectors_rotation_quaternions.cpp"


using namespace std;


#define N  10 //liczba osobnikow w populacji
#define Nr 30 //liczba rodzicow, zawsze parzysta 
#define Ni Nr/2
#define Ns N/2 // liczba zastepowanych osobnikow w populacji 
#define it_max 40 //max liczba iteracji
#define n 15 //liczba atomow w klastrze 
#define n1 n/2
#define n2 n-n1
#define mt 15
#define pmut 0.03
#define dmax 0.05

void mutacje(vector<vector<double>> &v);
vector<vector<vector<double>>> tworzenie_dzieci(vector<vector<vector<double>>> &rodzice);
void sortowanie_z(vector<vector<double>> &v);
void przesuniecie(vector<vector<double>> &v);
void obrot(vector<vector<double>> &v);
vector<double> il_wek(vector<double> x, vector<double> y);
void zapis_do_pliku(vector<vector<vector<double>>> v, ofstream &of);
void sortowanie_energii(vector<vector<vector<double>>> &v);
void turniejowy(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r);
void ruletka(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r);


void edges(vector<vector<double>> v);
double uniform();
void przesun_atomy(int ile,vector<vector<double>> & r);
void krzyzowanie(vector<vector<vector<double>>> r, vector<vector<vector<double>>> &d);
void ruletka_2(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r);
void turniejowy2(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r);
int random_index(int);

void mutacje(vector<vector<double>> &v){
        for(int i=0;i<n;i++){
			double u=(double)rand()/RAND_MAX;
			if(u<=pmut){
				for(int j=0;j<3;j++){
					double u1=(double)rand()/RAND_MAX;
					double dd=2*(u1-0.5)*dmax;
					v[i][j]+=dd;	
				}
				//v[n][0]=1000;
			}
        }
}



vector<vector<vector<double>>> tworzenie_dzieci(vector<vector<vector<double>>> &rodzice){
    vector<vector<vector<double>>> dzieci(Ni);
    for(int k=0;k<Ni;k++){
        dzieci[k] = vector<vector<double>>(n+1);
        for(int i = 0; i < dzieci[k].size()-1; i++){   
            dzieci[k][i] = vector<double>(3);
            for (int j = 0; j < n1; j++){    
                dzieci[k][j] =rodzice[2*k][j];  
            }
            for (int j = n1; j < n; j++){  
                dzieci[k][j]=rodzice[(2*k)+1][j];       
            }
        }
        dzieci[k][n] = vector<double>(1); 
    }
    return dzieci;
}

void sortowanie_z(vector<vector<double>> &v){
    vector<vector<double>> vec2(n+1);
    for(int i = 0; i < n; i++){   
        vec2[i] = vector<double>(3);
    }
    vec2[n]=vector<double>(1);
    for(int i=0;i<n-1;i++){
        for(int j=0;j<n-1-i;j++){
            if(v[j][2]<v[j+1][2]){
                double temp=v[j][2];
                for(int k=0;k<=n;k++){
                    vec2[j]=v[j];
                    v[j]=v[j+1];
                    v[j+1]=vec2[j];    
                }
            }
        }
    }
}
   
void przesuniecie(vector<vector<double>> &v){
    vector<double> vsr(3,0.);
    for(int i=0;i<n;i++){
        vsr[0]+=v[i][0];
        vsr[1]+=v[i][1];
        vsr[2]+=v[i][2];
    }
    vsr[0]/=n;
    vsr[1]/=n;
    vsr[2]/=n;

    for(int i=0;i<n;i++){
        v[i][0]-=vsr[0];
        v[i][1]-=vsr[1];
        v[i][2]-=vsr[2];
    }
}

void obrot(vector<vector<double>> &v){
    vector<double> u(3);
    vector<double> z{0,0,1};
    vector<double> t(3);

	/*
	 * rozklad normalny w 3D  - sferycznie konturowany (Box-Muller) - normowany do 1
	 */
	 
	double len;
	double x1,x2,x3,x4;
	double u1,u2,u3,u4;
	u1=(double)rand()/RAND_MAX;	
	u2=(double)rand()/RAND_MAX;	
	u3=(double)rand()/RAND_MAX;	
	u4=(double)rand()/RAND_MAX;	
	

	x1=sqrt(-2*log(u1))*sin(2.*M_PI*u2);
	x2=sqrt(-2*log(u1))*cos(2.*M_PI*u2);
	x3=sqrt(-2*log(u3))*sin(2.*M_PI*u4);
	len=sqrt(x1*x1+x2*x2+x3*x3);
	x1=x1/len;
	x2=x2/len;
	x3=x3/len;
	
	u[0]=x1;
	u[1]=x2;
	u[2]=x3;
    
    
    
    double kat= angle_between_vectors(u[0],u[1],u[2],z[0],z[1],z[2]);
    t=il_wek(u,z);
    for(int i=0;i<n;i++){
        rotation(t[0],t[1],t[2], kat,v[i][0],v[i][1],v[i][2]);
    }
}



vector<double> il_wek(vector<double> x, vector<double> y){
    vector<double> wynik(3);

    wynik[0]=x[1]*y[2]-y[1]*x[2];
    wynik[1]=x[2]*y[0]-y[2]*x[0];
    wynik[2]=x[0]*y[1]-y[0]*x[1];

    return wynik;
}



void zapis_do_pliku(vector<vector<vector<double>>> v, ofstream &of){
    for(int k=0;k<v.size();k++){
        for(int i = 0; i < v[k].size(); i++){
            for (int j = 0; j < v[k][i].size(); j++){  
                of << v[k][i][j] << " ";
            }
            of<<endl;
        }
        of<<endl;
    }
}







void sortowanie_energii(vector<vector<vector<double>>> &v){
    vector<vector<double>> vec2(n+1);
    for(int i = 0; i < vec2.size()-1; i++){   
        vec2[i] = vector<double>(3);
    }
    vec2[n]=vector<double>(1);
    
    for(int i=0;i<v.size()-1;i++){
        for(int j=0;j<v.size()-1;j++){
            if(v[j][n][0]>v[j+1][n][0]){
                double temp=v[j][n][0];
                for(int k=0;k<=n;k++){
                    vec2=v[j];
                    v[j]=v[j+1];
                    v[j+1]=vec2;   
                }
            }
        }
    }    
}

void sortowanie_energii2(vector<vector<vector<double>>> &v){
    vector<vector<double>> vec2(n+1);
    for(int i = 0; i < vec2.size()-1; i++){   
        vec2[i] = vector<double>(3);
    }
    vec2[n]=vector<double>(1);
    
    for(int i=0;i<v.size()-1;i++){
        for(int j=0;j<v.size()-1;j++){
            if(v[j][n][0]>v[j+1][n][0]){
                double temp=v[j][n][0];
                for(int k=0;k<n;k++){
                    vec2=v[j];
                    v[j]=v[j+1];
                    v[j+1]=vec2;   
                }
            }
        }
    }    
}



void turniejowy(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r){
    for(int i=0;i<Nr;i+=2){
        vector<int> indeksy(mt);
        double best1,best2, best_ind1, best_ind2;
        indeksy[0]= rand()%N;

        for(int j=1;j<mt;j++){
            int temp=0;
            while(temp==0){
                indeksy[j]= rand()%N;
                temp=1;
                for(int k=0;k<j;k++){
                    if(indeksy[j]==indeksy[k]){
                        temp=0;
                        break;
                    }
                }
            }
        }  
        best1=v[indeksy[0]][n][0];
        best2=v[indeksy[mt/2]][n][0];
        for(int j=1;j<mt/2;j++){
            if(best1>v[indeksy[j]][n][0]){
                best_ind1=indeksy[j];
                best1=v[indeksy[j]][n][0];
            }
        }
        for(int j=mt/2+1;j<mt;j++){
            if(best2>v[indeksy[j]][n][0]){
                best_ind2=indeksy[j];
                best2=v[indeksy[j]][n][0];
            }
        }
        r[i]=v[best_ind1];
        r[i+1]=v[best_ind2];
    }
}



/*************************************************************
 * przesuwamy atomy powyzej plaszczyzny x-y
 * 
 * okreslona liczbe atomow (ile) umieszczamy powyzej xy
 * 
 *************************************************************/
void przesun_atomy(int ile,vector<vector<double>> & r){
	double znak;
	vector<double>a(3);
	
	//sortowanie wzgledem wartosci z
	for(int i=0;i<n;i++){
		for(int j=0;j<(n-1);j++){
			if(r[j][2]<r[j+1][2]){
				a=r[j];
				r[j]=r[j+1];
				r[j+1]=a;
			}
		}
	}
	
	double z2=r[ile-1][2];   // sprawdzamy gdzie lezy interesujacy atom
	double dz=0.01;
	for(int i=0;i<n;i++){
		if(i<ile){
			znak=1.0;  //dz powyzej xy
		} else{
			znak=-1.0; //dz ponizej xy
		}
		r[i][2]=r[i][2]-z2+znak*dz; // przesuwamy na zero + znak*dz 
	}
}



/*************************************************************
 * krzyzowanie
 * 
 *************************************************************/

void krzyzowanie(vector<vector<vector<double>>> r, vector<vector<vector<double>>> &d){
	
	// przesuniecie srodka masy + losowy obrot w przestrzeni
	for(int i=0;i<Nr;i++){
            przesuniecie(r[i]);
		obrot(r[i]);
        }
        
     //   przesuniecie n1  atomow wzdluz kierunku osi "z" powyzej z=0
      for(int k=0;k<(Nr-1);k+=2){
		przesun_atomy(n1,r[k]);      //rodzic 1
		przesun_atomy(n1,r[k+1]);    //rodzic 2
	
		int indx1=n1;   //random_index(n); // punkt krzyzowania: rodzic_1_(0:indx1 ) +  rodzic_2_(indx1+1:n-1) 
		
		for(int i=0;i<n;i++){
			if(i<indx1){
					d[k/2][i]=r[k][i];
			}else{
					d[k/2][i]=r[k+1][i];
			}
			
		}
	}  
}

/*************************************************************
 * generator jednorodny
 * 
 *************************************************************/
double uniform(){
	return (double)rand()/RAND_MAX;
}


/*************************************************************
 * losowanie - ruletka
 * 
 *************************************************************/

void ruletka_2(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r){
	int kr=r.size();
	vector<double> wagi(N);
	vector<double> pr(N);
	
	vector<int> indx1(kr);
	vector<int> indx2(kr);
	
	double emin,emax,e,sum,u;
	int m1,m2,ile;
	
	emin=1.E+5;
	emax=-1.E+5;
	
	for(int i=0;i<N;i++){
		e=v[i][n][0];
		if(e<emin)emin=e;
		if(e>emax)emax=e;
	}
	
	sum=0.;
	for(int i=0;i<N;i++){
		e=v[i][n][0];
		wagi[i]=pow(fabs( (emax-e)/(emax-emin)),1.0);	
		sum+=wagi[i];
	}
	
	
	
	for(int i=0;i<N;i++){
		pr[i]=wagi[i]/sum;//normalizacja
		if(i>0)pr[i]=pr[i-1]+pr[i];	
	}
	
	
	pr[N-1]=1.0001; //zabezpieczenie
	
	// dobieramy pary rodzicow
	for(int k=0;k<kr;k+=2){
		
		do{
			
			// 1-rodzic
			for(int i=0;i<N;i++){
				for(int j=0;j<5;j++)u=(double)rand()/RAND_MAX;
				if(u<=pr[i]){
					m1=i;
					break;
				}
			}			
			
			
			// 2-rodzic	
			for(int i=0;i<N;i++){
				for(int j=0;j<5;j++)u=(double)rand()/RAND_MAX;
				if(u<=pr[i]){
					m2=i;
					break;
				}
			}
			
		//sprawdzamy powtarzajace sie pary - takich nie chcemy		
			ile=0;
			if(k>0){
				for(int i=0;i<k;i++){
					if(m1==indx1[i] && m2==indx2[i])ile=1;	
					if(m1==indx2[i] && m2==indx1[i])ile=1;			
				}
			}	
			indx1[k]=m1;	
			indx2[k]=m2;
			
			
		}while(m1==m2 || ile!=0); //drugiego pary az znajdziemy -ale musi byc rozny od pierwszego	
		
		r[k]=v[m1];
		r[k+1]=v[m2];
	
	}
}



/*************************************************************************************************************************
 * zapis krawedzi fullerenu do pliku
 * 
 *************************************************************************************************************************/
void edges(vector<vector<double>> v){
	
	FILE *fp;
	fp=fopen("edges.dat","w");
	
	for(int i=0;i<(n-1);i++){
		for(int j=i+1;j<n;j++){
		
			double xi=v[i][0];
			double yi=v[i][1];
			double zi=v[i][2];
			
			double xj=v[j][0];
			double yj=v[j][1];
			double zj=v[j][2];
			
			double dx=xi-xj;
			double dy=yi-yj;
			double dz=zi-zj;
			double r=sqrt(dx*dx+dy*dy+dz*dz);
			
			if(r>1.2 && r<1.7){
				fprintf(fp,"%15.5E   %15.5E   %15.5E   \n",xi,yi,zi);
				fprintf(fp,"%15.5E   %15.5E   %15.5E   \n\n\n",xj,yj,zj);
			}
		}
	}
	
	
	fclose(fp);
	
	
	
}





/*************************************************************************************************************************
 * losowanie na dwie grupy wykonujemy poprzez K-krotna permutacje tablicy indeksow 
 * 
 *************************************************************************************************************************/

void turniejowy_2(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r){
	int kmax=v.size();
	int kr=r.size();
	vector<int> indx(kmax);
	for(int k=0;k<kr-1;k+=2){
		double emin;
		int m1=-100;
		int m2=-100;
		int l;
		do{
			for(int i=0;i<kmax; i++)indx[i]=i; //pierwotna kolejnosc 	
			//permutacje na tablicy indeksow - mieszanie
			for(int i=0;i<kr*5;i++){
				int j1=random_index(kmax);
				int j2=random_index(kmax);
				int jj=indx[j2];
				indx[j2]=indx[j1];
				indx[j1]=jj;
			}
			
			// i<mt/2  - 1 rodzic (poczatek tablicy)
			emin=1.E+10;
			for(int i=0;i<mt/2;i++){
				l=indx[i];
				if(emin>v[l][n][0]){
					emin=v[l][n][0];
					m1=l;
				}
			}
			
			// i>=N-mt/2  - 2 rodzic (koniec tablicy)
			emin=1.E+10;
			for(int i=(N-mt/2);i<N;i++){
				l=indx[i];
				if(emin>v[l][n][0]){
					emin=v[l][n][0];
					m2=l;
				}
			}
				
			r[k]=v[m1];
			r[k+1]=v[m2];
		}while(m1==m2);   // zabezpieczenie na wypadek gdyby mt bylo za duze
	}
}

/*
 * generujemy losowy indeks calkowity z ciagu [0,1,2,3,...,n-1]
 */
int random_index(int index){
	double x=(double)rand()/(RAND_MAX+1.0);
	double dx=1./index;
	int j=(int)(x/dx);
	return j;
}








void ruletka(vector<vector<vector<double>>> v, vector<vector<vector<double>>> &r){
    vector<double> prob(Nr+1);
double sum;
for(int i=0; i<=Nr;i++){
    prob[i]=i;
    sum+=prob[i];
    }

for(int i=0; i<=Nr;i++){
    prob[i]=prob[i]/sum;
   // cout<<prob[i]<<endl;
}
//cout<<endl<<endl;
for(int j=0; j<Nr;j++){
    prob[j+1]=prob[j]+prob[j+1];
    //cout<<prob[j]<<endl;
}


	


//cout<<endl<<endl;
    for(int i=0;i<Nr;i++){
    double num=(double)rand()/RAND_MAX;
    for(int j=0;j<Nr;j++) {
        if(num>=prob[j]&&num<prob[j+1]){
            r[i]=v[Nr-j+1];
        }
          
    }
}
}
