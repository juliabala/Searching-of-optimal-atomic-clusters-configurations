		

/*************************************************************************************************
 * 
 * generujemy wielokaty dla GNUPLOTA
 * idziemy po sasiadach i wybieramy najmniejszy kat miedzy odcinkami skierowany clockwise
 * po przejsciu przez sasiada kasujemy go, aby wiecej go od tej strony nie odwiedzac
 * sasiedzi to te atomy ktore sa blizej niz rmax
 * tablica z polozeniami atomow znana
 * 
 * 	n - liczba atomow 
 * 	rmax - ponizej tej wartosci atomy sa traktowane jako sasiedzi
 * 	atom[n][3]:   xi=atom[i][0],  yi=atom[i][1],  zi=atom[i][2]  
 * 	plik - nazwa pliku do zapisu wielokatow 
 * 
 * 
 * w GNUPLOCIE rysunek wykonujemy tak:
 * set xyplane 0
 * set view equal
 * set pm3d depthorder border lw 2 
 * set style fill transparent solid 0.3
 * splot 'plik' u 1:2:3 w polygons fc "gold" 
 * 
 *************************************************************************************************/


#include <cstdlib>
#include<cmath>
#include<vector>
#include<stdio.h>
#include<stdlib.h>

using namespace std;

void write_polygons_from_atoms(double rmax, int nat, vector<vector<double>> atom, const char * plik){
	
	FILE *fp=fopen(plik,"w");
	int nemax=20;
	double alfa_min;
	int lmin;
	double r;
	double x0,y0,z0;
	double xi,yi,zi;
	double xj,yj,zj;
	double xa,ya,za;
	double xb,yb,zb;
	double xc,yc,zc;
	double xnb,ynb,znb;
	double xab,yab,zab;
	double xcb,ycb,zcb;
	double xab2,yab2,zab2;
	double xcb2,ycb2,zcb2;
	double rab,rab2;
	double rcb,rcb2;
	double alfa;
	double dxn,dyn,dzn;
	int ia,ib,ic,l,ll,j,jj,i,ipoz,last_poz;
	double cos_alfa,sin_alfa;
	double scalar_product;
		
	int ile,ic_min;
	
	vector<vector<int>> neigh;
	neigh.resize(nat,vector<int>(nemax,0));
	
	vector<vector<double>> plane_vector;
	plane_vector.resize(nat,vector<double>(3,0));
	
	/*
	 * szukamy sasiadow i zapisujemy do tablicy
	 */	
	for(i=0;i<nat;i++){
		xi=atom[i][0];
		yi=atom[i][1];
		zi=atom[i][2];
		for(j=0;j<nat;j++){
			if(j!=i){
				xj=atom[j][0];
				yj=atom[j][1];
				zj=atom[j][2];
				r=sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));
				if(r<=rmax){
					ile=neigh[i][0]+1;
					neigh[i][0]=ile;
					neigh[i][ile]=j;
				}
			}
		}
	}
	
	
	
	/*
	 * dla kazdego atomu tworzymy plaszczyzne na ktora bedziemy 
	 * rzutowac wiazania miedzyatomowe
	 * potrzebujemy tylko wektora normalnego (unormowanego) 
	 * skierowanego od plaszczyzny w strone i-tego atomu
	 * 
	 */	
	
	for(i=0;i<nat;i++){
		xi=atom[i][0];
		yi=atom[i][1];
		zi=atom[i][2];
		
		x0=0.;
		y0=0.;
		z0=0.;
		
		for(int jj=1;jj<=neigh[i][0];jj++){
			j=neigh[i][jj];
			xj=atom[j][0];
			yj=atom[j][1];
			zj=atom[j][2];
			/*
			 * punkt pod atomem z usrednienia wiazan z sasiadami
			 */
			x0=x0+(xj-xi);
			y0=y0+(yj-yi);
			z0=z0+(zj-zi);
		}
		/*
		 * wektor normalny (unormowany) skierowany do atomu i-tego
		 */
		dxn=xi-x0;
		dyn=yi-y0;
		dzn=zi-z0;
		r=sqrt(pow(dxn,2)+pow(dyn,2)+pow(dzn,2));
		plane_vector[i][0]=dxn/r;
		plane_vector[i][1]=dyn/r;
		plane_vector[i][2]=dzn/r;
		
	}
		
	
	/*
	 * tworzymy wielokaty i zapisujemy do pliku
	 */	
	
	for(i=0;i<nat;i++){
		
		/*
		 * tworzymy pojedyncza sciezke po wielokacie
		 * zaczynamy od ostatniego w kolejce sasiada
		 * po zakonczeniu scezki sasiada kasujemy
		 * 
		 */
		while(neigh[i][0]>0){
			ia=i;                        //punkt startowy A
			j=neigh[ia][0]; 
			ib=neigh[ia][j];             //punkt startowy B
			
			if(ib==i){
				ipoz=neigh[i][1];
				neigh[i][1]=ib;
				neigh[i][j]=ipoz;
				ib=ipoz;
			}
			
		
			
		
		
			/*
			 * przesuwamy punkt A do B, B do C itd
			 */
			while(neigh[ib][0]>0){
				
				/* atom A - aktualny bo go zmieniamy*/
					xa=atom[ia][0];
					ya=atom[ia][1];
					za=atom[ia][2];
				
				/* sasiad A: B  - lecimy po kolei dopoki sa sasiedzi*/
					xb=atom[ib][0];
					yb=atom[ib][1];
					zb=atom[ib][2];	
					
					/*wektor normalny do plaszczyzny dla atomu B*/
					xnb=plane_vector[ib][0];
					ynb=plane_vector[ib][1];
					znb=plane_vector[ib][2];
					
					
					/* sasiad B:  C   - szukamy tego ktorego wiazanie r_BC daje najmniejszy kat z wiazaniem r_AB*/	
						ic_min=0;
						alfa_min=1.99*M_PI;				
					for(int l=1; l<=neigh[ib][0];l++){
						ic=neigh[ib][l];
						if(ic==ia){
							continue; //nie wracamy do A
						}
						xc=atom[ic][0];
						yc=atom[ic][1];
						zc=atom[ic][2];	
					
						/*wektor do rab=ra-rb*/
						xab=xa-xb;
						yab=ya-yb;
						zab=za-zb;
						/* rzut na plaszczyzne scentrowana na atomie B*/
						scalar_product=xab*xnb+yab*ynb+zab*znb;
						xab2=xab-scalar_product*xnb;
						yab2=yab-scalar_product*ynb;
						zab2=zab-scalar_product*znb;
						rab2=sqrt(pow(xab2,2)+pow(yab2,2)+pow(zab2,2));
						
						/*wektor do rcb=rc-rb*/
						xcb=xc-xb;
						ycb=yc-yb;
						zcb=zc-zb;
						/* rzut na plaszczyzne scentrowana na atomie B*/
						scalar_product=xcb*xnb+ycb*ynb+zcb*znb;
						xcb2=xcb-scalar_product*xnb;
						ycb2=ycb-scalar_product*ynb;
						zcb2=zcb-scalar_product*znb;
						rcb2=sqrt(pow(xcb2,2)+pow(ycb2,2)+pow(zcb2,2));
						
						/* liczymy kat miedzy wektorami na plaszczyznie */
						cos_alfa=(xab2*xcb2+yab2*ycb2+zab2*zcb2)/rab2/rcb2; 
						sin_alfa=((yab2*zcb2-zab2*ycb2)*xnb+(zab2*xcb2-xab2*zcb2)*ynb+(xab2*ycb2-yab2*xcb2)*znb)/rab2/rcb2;
						
						if(sin_alfa>=0){
							alfa=acos(cos_alfa);
						}else{
							alfa=2.*M_PI-acos(cos_alfa);
						}
					
						if(alfa_min>alfa){
							alfa_min=alfa;
							ic_min=ic;
						}
						
					}
					
					
					
					/*
					* mamy dane 3 atomy: A,B,C
					* zapisujemy do pliku odcinek: AB
					* przesuwamy atomy: A->B i powtarzamy procedure do momentu az ic=i (poczatkowy atom A)
					* 
					*/
					fprintf(fp,"%15.5E   %15.5E   %15.5E \n",xa,ya,za);
					fprintf(fp,"%15.5E   %15.5E   %15.5E   \n",xb,yb,zb);
					
					/* z listy A kasujemy sasiada B*/
					for(int l=1;l<=neigh[ia][0];l++){
						if(ib==neigh[ia][l])ipoz=l;				
					}
						last_poz=neigh[ia][0];
						neigh[ia][ipoz]=neigh[ia][last_poz]; //przepisujemy ostatnia pozycje na wczesniejsze miejsce
						neigh[ia][0]=neigh[ia][0]-1;         //zmniejszamy liczbe sasiadow
					
					
					
					
					/* przesuwamy punkty: A->B, B->C */
						ia=ib;
						ib=ic_min;
					/* jesli wrocilismy do punktu i-tego to wychodzimy ze sciezki - jest juz zamknieta */
						if(ib==i){
							xc=atom[ic_min][0];
							yc=atom[ic_min][1];
							zc=atom[ic_min][2];
							fprintf(fp,"%15.5E   %15.5E   %15.5E \n",xb,yb,zb);
							fprintf(fp,"%15.5E   %15.5E   %15.5E \n",xc,yc,zc);
							fprintf(fp,"\n\n");
							
							/* z listy A kasujemy sasiada B*/
							for(int l=1;l<=neigh[ia][0];l++){
								if(ib==neigh[ia][l])ipoz=l;				
							}
								last_poz=neigh[ia][0];
								neigh[ia][ipoz]=neigh[ia][last_poz]; //przepisujemy ostatnia pozycje na wczesniejsze miejsce
								neigh[ia][0]=neigh[ia][0]-1;         //zmniejszamy liczbe sasiadow
							
							
							break;	
						}
		
			}
			
		}
		
	}
	
	fclose(fp);
	
	return;
}

