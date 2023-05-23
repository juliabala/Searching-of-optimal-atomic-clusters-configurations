#include <omp.h>
#include "energy_fullerene_brenner.cpp"
#include"fun.cpp"
#include"write_polygons_from_atoms.cpp"



/*  KOMPILACJA Z obsluga wielowatkowosci
 * 
 *                       g++ -Ofast -fopenmp  main.cpp -lgsl -lgslcblas -lm
 * 
 */




int main(){

    srand(time(NULL));

    if(Nr%2!=0){
        cout<<"liczba rodzicow musi byc nieparzysta"<<endl;
        return 0;
    }

    ofstream polozenia("polozenia.txt"); //losowe polozenia
    ofstream energie("en.txt"); //energie dla losowych polozen atomow
    ofstream best("best.txt");


    vector<vector<vector<double>>> rodzice(Nr,vector<vector<double>>(n+1,vector<double>(3))); 
    vector<vector<vector<double>>> dzieci(Ni,vector<vector<double>>(n+1,vector<double>(3))); 
    vector<vector<vector<double>>> temp(Ns,vector<vector<double>>(n+1,vector<double>(3))); 
    vector<vector<vector<double>>> vec(N);
    
#pragma omp parallel shared(vec,polozenia,energie)
{
#pragma omp for  nowait
    for(int k=0;k<N;k++){
        vec[k] = vector<vector<double>>(n+1);
        for(int i = 0; i < vec[k].size()-1; i++){   
            vec[k][i] = vector<double>(3);
            for (int j = 0; j < vec[k][i].size(); j++){  
                vec[k][i][j] =uniform();
                polozenia << vec[k][i][j] << " ";
            }    
            polozenia << endl;       
        } 
        polozenia << endl; 
        vec[k][n] = vector<double>(1); 

        local_optimization(n, vec[k], 0.005,30 , 0.1, 0.1, 1);
	  
        double pot;
        pot=compute_potential_all(n, vec[k], 1);
        vec[k][n][0]=pot;
        energie<<pot<<endl;
	  
    }
}
    

    
	  
	  
    int it=0;
    while(it!=it_max){
        if(n%2)
            sortowanie_energii2(vec);
        else    
            sortowanie_energii(vec);
	  
	  //turniejowy_2(vec, rodzice); // jesli malo atomow to turniejowy jest zbyteczny
	  ruletka_2(vec, rodzice);      // prawdopodobienstwo okreslane na podstawie szerokosci przedzialu energii pi=(Emax-Ei)/(Emax-Emin)/Suma(wag)
	  
	  
       //turniejowy(vec, rodzice);
      //ruletka(vec, rodzice);
       // zapis_do_pliku(rodzice, parents);

        krzyzowanie(rodzice,dzieci); // w procedurze przesuwamy rodzicow (srodek masy), obracamy i przesuwamy gora/dol, krzyzujemy 

	  
	#pragma omp parallel shared(dzieci)
	  {
	#pragma omp for  nowait
		for(int i=0; i<Ni;i++){
			local_optimization(n, dzieci[i], 0.005, 20 , 0.1, 0.1, 1);
			dzieci[i][n][0]=compute_potential_all(n, dzieci[i], 1);
		}
	}

	sortowanie_energii(dzieci);
	
	
	
//podmieniamy ostanie Ns osobnikow od konca na najlepsze osobniki z nowej generacji
         for(int i=0; i<Ns;i++)vec[N-1-i]=dzieci[i]; 
	   
	   
//mutujemy tylko stare osobniki poza najlepszym - nowa generacje zostawiamy w spokoju
         for(int i=1; i<(N-Ns);i++) {
		   mutacje(vec[i]);
		   double pot=compute_potential_all(n, vec[i], 1);
		   vec[i][n][0]=pot;
	}
        
// okreslamy energie po mutacji - zmienia sie
	for(int i=1; i<N;i++)	vec[i][n][0]=compute_potential_all(n, vec[i], 1);
        sortowanie_energii(vec);
	
// do pliku zapisujemy najlepsza i najgorsza energie - ich roznica bedzie sie zmniejszac w trakcie symulacji	
        best<<it<<"   "<<vec[0][n][0]<<"   "<<vec[N-1][n][0]<<endl;
        
// biezaca kontrola wyniku 	  
	if(it%10==0){
		printf("it, Emin=   %d  %15.5E\n",it,vec[0][n][0]);
		
		FILE *fp;
		fp=fopen("konf.dat","w");
		for(int i=0;i<n;i++){
			fprintf(fp,"%15.5E   %15.5E   %15.5E   %15.5E   %15.5E   %15.5E   \n",
				  vec[0][i][0],vec[0][i][1],vec[0][i][2],
				  vec[N-1][i][0],vec[N-1][i][1],vec[N-1][i][2]);
		}
		fclose(fp);
		//zapisujemy krawedzie najlepszego osobnika do pliku "egdes.dat" - do rysowania
		edges(vec[0]); 
		
		if(it>200){
			printf("liczymy polygony -  jesli sie zawiesi to trzeba policzyc je po samouzgodniueniu fullerenu\n");
			 /* w GNUPLOCIE rysunek wykonujemy tak:
			  * set xyplane 0
			  * set view equal
			  * set pm3d depthorder border lw 2 
			  * set style fill transparent solid 0.3
			  * splot 'plik' u 1:2:3 w polygons fc "gold"
			  */
			// write_polygons_from_atoms(1.7,n,vec[0],"polygons.dat");
		}
	}
		  
        it++;
    }
    write_polygons_from_atoms(1.7,n,vec[0],"polygons.dat");
    polozenia.close();
    energie.close();
    best.close();

    return 0;
    }




