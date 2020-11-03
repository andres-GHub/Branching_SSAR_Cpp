//~     SSAR_branching.cpp is a heavily modified version of simPL.cpp by Jordi
//~		(the former essentially maintains certain structure/skeleton of the latter
//~		but SSAR_branching.cpp is a completely different model)
//~
//~     Copyright 2017 Andres Zambrano
//~ 	This program is free software: you can redistribute it and/or modify
//~ 	it under the terms of the GNU General Public License as published by
//~ 	the Free Software Foundation, either version 3 of the License, or
//~ 	(at your option) any later version.

//~ 	This program is distributed in the hope that it will be useful,
//~ 	but WITHOUT ANY WARRANTY; without even the implied warranty of
//~ 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//~ 	GNU General Public License for more details.

//~ 	You should have received a copy of the GNU General Public License
//~ 	along with this program.  If not, see <https://www.gnu.org/licenses/>.
//~ 

//		--- Original Code copyright and license ---
//		simPL.cpp
//		
//      Copyright 2015 Jordi <jordi@jordi-laptop>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <list>
#include <stdarg.h>
#include <algorithm>
#include <sstream>
#include <stdio.h>
#include <cmath>
using namespace std;

const double pi = 3.14159265358979323846;

// ********* random seeds   ************* 
int iseed1= 12343, iseed2=67890;
int idJob=0;

//********** global variables ***********

//double Nb;
double tauE;

double Mmin;
list<pair<double,double> > avals; //(time, magnitude) to be rearenged
list<double> list_dM;

double period;

double lamb_max;
double lamb;

double lamb_back;
double k_prod, sigma;

double p;
double c_0;
double tau_0;
double b;
double g;
double z;
double b_as;
double c;
double tau;
double alpha;
double NewMag;

void New_Aval(double t, double M); // see function at end of code

double Back_Rate( double t);

void AS_Rate(double t, double M);

long N_as(double t, double M ); // sample the number of as for each mainshock
double Dt_AS(); // give a DT for each aftershcok since the mainshock (need to renormalise for scaling factor afterwards)        
double T_BG(); // give a DT for each aftershcok since the mainshock (need to renormalise for scaling factor afterwards)

// Branching process:
double New_Branch(double t, double M);

double RAND(), GAUSS();

//****************************************************************
//****************************************************************
    
    
int main(int argc, char** argv)
{
   Mmin      = atof(argv[1]); // Mmin de la distro
   lamb_back = atof(argv[2]); // produccio Background
   period    = atof(argv[3]); // periode de registre
   b         = atof(argv[4]);
   p         = atof(argv[5]); // omori exponet
   c_0       = atof(argv[6]);
   tau_0     = atof(argv[7]);
   g         = atof(argv[8]);
   z         = atof(argv[9]);
   idJob     = atoi(argv[10]);// simulation seed
   
   iseed1   += 2*(idJob);
 
   b_as  = g+z;
   
   //alpha = z+p*g;
   alpha = b_as + g*(p-1);
    
    //char fullPath[128];
 	ofstream fout;

   cout << "  Mmin = " << Mmin << "\nlamb_backgrnd = "<< lamb_back << "  Period = " << period;
   cout << "\nb= " << b << "\nb_as = " << b_as << "\np = " << p << "\nalpha = " << alpha << "\n";
   cout << "\nc_0 =" << c_0 << "\ntau_0 = " << tau_0<< "\ng = " << g << "\nz =" << z <<"\n";
   //sprintf(fullPath, "results/SSAR_N%4.2f.Seq", Nb ) ;
	
	//sprintf(fullPath, "results/SSAR_Sim.Seq" ) ;
	//fout.open(fullPath, ios::out);

    fout.precision(14);

    //k_prod = -c_0/(b_as * tau_0 *(p-1.0e0)*log(10.0e0));    
    //k_prod = Nb*pow(10.0e0,-b_as*Mmin);
    
    k_prod = pow(c_0, p)/tau_0;
    
    cout<<"\nk = " << k_prod;
  
//----------------------------------------------------------------------

// sample  the initial Background Poisson Process:

	long N = round(lamb_back*period);
	
	for (int i = 0; i<N; i++)
	{   	
		NewMag = Mmin-((1.e0/b)*log10(RAND()));
		New_Aval(T_BG(), NewMag);        
	}
	cout << "\nsampled background events: " << N;
	
// explore each event and generate it's own aftershock sequence   
    
   //list<tuple<double,double,double> >::iterator iter=avals.begin();
   list<pair<double,double> >::iterator iter=avals.begin();
   
   char fullPath2[128];

		sprintf(fullPath2, "results/SSAR_Sim_dm.Seq");
			fout.open(fullPath2, ios::out);

	
	while(iter != avals.end())
    {
		New_Branch(iter->first, iter->second);
        iter++;
    }
    fout << "dm"<<"\n";
    
    //------------Original form:
    //for(list<type>::iterator iter = list.begin(); iter != list.end(); iter++){
   //cout<<*iter<<endl;}
   //---------------------------
    
    for(list<double>::iterator iter = list_dM.begin(); 
    iter != list_dM.end(); iter++)
	{
		fout<<*iter<<endl;
	}
    fout.close();
    
  	cout << "\ntotal number: " << avals.size() << "  sorting ...";  
    
     char fullPath[128];
     
		sprintf(fullPath, "results/SSAR_Sim.Seq");
			fout.open(fullPath, ios::out);

    iter=avals.begin();
	cout << "\nrighting...";
	fout << "time	" << "M" << "\n";
	while(iter != avals.end())
    {
        fout << iter->first << "\t" << iter->second << "\n";
        iter++;
    }
    fout.close();
    return 0;
}
//*********************************************************************
//************************ Functions **********************************

long N_as(double t, double M ) // number of aftershocks foreach MS
{	
    long double N_0 =  k_prod*pow(10.0e0, alpha*M);
	cout <<"\n N = " << N_0 << "\n";

//  sample a Poisson number
    long double U = RAND();
    long k = -1;
    long double kk=N_0;
    while(U>=0)
    {
        k++;
        kk += log(max(k,long(1)));
        U -= exp(k*log(N_0) -kk );
    }
    return k;
}


double Dt_AS(double M) // time sampled from branching					
{
//----------------------------------
//Note on procedure for finding x:
// double PSI = (p-1)*(c**(p1-1))*(1/(T+c)**p)  FIND CDF THEN FIND t (WHICH IS x) t-> INF dt'
//----------------------------------
    
    //double x = (1-p)*pow(RAND()/(p-1),1/(1-p))-c;
	
	double x = c*(pow(RAND(),1/(1-p))-1);
	
	return x;
}

double T_BG() // time sampled from Poissonian background 	TIME FOR EACH BCKGRND EVENT
{
	
	return period*RAND() ;

	}

// generate all the aftershock sequence of a given MS:
double New_Branch(double t, double M)
{
	long N = N_as(t,M); // set number of aftershocks x mainshock
	double dm;
	
	for (int i = 0; i<N; i++)
	{	
	//-----------------------------------------
	// Note on finding m_as:	
		// P(m_as)=b_asln10 10**(-b_as(m_as-Mmin) => u = 10**(-b_as(m_as-Mmin)) => m_as = Mmin-log10(u)/b_as
	//-----------------------------------------	
		//double m_as = Mmin-(b_as*log(10)*pow(10, b_as*Mmin))*(log10(RAND()));

		//double m_as = Mmin+((1.0e0/b_as)*log10(RAND()));	
		
		// Using Defintion from
		// J. Davidsen and M. Baiesi 2015, "Self-similar Aftershock Rates"(Supplemental): 	
		//double m_as = Mmin-(pow(b_as,-1.0e0)*log10(RAND()));		// magnitude of the triggered event
		
		double m_as = Mmin-((1.0e0/b_as)*log10(RAND()));	
		
		double t_as = t + Dt_AS(M-m_as); // NOTE: dt depends on dm
		
		c    = c_0*pow(10.0e0, g*(M-m_as));
        tau  = tau_0*pow(10.0e0, -z*(M-m_as));

		dm = std::abs(M-m_as);
		
		//creating list of dm's :
		list_dM.push_back(dm);
		
		double Mmax = 9.0e0;

		if(t_as < period && m_as < Mmax)
	
		{
			New_Aval(t_as, m_as); //GENERATE NEW SOURCE OF AVALANCHE	 
		}
	}

	return dm;
}


void New_Aval(double t, double M)
{
    
    avals.push_back( make_pair (t,M));
   return;
}

double RAND()
{//  This is an adapted version of subroutine RANECU written by F. James
//  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
//  give a single random number at each call.
    int k,iz;

    k=iseed1/53668;
    iseed1=40014*(iseed1-k*53668)-k*12211;
    if (iseed1<0) iseed1=iseed1+2147483563;

    k=iseed2/52774;
    iseed2=40692*(iseed2-k*52774)-k*3791;
    if (iseed2<0) iseed2=iseed2+2147483399;

    iz=iseed1-iseed2;
    if(iz<1) iz=iz+2147483562;

    return(iz*4.656613e-10);
};

bool compare(const pair<double,double> first, const pair<double,double> second) // decides if an event is after another 
{
		return first.first < second.first;
}
