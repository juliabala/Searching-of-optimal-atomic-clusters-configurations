#include <cmath>
#include <vector>
#include <cstdlib>
#include <cstdio>

using namespace std;

inline double angle_between_vectors(double ax,double ay, double az, double bx,double by, double bz);
inline void rotation(double tx,double ty,double tz ,double angle,double & x ,double & y, double & z);
inline void multiply_quaternions(vector<double>& r,const vector<double>& p, const vector<double>& q, const int & itype);
inline void compute_uvec(const vector<double>& rotation,const vector<double>& ub,vector<double>& u);
void test_rotation();
	
	
/**********************************************************************
 * kat miedzy wektorami
 **********************************************************************/		
inline double angle_between_vectors(double ax,double ay, double az, double bx,double by, double bz)
{
	double r1,r2;
	double cs,angle;
	r1=sqrt(pow(ax,2)+pow(ay,2)+pow(az,2));
	r2=sqrt(pow(bx,2)+pow(by,2)+pow(bz,2));
	cs=(ax*bx+ay*by+az*bz)/r1/r2;
	angle=acos(cs);
	return angle;
}	
	
	
	
/*
 ************************   QUATERNION OF ROTATION   **************************************
 * 
 *  wykonujemy obrot wokol osi wyznaczonej przez wektor (tx,ty,tz)  na kwaternionie quater (v1)
 *  korzystamy z kwaternionow: v2 = t*v1*conjugate(t)
 * 
 * ****************************************************************************************
 */

inline void rotation(double tx,double ty,double tz ,double angle,double & x ,double & y, double & z)
{
	vector<double> p(4,0.);
	vector<double> q(4,0.);
	vector<double> quater(4,0.);
	double cs,sn;
	
	cs=cos(angle/2);
	sn=sin(angle/2);
	
	//zamieniamy wektor na kwaternion
	quater[0]=0.;
	quater[1]=x;
	quater[2]=y;
	quater[3]=z;
	
	// ---- normalizacja wektora osi obrotu -----
	double r;
	r=sqrt(tx*tx+ty*ty+tz*tz);
	tx=tx/r;
	ty=ty/r;
	tz=tz/r;
	
	p[0]=cs;
	p[1]=sn*tx;
	p[2]=sn*ty;
	p[3]=sn*tz;
	
	multiply_quaternions(q,p,quater,1);	
	multiply_quaternions(quater,q,p,-1);	
	
	//odczytujemy wspolrzedne obroconego wektora
	x=quater[1];
	y=quater[2];
	z=quater[3];

	return;
}


/*
 ************************   MULTIPLY QUATERNIONS   **************************************
 *  mnozenie kwaternionow:  r=p*q
 * itype=1,-1:   1-zwykly iloczyn, -1-drugi kwaternion sprzezony
 * 
 * **************************************************************************************
 */

inline void multiply_quaternions(vector<double>& r,const vector<double>& p, const vector<double>& q, const int & itype)
{
	
	if(itype!=(-1) && itype!=1){
		printf("bledna wartosc itype = %d\n",itype);
		exit(0);
	}

	r[0]=p[0]*q[0]+(-p[1]*q[1]-p[2]*q[2]-p[3]*q[3])*itype;
	r[1]=p[1]*q[0]+(p[0]*q[1]-p[3]*q[2]+p[2]*q[3])*itype;
	r[2]=p[2]*q[0]+(p[3]*q[1]+p[0]*q[2]-p[1]*q[3])*itype;
	r[3]=p[3]*q[0]+(-p[2]*q[1]+p[1]*q[2]+p[0]*q[3])*itype;
	
	return;
}


/*
 ************************   ACTUAL VECTORS U   ******************************************
 *  liczymy aktualne polozenie wektorow ui
 * ui na wejsciu to kwaterniony bazowe
 * na wyjsciu dostajemy wektor kierunkow
 * 
 * **************************************************************************************
 */
inline void compute_uvec(const vector<double>& rotation,const vector<double>& ub,vector<double>& u)
{
	vector<double> r(4,0.);
	vector<double> p(4,0.);
	multiply_quaternions(p,rotation,ub,1);
	multiply_quaternions(r,p,rotation,-1);
	u[0]=r[1];
	u[1]=r[2];
	u[2]=r[3];
	
	return;
}


void test_rotation()
{
	
	vector<double> rot(4,0.);
	vector<double> rin(4,0.);
	vector<double> rout(3,0.);
	double angle;
	
	angle=4*M_PI/3.;

	
	rot[0]=cos(angle/2);
	rot[1]=sin(angle/2)*1/sqrt(3.);
	rot[2]=sin(angle/2)*1/sqrt(3.);
	rot[3]=sin(angle/2)*1/sqrt(3.);
	
	rin[0]=0;
	rin[1]=0;
	rin[2]=0;
	rin[3]=1;
	
	//compute_uvec(rot,rin,rout);	
	
	printf("%10.3f %10.3f %10.3f\n",rout[0],rout[1],rout[2]);	
}
	
