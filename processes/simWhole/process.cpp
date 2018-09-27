#include "morph/world.h"
#include "morph/sockserve.h"
#include "morph/display.h"
#include "morph/tools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp> 
#define PI 3.14159265

using namespace std;

class Erm2009
{
public:
    int scale;
    int offset;
    int n;
    double ds;
    struct extremum {
    int radialIndex;
    double radialValue;
     };

    vector<vector<double> > X; //
    vector<vector<double> > H; // hex-grid info
    vector<vector<double> > N; // hex neighbourhood
    vector<int> C;
    vector<extremum> turnVal;

    vector<double> NN, CC;
  //class constructor
  Erm2009 (int scale, int offset, double z) {
 
    

        this->scale = scale;
        this->offset = offset;

        H.resize(6);
        double s = pow(2.0, scale-1); //make S bigger to incorporate the circle.
	ds = 29.0/s;
        n = 0;
	double x = 0;
	double y = 0;
	double radius = 0;
	double theta = 0;
	for (int r=-s/2; r<0; r++) {
	  for (int b=-s/2; b<=s; b++) {
	 	      x = (r-b)*sqrt(3.)/(2.*s);
	 	      y = (1.*(r+b))/(2.*s);
		      radius = 29.0*sqrt(x*x+y*y);
		       if (radius == 0)
			theta = 0;
		      else
		        if (y >= 0) 
			  theta = atan2(y,x);  
		        else
			  theta = 2*PI + atan2(y,x);
		      
	  
		      if (x*x + y*y <= 0.75) { // circular
                                 //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                 //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                 //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                 //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                                 //if (x*x < 0.4 && y*y < 0.4) { // square 1
                                 //if (x*x < 0.2 && y*y < 0.2) { // square 2
                                 //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                              H[0].push_back(x);                      //X
                              H[1].push_back(y);                      //Y
                              H[2].push_back(r);                      //r
                              H[3].push_back(b);                      //b
                              H[4].push_back(radius);                 //r
                              H[5].push_back(theta);                  //theta
                              n++;
                          } //if on shape determination
	  } //if on b
	} //if on r
	 

	for (int r=0; r<=s; r++) {
	  for (int b=-s; b<=s; b++) {
		      x = (r-b)*sqrt(3.)/(2.*s);
	 	      y = (1.*(r+b))/(2.*s);
		      radius = 29.0*sqrt(x*x+y*y);
		      if (radius == 0)
			theta = 0;
		      else
		        if (y >= 0) 
			  theta = atan2(y,x);  
		        else
			  theta = 2*PI + atan2(y,x);
		      if (x*x + y*y <= 0.75) { // circular
                                  //if (x*x/0.1 + y*y/0.8 <= 1) { // ellipse 1
                                  //if (x*x/0.8 + y*y/0.1 <= 1) { // ellipse 2
                                  //if (x*x/0.8 + y*y/0.01 <= 0.9) { // ellipse 3
                                  //if (x*x/0.01 + y*y/0.8 <= 0.9) { // ellipse 4
                                  //if (x*x < 0.4 && y*y < 0.4) { // square 1
                                  //if (x*x < 0.2 && y*y < 0.2) { // square 2
                                  //if (x*x/0.05 - y*y/0.05 <= 0.7) { // hyperbolic 1
                               H[0].push_back(x);                      //X
                               H[1].push_back(y);                      //Y
                               H[2].push_back(r);                      //r
                               H[3].push_back(b);                      //b
                               H[4].push_back(radius);                      //r
                               H[5].push_back(theta);                  //theta
                               n++;
                           } //if on shape determination
	   	} //if on b
	    } //if on r
	 
	
       

        

        // get neighbours
        N.resize(n);
        C.resize(n,0);
        for(int i=0;i<n;i++){
            N[i].resize(6,i); // CONNECT ALL TO 'boundary' UNIT AT N+1
            for(int j=0;j<n;j++){
                int dr = H[2][j]-H[2][i];
		// int dg = H[3][j]-H[3][i];
                int db = H[3][j]-H[3][i];
		// if(max(max(abs(dr),abs(dg)),abs(db))==1){
		   if(max(abs(dr),abs(db))==1){
                    //anticlockwise from east
                      if(db==0&&dr==+1){N[i][0]=j;C[i]++;}
                      if(db==+1&&dr== +1){N[i][1]=j;C[i]++;}
                      if(db== +1&&dr==0){N[i][2]=j;C[i]++;}
                      if(db==0&&dr==-1){N[i][3]=j;C[i]++;}
                      if(db==-1&&dr==-1){N[i][4]=j;C[i]++;}
                      if(db== -1&&dr==0){N[i][5]=j;C[i]++;}
		      
                     //previous definition
		     // if(dg==0&&dr==+1){N[i][0]=j;C[i]++;}
                     // if(dg==1&&dr== 0){N[i][1]=j;C[i]++;}
                     // if(dg==0&&dr==-1){N[i][2]=j;C[i]++;}
                     // if(dg==-1&&dr==-1){N[i][3]=j;C[i]++;}
                     // if(dg==-1&&dr== 0){N[i][4]=j;C[i]++;}
                     // if(dg==-1&&dr==+1){N[i][5]=j;C[i]++;}
		    }
            }
        }
        X.resize(n);
        for(int i=0;i<n;i++){
            X[i].resize(3,0.);
            X[i][0] = H[0][i];
            X[i][1] = H[1][i];
        }
        NN.resize(n);
        CC.resize(n);
  } //end of constructor

    vector<double> getLaplacian(vector<double> Q, double dx) {
      double overdxSquare = 1./(1.5*dx*dx);
        vector<double> L(n,0.);
        for(int i=0;i<n;i++){
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overdxSquare;
        }
        return L;
    }

    void step(double dt, double Dn, double Dc) {

       
       // Set up boundary conditions with ghost points
       
        for(int i=0;i<n;i++){
            if(C[i]<6){
                for(int j=0;j<6;j++){
		  if(N[i][j] == i) {
		    // int temp = (j+3) % 6;
		    //   NN[N[i][j]] = NN[N[i][temp]];
		    //   CC[N[i][j]] = CC[N[i][temp]];
		    // I think this is the correct ghost point to enforce no-flux
		     NN[N[i][j]] = NN[i];
		     CC[N[i][j]] = CC[i];
		  }
                }
	    }
	}
    
      

        double beta = 5.;
        double a = 1., b = 1., mu = 1., chi = Dn;
        vector<double> lapN = getLaplacian(NN,ds);
        vector<double> lapC = getLaplacian(CC,ds);

       
        vector<double> G(n,0.);

        for (int i=0;i<n;i++) {
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
	  double dr0N = (NN[N[i][0]]+NN[i])/2.;
	  double dg0N = (NN[N[i][1]]+NN[i])/2.;
	  double db0N = (NN[N[i][2]]+NN[i])/2.;
	  double dr1N = (NN[N[i][3]]+NN[i])/2.;
	  double dg1N = (NN[N[i][4]]+NN[i])/2.;
	  double db1N = (NN[N[i][5]]+NN[i])/2.;
	  //double ncentre = NN[i];

          double dr0C = CC[N[i][0]]-CC[i];
          double dg0C = CC[N[i][1]]-CC[i];
          double db0C = CC[N[i][2]]-CC[i];
	  double dr1C = CC[N[i][3]]-CC[i];
          double dg1C = CC[N[i][4]]-CC[i];
          double db1C = CC[N[i][5]]-CC[i];
	  

	  //finite volume for NdelC, h = s/2
	  G[i] = (dr0N*dr0C+dg0N*dg0C+db0N*db0C+dr1N*dr1C+dg1N*dg1C+db1N*db1C)/(1.5*ds*ds);
	  
	  //G[i] = (drN*drC+dgN*dgC+dbN*dbC)/(6.*ds*ds) + NN[i]*lapC[i];
	 
        } //matches for on i

        // step N
        for(int i=0;i<n;i++){
            NN[i]+=dt*( a-b*NN[i] + Dn*lapN[i] - chi*G[i]);
        }

        // step C
        double N2;
        for(int i=0;i<n;i++){
            N2 = NN[i]*NN[i];
            CC[i]+=dt*( beta*N2/(1.+N2) - mu*CC[i] + Dc*lapC[i] );
        }
    }//end step

  

  //function find_max to find turning points both values and indices.

  int find_maxR(vector<double> ray) {
    
    int iend = ray.size();
    //int iend = 0;
    turnVal.resize(1000);
    cout <<"ray size ="<<iend<<iend<<flush;
    double old_slope = 0;
    double new_slope = 0;
    int count = 0;
    old_slope = ray[1] - ray[0];
    for (int i =2; i<iend;i++){
      new_slope = ray[i]-ray[i-1];
      if (new_slope*old_slope < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = ray[i]; //should really interpolate
        count++;

      old_slope = new_slope;
      }
    }
    return count;
  }

  int find_maxTh(vector<double> ray) {
    ofstream cfile ( "logsWhole/angleTurn.txt" );
    int iend = ray.size();
    turnVal.resize(1000);
    cout <<"ray size ="<<iend<<iend<<flush;
    double old_slope = 0;
    double new_slope = 0;
    int count = 0;
    old_slope = (ray[1] - ray[0]) / fabs(ray[1] -ray[0]);
    for (int i =2; i<=iend+1;i++){
      
      new_slope = (ray[i%iend]-ray[(i-1)%iend]) / fabs((ray[i%iend]-ray[(i-1)%iend]));
      if (new_slope*old_slope < 0.000001) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = ray[i%iend]; //should really interpolate
 cfile << " index "<< turnVal[count].radialIndex << " value " << turnVal[count].radialValue <<endl;
       cfile << "new slope " << new_slope << " old slope" << old_slope<<endl; 
       count++;
      }
	old_slope = new_slope;
    }
      
    return count;
  }
					  
    
//sectorize over radius
  vector <double> sectorize_radius (int numSectors) {
    ofstream dfile ( "logsWhole/sectorRadius.txt" );
    double maxRadius = sqrt(0.75)*29.0;
    vector <double>  radiusCC;
    radiusCC.resize(numSectors);
    vector <int> radiusCount;
    radiusCount.resize(numSectors);
    double startRadius, endRadius, radiusInc; //sector radii
    radiusInc = maxRadius /(1.*numSectors);
    double startAngle, endAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    startAngle = 6*angleInc;
    endAngle = 7*angleInc;
    
    for (int k=0;k<numSectors;k++) {
      //startAngle = fmod((k*angleInc),2*PI);
      //endAngle = fmod(((k+1)*angleInc),2*PI);
      startRadius = (k*radiusInc);
      endRadius = (k+1)*radiusInc;
      
      for (int i=0; i< this->n;i++) {
	if (this->H[5][i]>=startAngle && this->H[5][i]<endAngle) {
	  if (this->H[4][i]>=startRadius && this->H[4][i]<endRadius) {
	      radiusCount[k]++;
	      radiusCC[k] += this->CC[i];
	  } //end of if on radius
	} //end of if on angleSector
      } //end of loop over i
     
      dfile <<" radiusInc "<<radiusInc <<"  startRadius "<<startRadius<<"  endRadius "<<endRadius<<endl;
      if (radiusCount[k] != 0)
	radiusCC[k] = radiusCC[k] / (1.*radiusCount[k]) - 2.5;
      else
	radiusCC[k] += -999.999;
    }//end loop on k
     return radiusCC;
    
  } //end of function sectorize_radius

    // sectorize over angle
    vector <double> sectorize_angle (int numSectors) {
    ofstream efile ( "logsWhole/sectorAngle.txt" );
    double maxRadius = sqrt(0.75)*29.0;
    vector <double>  angleCC;
    angleCC.resize(numSectors);
    vector <int> angleCount;
    angleCount.resize(numSectors);
    double startRadius, endRadius, radiusInc; //sector radii
    radiusInc = maxRadius /(1.*numSectors);
    startRadius = (numSectors-1)*radiusInc;
    endRadius = numSectors*radiusInc;
      
    double startAngle, endAngle, angleInc; //sector angles
    angleInc = 2*PI/(1.*numSectors);
    for (int k=0;k<numSectors;k++) {
      //startAngle = fmod((k*angleInc),2*PI);
      //endAngle = fmod(((k+1)*angleInc),2*PI);
      startAngle = (k*angleInc);
      endAngle = (k+1)*angleInc;
      
      for (int i=0; i< this->n;i++) {
	if (this->H[4][i]>=startRadius && this->H[4][i]<endRadius) {
	if (this->H[5][i]>=startAngle && this->H[5][i]<endAngle) {
	  angleCount[k]++;
	  angleCC[k] += this->CC[i];
	} //end of if on angle
	} //end of if on radiusSector
      } // end of loop over i
      efile <<" angleinc "<<angleInc<< "  start angle "<< startAngle<<"  end angle "<<endAngle <<endl;
      if (angleCount[k] != 0)
	angleCC[k] = angleCC[k] / (1.*angleCount[k]) - 2.5;
      else
	angleCC[k] += -999.999;
    }//end loop on k
     return angleCC;
    
  } //end of function sectorize_angle

  
  //function bessel_ray for computing bessel functions along a ray
  // takes a vector of doubles representing the radius
  // returns the Bessel function values
  // vector<double> bessel_ray (int v, vector<double> ray) {
  //   vector <double> result(n,0.);
  //   result = cyl_bessel_j(v, ray);
  //   return result;
  // }
    
       
}; // Erm2009


int main (int argc, char **argv)
{
    srand(atoi(argv[3]));

    double dt = .0001;            // integration timestep // SHOULD BE 0.001

    // INITIALIZATION
    morph::World W(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),dt);
    //W.waitForConnection();

    // DISPLAYS
    vector<morph::Gdisplay> displays;
    vector<double>fix(3,0.);
    vector<double>eye(3,0.);
    eye[2] = -0.4;
    vector<double>rot(3,0.);
    displays.push_back (morph::Gdisplay (600, "NN field", 0., 0., 0.));
    displays[0].resetDisplay(fix,eye,rot);
    displays[0].redrawDisplay();
    displays.push_back (morph::Gdisplay (600, "CC field", 0., 0., 0.));
    displays[1].resetDisplay(fix,eye,rot);
    displays[1].redrawDisplay();


    Erm2009 M(8,0,0.);
   double  value;
   double factor = 7.0 / 29.0;
  vector <double> ray;
   for (int i=0;i<M.n;i++) {
     value  = boost::math::cyl_bessel_j(0,M.H[4][i]*factor);
       ray.push_back(value);
  }

 

   
   ofstream bfile ("logsWhole/sectorRadius.out");
   ofstream afile ("logsWhole/sectorAngle.out");
   
 
   for (int i=0;i<M.n;i++) {
  //    M.NN[i]=(morph::Tools::randDouble())*0.01 +1.;
    //  M.CC[i]=(morph::Tools::randDouble())*0.01+2.5;
     // M.CC[i]=sin(i*3.142)*0.1+2.5;
     //	 M.CC[i]=sin(i*3.142)*0.1+2.5;
     }

      for (int i=0;i<M.n;i++) {
     //   // M.NN[i]= 1. + sin(i*3.142)*sin(i*3.142)*0.001;
        M.NN[i]= 1. + (morph::Tools::randDouble())*0.1;
        M.CC[i]=ray[i]*0.05+2.5;
       }
   
    unsigned int frameN = 0;
    unsigned int frameM = 0;
    //initial values for Dn Dc
    double Dn = 100.;
    double Dc = Dn*0.3;
    double TIME=0.;
    vector<double*> f;
    f.push_back(&TIME);
    
    int vdegree = 0;
     int mdegree = 0;
     vector <double> rayv;
     vector <double> ringv;
     // double rad; //radius for ringv
     //find the radius for the ring
     // for (int i=0;i<M.n;i++)
     //   if ((M.H[3][i] == 0) && (M.H[2][i] == 50)) //boundary is at 55
     //     rad = M.H[4][i];
  
   
 
  
  

    vector<vector<double*> > S;
    {
        vector<double*> s;
        for(unsigned int i=0;i<M.CC.size();i++){
            s.push_back(&M.CC[i]);
        } S.push_back(s);
    }
    {
        vector<double*> s;
        for(unsigned int i=0;i<M.NN.size();i++){
            s.push_back(&M.NN[i]);
        } S.push_back(s);
    }

    bool doing = true;
    while (doing) {

        std::stringstream TIMEss;
        TIMEss<<setw(10)<<setfill('0')<<TIME;
        const char* TIMEcs = TIMEss.str().c_str();

        std::stringstream out;
        out.clear();
        out.setf(ios::fixed,ios::floatfield);

        // DEFINE OUTPUT MESSAGE
        out<<TIME<<",";

        vector <string> command;
        string messageI=W.master.exchange(out.str().c_str());
        stringstream ss(messageI);
        while (ss.good()){
            string substr;
            getline(ss,substr,',');
            command.push_back(substr);
        } ss.clear();

        // Interpret commands:
        switch (stoi(command[0])) {

        case 0: // *** QUIT ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 0=QUIT"<<endl<<flush;

            W.logfile.close();
            for(unsigned int i=0;i<displays.size();i++){
                displays[i].closeDisplay();
            }
            W.master.closeSocket();
            for(unsigned int i=0;i<W.ports.size();i++){
                W.ports[i].closeSocket();
            }
            doing=false;
            break;
        }

        case 1: // *** STEP ***
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 1=STEP"<<endl<<flush;

            // Exchange comms
            vector< vector <string> > commands(W.ports.size());
            for(unsigned int i=0;i<W.ports.size();i++){
                string messageI=W.ports[i].exchange(out.str().c_str());
                stringstream ss(messageI);
                while (ss.good()){
                    string substr;
                    getline(ss,substr,',');
                    commands[i].push_back(substr);
                } ss.clear();
            }

            // DO STUFF HERE

            M.step(dt, Dn, Dc);
	   
  
            // END STEP
            TIME++;
            break;
        }

        case 2: // *** PLOT ***
        {
	  //W.logfile<<W.processName<<"@"<<TIMEcs<<": 2=PLOT"<<endl<<flush;
            displays[0].resetDisplay(fix,eye,rot);
	    displays[1].resetDisplay(fix,eye,rot);

            vector<double> plt = M.NN;
            double maxV = -1e7;
            double minV = +1e7;
            for (int i=0;i<M.n;i++) {
                if (M.C[i]==6) {
                    if (plt[i]>maxV) { maxV = plt[i]; }
                    if (plt[i]<minV) { minV = plt[i]; }
                }
            }
            double scaleV = 1./(maxV-minV);
            vector<double> P(M.n,0.);
            for (int i=0;i<M.n;i++) {
                P[i] = fmin (fmax (((plt[i])-minV)*scaleV,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<M.n;i++) {
                vector <double> cl = morph::Tools::getJetColor(P[i]);
                displays[0].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl);
                displays[0].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl);
            }
            displays[0].redrawDisplay();
	    
	    vector<double> plu = M.CC;
            double maxV1 = -1e7;
            double minV1 = +1e7;
            for (int i=0;i<M.n;i++) {
                if (M.C[i]==6) {
                    if (plu[i]>maxV1) { maxV1 = plu[i]; }
                    if (plu[i]<minV1) { minV1 = plu[i]; }
                }
            }
            double scaleV1= 1./(maxV1-minV1);
            vector<double> P1(M.n,0.);
            for (int i=0;i<M.n;i++) {
                P1[i] = fmin (fmax (((plu[i])-minV1)*scaleV1,0.),1.);
                // M.X[i][2] = P[i];
            }

            for (int i=0;i<M.n;i++) {
                vector <double> cl1 = morph::Tools::getJetColor(P1[i]);
                displays[1].drawTriFill(M.X[i],M.X[M.N[i][0]],M.X[M.N[i][1]],cl1);
                displays[1].drawTriFill(M.X[i],M.X[M.N[i][3]],M.X[M.N[i][4]],cl1);
            }
            displays[1].redrawDisplay();
            break;
        }

        case 3:
        {
	  
            std::stringstream frameFile1;
            frameFile1<<"logsWhole/"<<W.processName<<"frameNN";
            frameFile1<<setw(5)<<setfill('0')<<frameN;
            frameFile1<<".png";
            displays[0].saveImage(frameFile1.str());
            frameN++;
	    std::stringstream frameFile2;
            frameFile2<<"logsWhole/"<<W.processName<<"frameCC";
            frameFile2<<setw(5)<<setfill('0')<<frameM;
            frameFile2<<".png";
            displays[1].saveImage(frameFile2.str());
            frameM++;

	     
	    


            break;
        }

        case 4:
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 4=RECD"<<endl<<flush;
            switch(stoi(command[1]))
            {
            case 0: Dn = stod(command[2]); break;
            case 1: Dc = stod(command[2]); break;
            }
            break;
        }

        case 5: // *** SAVE ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 5=SAVE"<<endl<<flush;
	    

            if (command.size()==2) {
                ofstream outFile;
                std::stringstream oFile; oFile<<command[1];
                outFile.open(oFile.str().c_str(),ios::out|ios::binary);
                double n=(double)S.size();double* N=&n;
                outFile.write((char*)N,sizeof(double));
                for (unsigned int i=0;i<S.size();i++) {
                    double m=(double)S[i].size();double* M=&m;
                    outFile.write((char*)M,sizeof(double));
                    for (unsigned int j=0;j<S[i].size();j++) {
                        outFile.write((char*)S[i][j],sizeof(double));
                    }
                }
                outFile.close();
            } else {
                W.logfile<<"No output filename."<<endl<<flush;
            }

	    rayv.resize(0);
      
      
	    ringv.resize(0);
	    int numSectors = 12;
	    
            rayv = M.sectorize_radius(numSectors);
	    for (int i=0;i < (int) rayv.size();i++) {
	      bfile << " " << rayv[i] << endl;
	    }
   
            ringv = M.sectorize_angle(numSectors);
	    for (int i=0;i < (int) ringv.size();i++) {
	      afile << " " << ringv[i] << endl;
	    }
            vdegree = M.find_maxR(rayv);

            mdegree = M.find_maxTh(ringv);
	    W.logfile<<W.processName<<"  vdegree= "<<vdegree<<endl<<flush;
            W.logfile<<W.processName<<"  mdegree= "<<mdegree<<endl<<flush;
	    // W.logfile<<W.processName<<"  rad  "<<rad<<endl<<flush;
            // W.logfile<<W.processName<<"  ringvsize= "<<temp<<endl<<flush;
	    // W.logfile<<W.processName<<"  ds "<<dstemp<<M.ds<<endl<<flush;

            break;
        }

        case 6: // *** LOAD ***
        {
            W.logfile<<W.processName<<"@"<<TIMEcs<<": 6=LOAD"<<endl<<flush;
	    W.logfile<<W.processName<<"value of n "<<M.n<<endl<<flush;
	    
            if (command.size()==2) {
                double dummy;
                ifstream inFile;
                std::stringstream iFile;
                iFile<<command[1];
                inFile.open(iFile.str().c_str(),ios::in|ios::binary);
                inFile.read((char*)&dummy,sizeof(double));
                int I=(int)dummy;
                if (I==static_cast<int>(S.size())) {
                    for (int i=0;i<I;i++) {
                        inFile.read((char*)&dummy,sizeof(double));
                        int J=(int)dummy;
                        if (J==static_cast<int>(S[i].size())) {
                            for (int j=0;j<J;j++) {
                                inFile.read((char*)S[i][j],sizeof(double));
                            }
                        } else {
                            W.logfile<<"Wrong dims I."<<endl<<flush;
                        }
                    }
                } else {
                    W.logfile<<"Wrong dims J."<<endl<<flush;
                }
                inFile.close();
            } else {
                W.logfile<<"No input filename."<<endl<<flush;
            }
	
	
	    
            break;
        }
	}
    }

    return 0;
};
