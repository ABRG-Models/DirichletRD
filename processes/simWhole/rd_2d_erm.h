#include "morph/display.h"
#include "morph/tools.h"
#include "morph/ReadCurves.h"
#include "morph/HexGrid.h"
#include "morph/HdfData.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#ifdef __GLN__ // Currently I've only tested OpenMP on Linux
#include <omp.h>
#endif
#include <hdf5.h>
#include <unistd.h>

#define DEBUBG 1
#define SAVE_PNGS 1
#define DBGSTREAM std::cout
#include <morph/MorphDbg.h>

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;

using morph::HexGrid;
using morph::ReadCurves;
using morph::HdfData;
using namespace std;

/*!
 * Reaction diffusion system; Ermentrout 2009.
 */
class RD_2D_Erm
{
public:

    /*!
     * Constants
     */
    //@{
    //! Square root of 3 over 2
    const double R3_OVER_2 = 0.866025403784439;
    //! Square root of 3
    const double ROOT3 = 1.73205080756888;
    //@}

    /*!
     * The logpath for this model. Used when saving data out.
     */
    string logpath = "logs";

    /*!
     * Setter which attempts to ensure the path exists.
     */
    void setLogpath (const string p) {
        this->logpath = p;
        // Ensure log directory exists
        morph::Tools::createDir (this->logpath);
    }

    /*!
     * Frame number, used when saving PNG movie frames.
     */
    unsigned int frameN = 0;

    /*!
     * Holds the number of hexes in the populated HexGrid
     */
    unsigned int nhex = 0;

    /*!
     * Set N>1 for maintaing multiple expression gradients
     */
    unsigned int N = 1;

    /*!
     * The c_i(x,t) variables from the Ermentrout paper (chemoattractant concentration)
     */
    vector<vector<double> > c;

    /*!
     * The n_i(x,t) variables from the Ermentrout paper (density of tc axons)
     */
    vector<vector<double> > n;

    /*!
     * Holds the Laplacian
     */
    vector<vector<double> > lapl;

    /*!
     * Holds the Poisson terms (final non-linear term in Ermentrout equation 1)
     */
    vector<vector<double> > poiss;

    /*!
     * Our choice of dt.
     */
    double dt = 0.0001;

    /*!
     * Compute half and sixth dt in constructor.
     */
    //@{
    double halfdt = 0.0;
    double sixthdt = 0.0;
    //@}

    /*!
     * The HexGrid "background" for the Reaction Diffusion system.
     */
    HexGrid* hg;

    /*!
     * Store Hex positions for saving.
     */
    //@{
    vector<float> hgvx;
    vector<float> hgvy;
    //@}

    /*!
     * Hex to hex distance. Populate this from hg.d after hg has been
     * initialised.
     */
    double d = 1.0;
    /*!
     * radius of circular boundary, set from class constructor
     */
    double radius = 0.65;
    

    /*!
     * Parameters of the Ermentrout model - default values.
     */
    //@{
    //! Diffusion constant for n
    double Dn = 0.3;
    //! Diffusion constant for c
    double Dc = Dn * 0.3;
    //! saturation term in function for production of c
    double beta = 5.0;
    //! production of new axon branches
    double a = 1.0;
    //! pruning constant
    double b = 1.0;
    //! decay of chemoattractant constant
    double mu = 1.0;
    //! degree of attraction of chemoattractant
    double chi = Dn;
    //@}

    /*!
     * Track the number of computational steps that we've carried
     * out. Only to show a message saying "100 steps done...", but
     * that's reason enough.
     */
    unsigned int stepCount = 0;

    /*!
     * Simple constructor; no arguments.
     */
    RD_2D_Erm (void) {
        this->halfdt = this->dt/2.0;
        this->sixthdt = this->dt/6.0;
    }

     /*!
     * Constructor with radius argument
     */
    RD_2D_Erm (double rad) {
        this->radius = rad;
        this->halfdt = this->dt/2.0;
        this->sixthdt = this->dt/6.0;
    }

    /*!
     * Destructor required to free up HexGrid memory
     */
    ~RD_2D_Erm (void) {
        delete (this->hg);
    }

    /*!
     * A utility function to resize the vector-vectors that hold a
     * variable for the N different thalamo-cortical axon types.
     */
    void resize_vector_vector (vector<vector<double> >& vv) {
        vv.resize (this->N);
        for (unsigned int i =0; i<this->N; ++i) {
            vv[i].resize (this->nhex, 0.0);
        }
    }

    /*!
     * Resize a variable that'll be nhex elements long
     */
    void resize_vector_variable (vector<double>& v) {
        v.resize (this->nhex, 0.0);
    }

    /*!
     * Resize a parameter that'll be N elements long
     */
    void resize_vector_param (vector<double>& p) {
        p.resize (this->N, 0.0);
    }

    /*!
     * Initialise this vector of vectors with noise.
     */
    void noiseify_vector_vector (vector<vector<double> >& vv, double off, double sig) {
        for (unsigned int i = 0; i<this->N; ++i) {
            for (auto h : this->hg->hexen) {
                vv[i][h.vi] = morph::Tools::randDouble() *sig + off;
            }
        }
    }

    /*!
     * Initialise HexGrid, variables. Carry out any one-time
     * computations of the model.
     */
    void init (vector<morph::Gdisplay>& displays, bool useSavedGenetics = false) {

        DBG ("called");
        // Create a HexGrid
        this->hg = new HexGrid (0.01, 3, 0, morph::HexDomainShape::Boundary);
        // Read the curves which make a boundary
	// ReadCurves r("./trial.svg");
	CirclePath r(this->radius);
	
        // Set the boundary in the HexGrid
        this->hg->setBoundary (r);
        // Compute the distances from the boundary
        this->hg->computeDistanceToBoundary();
        // Vector size comes from number of Hexes in the HexGrid
        this->nhex = this->hg->num();
        // Spatial d comes from the HexGrid, too.
        this->d = this->hg->getd();
        // Save hex positions in vectors for datafile saving
        for (auto h : this->hg->hexen) {
            this->hgvx.push_back (h.x);
            this->hgvy.push_back (h.y);
        }

        // Resize and zero-initialise the various containers
        this->resize_vector_vector (this->c);
        this->resize_vector_vector (this->n);
        this->resize_vector_vector (this->lapl);
        this->resize_vector_vector (this->poiss);

        // Initialise a with noise
        this->noiseify_vector_vector (this->n, 1., 0.1);
        this->noiseify_vector_vector (this->c, beta*0.5, 0.1);
	DBG ("value of n(100)" << n[100] << c[100]);
    }

    /*!
     * Computations
     */
    //@{

    /*!
     * Compute one step of the model
     */
    void step (void) {

        this->stepCount++;

        if (this->stepCount % 100 == 0) {
            DBG ("System computed " << this->stepCount << " times so far...");
        }

        for (unsigned int i=0; i<this->N; ++i) {

            this->compute_poiss (n[i],c[i],i);  // compute the non-linear Poission term in Eq1
            this->compute_lapl (n[i], i);       // populate lapl[i] with laplacian of n

            // integrate n
            for (unsigned int h=0; h<this->nhex; ++h) {
                n[i][h] += (a - b*n[i][h] + Dn*lapl[i][h] - chi*poiss[i][h])*dt;
            }

            this->compute_lapl (c[i], i);       // populate lapl[i] with laplacian of c

            // integrate c
            double n2;
            for (unsigned int h=0; h<this->nhex; ++h) {
                n2 = n[i][h]*n[i][h];
                c[i][h] += (beta*n2/(1.+n2) - mu*c[i][h] +Dc*lapl[i][h])*dt;
            }
        }
    }

    /*!
     * Computes the Laplacian
     * Stable with dt = 0.0001;
     */
    void compute_lapl (vector<double>& fa, unsigned int i) {

        double norm  = (2) / (3 * this->d * this->d);

#pragma omp parallel for schedule(dynamic,50)
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            Hex* h = this->hg->vhexen[hi];
            // 1. The D Del^2 term

            // Compute the sum around the neighbours
            double thesum = -6 * fa[h->vi];
            if (h->has_ne) {
                thesum += fa[h->ne->vi];
            } else {
                thesum += fa[h->vi]; // A ghost neighbour-east with same value as Hex_0
            }
            if (h->has_nne) {
                thesum += fa[h->nne->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nnw) {
                thesum += fa[h->nnw->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nw) {
                thesum += fa[h->nw->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nsw) {
                thesum += fa[h->nsw->vi];
            } else {
                thesum += fa[h->vi];
            }
            if (h->has_nse) {
                thesum += fa[h->nse->vi];
            } else {
                thesum += fa[h->vi];
            }

            this->lapl[i][h->vi] = norm * thesum;
        }
    }

    /*!
     * Computes the Poisson term
     *
     * Stable with dt = 0.0001;
     */
    void compute_poiss (vector<double>& fa1, vector<double>& fa2, unsigned int i) {

        // Compute non-linear term

#pragma omp parallel for schedule(dynamic,50)
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            Hex* h = this->hg->vhexen[hi];

            vector<double> dum1(6,fa1[h->vi]);
            vector<double> dum2(6,fa2[h->vi]);

            if (h->has_ne) {
                dum1[0] = fa1[h->ne->vi];
                dum2[0] = fa2[h->ne->vi];
            }
            if (h->has_nne) {
                dum1[1] = fa1[h->nne->vi];
                dum2[1] = fa2[h->nne->vi];
            }
            if (h->has_nnw) {
                dum1[2] = fa1[h->nnw->vi];
                dum2[2] = fa2[h->nnw->vi];
            }
            if (h->has_nw) {
                dum1[3] = fa1[h->nw->vi];
                dum2[3] = fa2[h->nw->vi];
            }
            if (h->has_nsw) {
                dum1[4] = fa1[h->nsw->vi];
                dum2[4] = fa2[h->nsw->vi];
            }
            if (h->has_nse) {
                dum1[5] = fa1[h->nse->vi];
                dum2[5] = fa2[h->nse->vi];
            }

#ifdef SECOND_REPORT_SOLUTION
            // John Brooke's 'second interim report' solution
            double val =
            (dum1[0]+dum1[1])*(dum2[0]-fa2[h->vi])+
            (dum1[1]+dum1[2])*(dum2[1]-fa2[h->vi])+
            (dum1[2]+dum1[3])*(dum2[2]-fa2[h->vi])+
            (dum1[3]+dum1[4])*(dum2[3]-fa2[h->vi])+
            (dum1[4]+dum1[5])*(dum2[4]-fa2[h->vi])+
            (dum1[5]+dum1[0])*(dum2[5]-fa2[h->vi]);
            this->poiss[i][h->vi] = val / (ROOT3 * this->d * this->d);
#endif

            // John Brooke's final thesis solution (based on 'finite volume method'
            // of Lee et al. https://doi.org/10.1080/00207160.2013.864392
            double val =
            (dum1[0]+fa1[h->vi])*(dum2[0]-fa2[h->vi])+
            (dum1[1]+fa1[h->vi])*(dum2[1]-fa2[h->vi])+
            (dum1[2]+fa1[h->vi])*(dum2[2]-fa2[h->vi])+
            (dum1[3]+fa1[h->vi])*(dum2[3]-fa2[h->vi])+
            (dum1[4]+fa1[h->vi])*(dum2[4]-fa2[h->vi])+
            (dum1[5]+fa1[h->vi])*(dum2[5]-fa2[h->vi]);
            this->poiss[i][h->vi] = val / (3 * this->d * this->d);
        }
    }
    //@} // computations

    /*!
     * Analysis functions
     */
    //@{
 
  struct extremum {
    int radialIndex;
    double radialValue;
  };

   vector<extremum> turnVal;
  
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
   vector<double> bessel_ray (int v, vector<double> ray) {
     vector <double> result(n,0.);
     result = cyl_bessel_j(v, ray);
     return result;
   }
    //@} //analysis functions

   
    /*!
     * Plotting functions
     */
    //@{

    /*!
     * Plot the system on @a disps
     */
    void plot (vector<morph::Gdisplay>& disps) {

        this->plot_f (this->n, disps[0], true);
        this->plot_f (this->c, disps[1], true);

	//#ifdef SAVE_PNGS
        std::stringstream frameFile1;
        frameFile1<<"./logs/tmp/demo";
	// frameFile1<<std::setw(5)<<setfill('0')<<frameN;
        frameFile1<<".png";
        disps[0].saveImage(frameFile1.str());
        frameN++;
	//#endifting functions
	//#endif
    }  
    

    /*!
     * Plot a or c
     */
    void plot_f (vector<vector<double> >& f, morph::Gdisplay& disp, bool combinedNorm) {
        vector<double> fix(3, 0.0);
        vector<double> eye(3, 0.0);
        vector<double> rot(3, 0.0);

        vector<double> maxa (this->N, -1e7);
        vector<double> mina (this->N, +1e7);
        // Copies data to plot out of the model
        if (combinedNorm) {
            double maxb = -1e7;
            double minb = +1e7;
            for (auto h : this->hg->hexen) {
                if (h.onBoundary() == false) {
                    for (unsigned int i = 0; i<this->N; ++i) {
                        if (f[i][h.vi]>maxb) { maxb = f[i][h.vi]; }
                        if (f[i][h.vi]<minb) { minb = f[i][h.vi]; }
                    }
                }
            }
            for (unsigned int i = 0; i<this->N; ++i) {
                mina[i] = minb;
                maxa[i] = maxb;
            }
        } else {
            for (auto h : this->hg->hexen) {
                if (h.onBoundary() == false) {
                    for (unsigned int i = 0; i<this->N; ++i) {
                        if (f[i][h.vi]>maxa[i]) { maxa[i] = f[i][h.vi]; }
                        if (f[i][h.vi]<mina[i]) { mina[i] = f[i][h.vi]; }
                    }
                }
            }
        }
        vector<double> scalea (this->N, 0);
        for (unsigned int i = 0; i<this->N; ++i) {
            scalea[i] = 1.0 / (maxa[i]-mina[i]);
        }

        // Determine a colour from min, max and current value
        vector<vector<double> > norm_a;
        this->resize_vector_vector (norm_a);
        for (unsigned int i = 0; i<this->N; ++i) {
            for (unsigned int h=0; h<this->nhex; h++) {
                norm_a[i][h] = fmin (fmax (((f[i][h]) - mina[i]) * scalea[i], 0.0), 1.0);
            }
        }

        // Create an offset which we'll increment by the width of the
        // map, starting from the left-most map (f[0])
        float hgwidth = this->hg->getXmax()-this->hg->getXmin();
        array<float,3> offset = {{0.0f,0.0f,0.0f}};//{ 2*(-hgwidth-(hgwidth/20)), 0.0f, 0.0f };

        // Draw
        disp.resetDisplay (fix, eye, rot);
        for (unsigned int i = 0; i<this->N; ++i) {
            for (auto h : this->hg->hexen) {
                array<float,3> cl_a = morph::Tools::getJetColorF (norm_a[i][h.vi]);
                disp.drawHex (h.position(), offset, (h.d/2.0f), cl_a);
            }
            offset[0] += hgwidth + (hgwidth/20);
        }
        disp.redrawDisplay();
    }
    //@} // plotting

    /*!
     * HDF5 file saving/loading methods
     */
    //@{

    /*!
     * Save positions of the hexes - note using two vector<floats>
     * that have been populated with the positions from the HexGrid,
     * to fit in with the HDF API.
     */
    void saveHexPositions (HdfData& dat) {
        dat.add_float_vector ("/x", this->hgvx);
        dat.add_float_vector ("/y", this->hgvy);
        // And hex to hex distance:
        dat.add_double ("/d", this->d);
    }

    /*!
     * Save some data like this.
     */
    void saveState (void) {
        string fname = this->logpath + "/2Derm.h5";
        HdfData data (fname);

        // Save some variables
        for (unsigned int i = 0; i<this->N; ++i) {

            stringstream vss;
            vss << "c_" << i;
            string vname = vss.str();
            data.add_double_vector (vname.c_str(), this->c[i]);
            vname[0] = 'n';
            data.add_double_vector (vname.c_str(), this->n[i]);
        }

        // Parameters
        data.add_double ("/Dn", this->Dn);
        data.add_double ("/Dc", this->Dc);
        data.add_double ("/beta", this->beta);
        data.add_double ("/a", this->a);
        data.add_double ("/b", this->b);
        data.add_double ("/mu", this->mu);
        data.add_double ("/chi", this->chi);

        // HexGrid information
        this->saveHexPositions (data);
    }
    //@} // HDF5

}; // RD_2D_Erm
