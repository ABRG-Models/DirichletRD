#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;
using namespace morph;
using namespace std;

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    string logpath = argv[1];
	int scale = stoi(argv[2]);
	double xspan = stod(argv[3]);
    bool LfixedSeed = stoi(argv[4]);
    bool Lgraphics = stoi(argv[5]);
    cout << "after reading input parameters" << endl;
    /*
     * open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
     */
/*
    double dt = conf.getDouble("dt",0.0001);
    double Dn = conf.getDouble("Dn",5.0);
    double Dchi = conf.getDouble("DChi",5.0);
    double Dc = conf.getDouble("Dc",1.5);
    int scale = conf.getInt("scale",8);
    double xspan = conf.getDouble("xspan",5.0);
    int numsteps = conf.getInt("numsteps",100);
    int numAdjust = conf.getInt("numAdjust",1000000);
    int numprint = conf.getInt("numprint",95);
    string logpath = conf.getString("logpath","../logs");
    //string svgpath = conf.getString("svgpath","./rat.svg");
    double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
    double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
    int numSectors = conf.getInt("numsectors",12);
    bool Lcontinue = conf.getBool("Lcontinue",false);
*/
    bool overwrite_logs = true;

    unsigned int seed;
    if (LfixedSeed) {
        seed = 1;
    }
    else {
        seed = time(NULL);
    }

    // A ra2yyndo2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);

    /*
     * NOW create a log directory if necessary, and exit on any
     * failures.
     */
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            cerr << "Failed to create the logpath directory "
                 << logpath << " which does not exist."<< endl;
            return 1;
        }
	}


// initialise DRegion class setting scale
    DRegion M(scale,xspan,logpath);
    cout << "before dissect_boundary " << endl;
    vector<std::pair<double,double>> cGravity;
    cGravity = M.dissectBoundary(); //dissect region boundary
    cout << "Edges size = " << M.edges.size() << endl;
	M.setRadialSegments(); //set the radial segments for regions
    cout << "after first setRadialSegments " << endl;


    double rhoInit = 2.5;
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
    array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
    array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
    array<float,3> cl_d = morph::Tools::getJetColorF (0.00);
    array<float,3> offset = {{0, 0, 0}};
    float hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;// section for solving on the curved boundaries
     int boundaryCount = 0;
     if (Lgraphics) {
        morph::Gdisplay isp(900, 900, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
        cout << "after setting display" << endl;
        isp.resetDisplay (fix, eye, rot);
        cout << "after setting display" << endl;
    //plot stuff here.
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                isp.drawHex (h.position(), (h.d/2.0f), cl_a);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                isp.drawHex (h.position(), (h.d/2.0f), cl_c);
            }
        }
        for (int regcnt = 0; regcnt < NUMPOINTS;regcnt++) {
            std::array<float,3> pos = {M.centres[regcnt].first, M.centres[regcnt].second, 0};
            cout << " drawing centre of region x " << M.centres[regcnt].first << " y " << M.centres[regcnt].second << endl;
            isp.drawHex (pos,offset,1.0*hexWidth,cl_a);
            for (unsigned int i=0; i< M.vCoords[regcnt].size();i++) {
                std::array<float,3> pos = {M.vCoords[regcnt][i].first, M.vCoords[regcnt][i].second, 0};
                isp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }
        }
        cout << "boundaryCount "<<boundaryCount<<endl;
        usleep (10000000); // ten seconds
        isp.redrawDisplay();
        usleep (10000000); // ten seconds
        isp.saveImage(logpath + "/Tesselation.png");
        isp.closeDisplay();
     }
    /*
     * Now we create a vector of ksSolvers
     * each region with a separate ksSolver
     */
    vector<ksSolver> S;
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    M.populateBoundPolygon(1);
    cout << "just after setting polygonal boundaries " << M.curvedBoundary.size()<<endl;
    for (int j = 0;j<NUMPOINTS;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector "<< j <<endl;
    }


    /*
     * now draw the intial tesselation
     */
    int internalCount = 0;
    boundaryCount = 0;
    if (Lgraphics) {
        morph::Gdisplay idisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        idisp.resetDisplay (fix, eye, rot);
        for (int j=0;j<NUMPOINTS;j++) {
        // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
// cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                    idisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                    boundaryCount++;
                }
                else
                {
                    cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
//cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                //idisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                idisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                idisp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }

        } // end of loop over regions
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                idisp.drawHex (h.position(), (h.d/2.0f), cl_b);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                idisp.drawHex (h.position(), (h.d/2.0f), cl_c);
            }
        }
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        idisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        idisp.saveImage(logpath + "/Tesselation1.png");
        idisp.closeDisplay();
    } //end of graphics section

    /*
     * write out the HexGrids
     */
    cout << " write out the HexGrids" << endl;
    for (int j=0; j<NUMPOINTS;j++){
        string str = to_string(j);
        string ename = logpath + "/" + str + ".h5";
        S[j].Hgrid->save(ename);
    }
    return 0;
};
