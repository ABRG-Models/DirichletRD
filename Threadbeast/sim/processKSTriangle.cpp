/*
using morph::HexGrid;
using morph::HdfData;:
using morph::Tools;
using morph::Display
*/
#include "ksRegion.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;
using namespace morph;
using namespace std;
using namespace std::chrono;

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    string jsonfile = argv[1];
    vector<morph::BezCurvePath<float>> triangleBound;
    triangleBound.resize(0);
    hexGeometry* hGeom;
    double dTress = 0.5;
    // open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
    cout << "just before reading parameters" << endl;
    double dt = conf.getDouble("dt",0.0002);
    cout << "dt " << dt << endl;
    double Dn = conf.getDouble("Dn",36.0);
    cout << "Dn " << Dn << endl;
    double Dchi = conf.getDouble("DChi",0.0);
    cout << "Dchi " << Dchi << endl;
    double Dc = conf.getDouble("Dc",12.0);
    cout << "Dc  " << Dc << endl;
    int scale = conf.getInt("scale",7);
    cout << "scale " << scale << endl;
    double xspan = conf.getDouble("xspan",4.0);
    cout << "xspan " << xspan << endl;
    int numsteps = conf.getInt("numsteps",200);
    cout << "numsteps " << numsteps << endl;
    int numAdjust = conf.getInt("numAdjust",101);
    cout << "numAdjust " << numAdjust << endl;
    int numprint = conf.getInt("numprint",85);
    cout << "numprint " << numprint << endl;
    string logpath = conf.getString("logpath","../logs");
    cout <<  logpath << endl;
    bool Lgraphics = conf.getBool("Lgraphics", true);
    cout << "Lgraphics " << Lgraphics << endl;
    bool LDn = conf.getBool("LDn", true);
    cout << "LDn " << LDn << endl;
    double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
    cout << "boundarFalloffDist " << boundaryFalloffDist << endl;
    double aNoiseGain = conf.getDouble("aNoiseGain",0.2);
    cout << "aNoiseGain " << aNoiseGain << endl;
    int numSectors = conf.getInt("numsectors",14);
    cout << "numSectors " << numSectors << endl;
    bool Lcontinue = conf.getBool("Lcontinue",false);
    cout << "Lcontinue " << Lcontinue << endl;
    double nnInitialOffset = conf.getDouble("nnInitialOffset", 1.0);
    cout << "nnInitialOffset " << nnInitialOffset << endl;
    double ccInitialOffset = conf.getDouble("ccInitialOffset",2.5);
    cout << "ccInitialOffset " << ccInitialOffset << endl;
    bool overwrite_logs = conf.getBool("overwrite_logs",true);
    cout << "overwrite_logs " << overwrite_logs << endl;
    bool skipMorph  = conf.getBool("skipMorph",true);
    cout << "skipMorph " << skipMorph << endl;
    bool lPerturb = conf.getBool("lPerturb",true);
    cout << "lPerturb " << lPerturb << endl;
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    cout << "LfixedSeed " << LfixedSeed << endl;
    double ratio = conf.getDouble("ratio", 1.0);
    cout << "ratio  " << ratio  << endl;
    double pRatio = conf.getDouble("pRatio",0.0);
    cout << "pRatio " << pRatio << endl;
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << " pRatio " << pRatio << " ratio " << ratio << "  " << logpath <<  endl;
    ofstream afile (logpath + "/centroids.out",ios::app);
    /*
     * decide if we are using a fixed or
     * run-time dependent seed. Use milliseconds to
     * allow parallel launching
     */

    unsigned int seed;
    chrono::milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    if (LfixedSeed)
        seed = 1;
    else
        seed = static_cast<unsigned int> (ms1.count());

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
     else {
        // Directory DOES exist. See if it contains a previous run and
        // exit without overwriting to avoid confusion.
        if (overwrite_logs == false
            && (morph::Tools::fileExists (logpath + "/params.json") == true
                || morph::Tools::fileExists (logpath + "/guidance.h5") == true
                || morph::Tools::fileExists (logpath + "/positions.h5") == true)) {
            cerr << "Seems like a previous simulation was logged in " << logpath << ".\n"
                 << "Please clean it out manually, choose another directory or set\n"
                 << "overwrite_logs to true in your parameters config JSON file." << endl;
            return 1;
        }
    }
    ofstream gfile ( logpath + "/edges.out");
    ofstream jfile ( logpath + "/results.txt");
	ofstream degfile1 (logpath + "/degree1.data");
	ofstream degfile2 (logpath + "/degree2.data");
	ofstream degfile3 (logpath + "/degree3.data");
    ifstream vfile ("./triangle.inp");
    /*
     * analysis tools
     */
    Analysis L;

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate ksSolver
     */
    vector<ksSolver> S;
    morph::BezCurvePath<float> outer;
    const int rowX = 2; const int rowY = 2;
    //next line for creating a tessellation of equilateral triangles
    //triangleBound = h->eqTriangleTess(dTress, M.centres, M.outerBound);
    //next lines for creating a tessellation of isosceles triangles.
    //if ratio is set to 1.0 we get equilateral triangles
    //if pRatio is > 0.0 and lPerturb == true then we get a tessellation of
    //perturbed triangles.
    vector<vector<hexGeometry::point>> vtxs;
    vtxs.resize(NUMPOINTS);
    vector<hexGeometry::point> vertices;
    vertices = hGeom->isosVertices(ratio, rowX, rowY, pRatio, lPerturb);
    vector<vector<int>> vIndices;
    vector<vector<int>> nbrList =hGeom->triangleNeighbors(rowX, rowY, vIndices);
    triangleBound = hGeom->genTriangleTess(rowX, rowY, vertices, vIndices, outer);
    for (int j=0; j<NUMPOINTS; j++) {
        vtxs[j].push_back(vertices[vIndices[j][0]]);
        vtxs[j].push_back(vertices[vIndices[j][1]]);
        vtxs[j].push_back(vertices[vIndices[j][2]]);
    }
    //line below specific to perfect isosceles triangles
    //triangleBound = hGeom->isosTriangleTess(ratio, rowX, rowY, vtxs, outer);
    cout << "size of triangleBound " << triangleBound.size() << endl;
    /*
     * Tessellation
     */
    Tessellation M(scale, xspan, logpath, outer, vtxs);
    cout << "after Tessellation" << endl;
    /*
     * Solve the KS equations over the regions
     */
    for (int j = 0;j<NUMPOINTS;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, triangleBound[j], M.centres[j]));
        //S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], centroid));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }



	cout << "just after populating the regions from the ksSolver vector"<<endl;
// now draw the intial tesselation
    int internalCount = 0;
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
    if (Lgraphics) {
        morph::Gdisplay idisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        idisp.resetDisplay (fix, eye, rot);
        // plot stuff here.
        for (int j=0; j<NUMPOINTS; j++) {
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    idisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                }
            }
            std::array<float,3> pos = {M.centres[j].first, M.centres[j].second, 0};
            cout << " drawing centre of region x " << M.centres[j].first << " y " << M.centres[j].second << endl;
            idisp.drawHex (pos,offset,1.0*hexWidth,cl_a);
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                 std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                 idisp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }
        }

        idisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        idisp.saveImage(logpath + "/Tesselation1.png");
        idisp.closeDisplay();
    } //end of graphics section
// initialise the fields
    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
        for (unsigned int j=0;j<NUMPOINTS;j++)
        {
		    std::string ccstr = "c" + to_string(j);
		    cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
			char * ccst = new char[ccstr.length()+1];
			std::strcpy(ccst,ccstr.c_str());
		    std::string nstr = "n" + to_string(j);
			char * nst = new char[nstr.length()+1];
			std::strcpy(nst,nstr.c_str());
			cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
	        ginput.read_contained_vals(ccst,S[j].CC);
	        ginput.read_contained_vals(nst,S[j].NN);
		}
	  }
      else
	  {
        for (unsigned int j=0;j<NUMPOINTS;j++)
		{
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        double choice = ruf.get();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    int boundaryCount = 0;
    // begin first time stepping loop
    for (int i=0;i<numsteps;i++)
    {
	    int countHex = 0;
   	    for (int j = 0;j<NUMPOINTS;j++) //loop over regions
        {
     	    S[j].step(dt, Dn, Dchi, Dc);
		    if (i%100 == 0)
		    cout << "just after morph 0 squid i "<< j << " NN[5] " << S[j].NN[5] << endl;
	    }
	    if (i%numprint == 0 && Lgraphics)
	    {
          //set up display
            morph::Gdisplay iidisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            iidisp.resetDisplay (fix, eye, rot);
   	        for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	        {
	 	        cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
		        cout << "in loop over regions " << " size is " << regsize << endl;
			    countHex += regsize;
		        vector<double> regionNN;
			    vector<double> tempNN;
		      //for (auto h : S[j].Hgrid->hexen)
			    for (unsigned int k=0;k<regsize;k++)
			    {
                    tempNN.push_back(S[j].NN[k]);
                    cout << "tempNN " << S[j].NN[k] << " k " << k  << endl;
                }
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(tempNN);
                cout << "total number of hexes counted " << countHex << endl;
                int idx = 0;
                for (auto &h : S[j].Hgrid->hexen)
                {
                    array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                    if (!h.boundaryHex()) {
                        iidisp.drawHex(h.position(),(h.d/2.0f),colour);
                         cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                    }
                    idx++;
                }
                cout << " just befoere draw vCoords " << endl;
                for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                    std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                    iidisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                }
                cout << " just after draw vCoords " << endl;
            }//end of loop over regions
            cout << "just before redraw display 1" << endl;
            iidisp.redrawDisplay();
		    cout << "just after redraw display 1" << endl;
		    cout << "just after to_string"<<endl;
		    iidisp.saveImage(logpath + "/nnField1.png");
            usleep (1000000); // one hundred seconds
		    cout << "just after saveImage 1" << endl;
		    iidisp.closeDisplay();
		    cout << "just after close display 1 i " << i << endl;
        }//end of loop on numprint drawing fields
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
    } // end of second time stepping
//code run at end of timestepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData fdata(fname);
	for (unsigned int j=0;j<NUMPOINTS;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        fdata.add_contained_vals(ccst,S[j].CC);
        fdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
    }
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);

    /*
     * Now populate the Tessellation regions with NN
     */
    cout << "first setting of centroids" << endl;
    // repopulate the regions
    for (int j=0;j<NUMPOINTS;j++)
    {
        double NNNorm;
        cout << "jsut before renewRegion J " << j << endl;
        NNNorm = M.renewRegion(j,S[j].Hgrid->hexen, S[j].NN);
        cout << "nnNorm after M.renewRegion " << NNNorm << endl;
    }
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewBoundary(j,S[j].Hgrid->hexen);
    }
    M.setRadialSegments();
    cout << "just before renewDissect first time" << endl;
    M.edges_clear();
    cout << "just before setting regionList 1" << endl;
    cout << "just before setting regionList 2" << endl;
    M.setregionList(nbrList);
    cout << "just before renewDissect " << endl;
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewDissect(j,0);
    }
    cout << "Edges size " << M.edges.size() << endl;

    afile << "centroids for zero morphin " << endl;
    for (int j=0; j<NUMPOINTS;j++){
        afile << " Region size " << M.regionHex[j].size() << endl;
    }
    vector <int> radiusDVector;
    vector <int> angleDVector;
	vector <double> angleVector;
    vector <double> radiusVector;
    int degreeRadius;
    int degreeAngle;
	double tempArea = 0;
	double tempPerimeter = 0;
    int angleOffset = 0;
    int radiusOffset = 0;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    degreeRadius = 0;
    degreeAngle = 0;
	double avDegreeAngle = 0;
	double avDegreeRadius = 0;
    double avAbsCorrelation = 0;
    double occupancy = 0;
    int countRegions = 0;
    int max_comp = NUMPOINTS * 3;

    for (int j=0;j<NUMPOINTS;j++) {
        if (M.regArea(j) != 0){
			int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            gfile << "before sectorize Dangle" << endl;
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            // radiusVector = M.meanzero_vector(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
			    radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			    newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			    if (newdegreeRadius > degreeRadius) {
                    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  cout << "just before randomcorrelate_edges morph1 " << endl;
		  M.random_correlate(max_comp,0);
          cout << "just before correlate_edges morph1 " << endl;
          M.correlate_edges();
          cout << "just after  correlate_edges morph1 " << endl;
          for (int j=0;j<NUMPOINTS;j++) {
              countRegions++;
              occupancy += M.regNNfrac(j);
              cout << "just after regNNfrac " << endl;
              //tempArea = M.regArea(j);
              //tempPerimeter = M.regPerimeter(j);
              cout << "just before renewcorrelate_edges morph1 " << endl;
              avAbsCorrelation += M.renewcorrelate_edges(j,0);
              cout << "just after renewcorrelate_edges morph1 " << endl;
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
        cout << "Regions counted in first analysis loop " << countRegions << endl;
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after solving on polygonal regions via ksSolver

    /*
     * now create an array of morphed regions, first morph
     */

    if (skipMorph) return 0;
    cout << "just before setting curved boundaries first morph" <<endl;
// section for solving on the curved boundaries
    S.resize(0);
    //set radius for creating circular regions also calculate the adjusted Dn values
    vector<double> DnVal;
    DnVal.resize(NUMPOINTS,0.0);
    vector<double> DchiVal;
    DchiVal.resize(NUMPOINTS,0.0);
    vector<double> DcVal;
    DcVal.resize(NUMPOINTS,0.0);
    vector<double> morph0Area;
    morph0Area.resize(NUMPOINTS,0.0);
	for (int j = 0;j<NUMPOINTS;j++) {
        morph0Area[j] = M.hexArea*M.regArea(j);
    }
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(1);
    cout << "just after setting curved boundaries " << M.curvedBoundary.size()<<endl;
    for (int j = 0;j<NUMPOINTS;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }
    afile << "first morph  setting of centroids" << endl;
    for (int j=0; j<NUMPOINTS;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }

    for (int j = 0;j<NUMPOINTS;j++) {
	//	for (auto h : S[j].Hgrid->hexen) {
	//	    cout << "hexNumber in region " << j <<  " h.vi " << h.vi << endl;
    //    }
        double area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] = Dchi * morph0Area[j] / area;
            DnVal[j] = Dn *  morph0Area[j] / area;
            DcVal[j] = Dc * morph0Area[j] / area;
        }
        else {
            DchiVal[j] = Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }

        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
// now draw the intial tesselation
    internalCount = 0;
    boundaryCount = 0;
    if (Lgraphics) {
        morph::Gdisplay mdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        mdisp.resetDisplay (fix, eye, rot);
        for (int j=0;j<NUMPOINTS;j++) {
        // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
// cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                    mdisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                    boundaryCount++;
                }
                else
                {
                    cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
//cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                //mdisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                mdisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                mdisp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }

        } // end of loop over regions
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                mdisp.drawHex (h.position(), (h.d/2.0f), cl_b);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                mdisp.drawHex (h.position(), (h.d/2.0f), cl_c);
            }
        }
       //  else
        // {
          //   isp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
         //}
      //end of for loop on j as region number
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        mdisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        mdisp.saveImage(logpath + "/Tesselation2.png");
        mdisp.closeDisplay();
    } //end of graphics section
// initialise the fields
    string gname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(gname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<NUMPOINTS;j++)
        {
		    std::string ccstr = "c" + to_string(j);
		    cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
			char * ccst = new char[ccstr.length()+1];
			std::strcpy(ccst,ccstr.c_str());
		    std::string nstr = "n" + to_string(j);
			char * nst = new char[nstr.length()+1];
			std::strcpy(nst,nstr.c_str());
			cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
	        ginput.read_contained_vals(ccst,S[j].CC);
	        ginput.read_contained_vals(nst,S[j].NN);
		}
	  }
      else
	  {
        for (unsigned int j=0;j<NUMPOINTS;j++)
		{
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        double choice = ruf.get();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    // begin second time stepping loop
    for (int i=0;i<numsteps;i++)
    {
	    int countHex = 0;
   	    for (int j = 0;j<NUMPOINTS;j++) //loop over regions
        {
     	    S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
		    if (i%100 == 0)
		    cout << "just after morph 1 octopus i "<< j << " NN[5] " << S[j].NN[5] << endl;
	    }
	    if (i%numprint == 0 && Lgraphics)
	    {
          //set up display
            morph::Gdisplay mmdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            mmdisp.resetDisplay (fix, eye, rot);
   	        for (int j = 0;j<NUMPOINTS;j++) //loop over regions
	        {
	 	        cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
		        cout << "in loop over regions " << " size is " << regsize << endl;
			    countHex += regsize;
		        vector<double> regionNN;
			    vector<double> tempNN;
		      //for (auto h : S[j].Hgrid->hexen)
			    for (unsigned int k=0;k<regsize;k++)
			    {
                    tempNN.push_back(S[j].NN[k]);
                    cout << "tempNN " << S[j].NN[k] << " k " << k  << endl;
                }
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(tempNN);
                cout << "total number of hexes counted " << countHex << endl;
                int idx = 0;
                for (auto &h : S[j].Hgrid->hexen)
                {
                    array<float,3> colour = morph::Tools::getJetColorF(regionNN[idx]);
                    if (!h.boundaryHex()) {
                        mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                    //    cout << "just after drawHex 2 colour "<< " value " << regionNN[idx] << endl;
                    }
                    idx++;
                }
                for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                    std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                    mmdisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                }
            }//end of loop over regions
    for (auto h : M.Hgrid->hexen) {
         if (M.Creg[h.vi] ==  1 ) {
             mmdisp.drawHex (h.position(), (h.d/2.0f), cl_b);
             boundaryCount++;
         }
         else if (M.Creg[h.vi] > 1) {
              mmdisp.drawHex (h.position(), (h.d/2.0f), cl_c);
         }
    }
            cout << "just before redraw display 1" << endl;
            mmdisp.redrawDisplay();
		    cout << "just after redraw display 1" << endl;
		    cout << "just after to_string"<<endl;
		    mmdisp.saveImage(logpath + "/nnField2.png");
            usleep (1000000); // one hundred seconds
		    cout << "just after saveImage 1" << endl;
		    mmdisp.closeDisplay();
		    cout << "just after close display 1 i " << i << endl;
        }//end of loop on numprint drawing fields
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
    } // end of second time stepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData gdata(gname);
	for (unsigned int j=0;j<NUMPOINTS;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        gdata.add_contained_vals(ccst,S[j].CC);
        gdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
     }
     gdata.add_val ("/Dchi", Dchi);
     gdata.add_val ("/Dn", Dn);
     gdata.add_val ("/Dc",Dc);
     cout << " just after writing data "  << endl;
     // write the NN and CC vals for each region
      gfile << endl << "analysis on first morphing iteration " << endl;   // repopulate the regions
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewRegion(j,S[j].Hgrid->hexen,S[j].NN);
    }
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewBoundary(j,S[j].Hgrid->hexen);
    }
    afile << "second setting of centroids" << endl;
    for (int j=0; j<NUMPOINTS;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
	cout << "just after populating the ksVector"<<endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
      for (int j=0;j<NUMPOINTS;j++) {
              if (M.regArea(j) != 0){
			      int regionCount = 0;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M.regArea(j);
                  tempPerimeter = M.renewRegPerimeter(j);

                  // digital version
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
                  gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                  angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                  degreeAngle = L.find_zeroAngle(angleVector,3);
                  gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  degreeRadius = L.find_zeroDRadius(radiusDVector);
                  gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
				  {
			          radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				      degreeRadius = newdegreeRadius;
		          }


                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                  regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT


    //why does this have no effect?
    // redissect the boundaries
    cout << "just before renewDissect second time" << endl;
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewDissect(j,1);
    }
    cout << "Edges size " << M.edges.size() << endl;


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  cout << "just after renewcorrelate_edges morph1 " << endl;
		  M.random_correlate(max_comp,1);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (int j=0;j<NUMPOINTS;j++) {
	        if (M.regArea(j) != 0)
			{
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);
              avAbsCorrelation += M.renewcorrelate_edges(j,1);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile2 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after first morphing
//begin second morphing
    if (skipMorph) return 0;
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(0);
    S.resize(0);
    cout << "just after setting curved boundaries second iteration" << M.curvedBoundary.size()<<endl;
    for (int j = 0;j<NUMPOINTS;j++)
    {
        std::pair<float, float> centroid =  M.baryCentre(j);
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }
    cout << "just after populating the ksVector"<<endl;
    for (int j = 0;j<NUMPOINTS;j++) {
        double area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] =  Dchi * morph0Area[j] /  area;
            DnVal[j] = Dn * morph0Area[j] / area;
            DcVal[j] = Dc * morph0Area[j] / area;
        }
        else {
            DchiVal[j] =  Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }
        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
    boundaryCount = 0;
    internalCount = 0;
// now draw the intial tesselation
    if (Lgraphics) {
        morph::Gdisplay ndisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        ndisp.resetDisplay (fix, eye, rot);
		for (int j=0;j<NUMPOINTS;j++)
		{
          // plot stuff here.
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
		      //cout << "h.boundaryHex in region " << j <<  " h.vi " << h.vi << endl;
		     // cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_a);
                boundaryCount++;
                }
                else
                {
		      //cout << "h.internal hex in region " << j <<  " h.vi " << h.vi << endl;
		      //cout << "region " << j <<  " x " << S[j].Hgrid->d_x[h.vi] << " y " << S[j].Hgrid->d_y[h.vi] << endl;
//                ndisp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
                    internalCount++;
                }
            }
            for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                ndisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                pos = {M.mCoords[j][i].first, M.mCoords[j][i].second, 0};
                ndisp.drawHex (pos,offset,4.0*hexWidth,cl_b);
            }
		//cout << "boundaryCount 2 "<<boundaryCount<< " internalCount 2 "<< internalCount <<endl;


        // Draw small hex at boundary centroid
        array<float,3> c;
        c[2] = 0;
        c[0] = S[j].Hgrid->boundaryCentroid.first;
        c[1] = S[j].Hgrid->boundaryCentroid.second;
        cout << "d/2: " << S[j].Hgrid->hexen.begin()->d/4.0f << endl;
        //ndisp.drawHex (c, offset, (S[j].Hgrid->hexen.begin()->d/2.0f), cl_a);
        cout << "boundaryCentroid x,y: " << c[0] << "," << c[1] << endl;

        }//end of loop on NUMPOINTS to plot boundaries
        for (auto h : M.Hgrid->hexen) {
            if (M.Creg[h.vi] ==  1 ) {
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_b);
                boundaryCount++;
            }
            else if (M.Creg[h.vi] > 1) {
                ndisp.drawHex (h.position(), (h.d/2.0f), cl_c);
            }
        }
      // red hex at zero
      //ndisp.drawHex (pos, 0.005, cl_aa);

        usleep (100000);
        cout << "before redrawDisplay 3 " << endl;
        ndisp.redrawDisplay();
        cout << "after redrawDisplay 3" << endl;
        usleep (100000); // one hundred seconds
        ndisp.saveImage(logpath + "/Tesselation3.png");
        ndisp.closeDisplay();
    }
// initialise the fields
        cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
        string hname = logpath + "/third.h5";
        if (Lcontinue) {
             morph::HdfData hinput (hname,1);
             for (int j=0;j<NUMPOINTS;j++){
                std::string nstr = "n" + to_string(j);
                char * nst = new char[nstr.length()+1];
                std::strcpy(nst,nstr.c_str());
                std::string ccstr = "c" + to_string(j);
                char * ccst = new char[ccstr.length()+1];
                cout << "labels "<< nstr <<" , " << ccstr<<endl;
                std::strcpy(ccst,ccstr.c_str());
                hinput.read_contained_vals(nst,S[j].NN);
                hinput.read_contained_vals(ccst,S[j].CC);
                cout<< "just after input of NN and CC1"<< endl;
       //   input.close();
            }
	 }
     else {
        for (int j=0;j<NUMPOINTS;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                double choice = ruf.get();
                if (choice > 0.5)
                {
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
                }
                else {
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
                }

                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                } //end of if on boundary distance
            }//end of loop over region
        }//end of loop over all regions
     } //end of else on Lcontinue
     cout <<  "just after field creation" << endl;


     //begin second time stepping loop
    int loopsteps = 0;
    for (int i=0;i<numsteps;i++)
	{
	  //cout << " head of second time stepping loop i " << i << endl;
        for (int j = 0;j<NUMPOINTS;j++) {
     	    S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
     	 //  S[j].step(dt, Dn, Dchi, Dc);
	    }

	    if ((i%numprint == 0) && Lgraphics) {
            int countHex = 0;
          //set up display
            morph::Gdisplay nndisp(900, 900, 0, 0, "morph 2 u run", rhoInit, 0.0, 0.0);
            nndisp.resetDisplay (fix, eye, rot);
            for (int j = 0;j<NUMPOINTS;j++) {
                cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
                cout << "in loop over regions " << " size is " << regsize << endl;
                countHex += regsize;
                vector<double> regionNN(regsize);
                vector<double> tempNN(regsize);
                int k = 0;
                for (auto h : S[j].Hgrid->hexen) {
                    tempNN[k] = S[j].NN[h.vi];
                    k++;
                }
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(tempNN);
                cout << "total number of hexes counted " << countHex << endl;
		      //cout << "just before drawHex 2 "<< h.vi << "normalNN " << normalNN[h.vi] << endl;
                for (auto h : S[j].Hgrid->hexen) {
		      //if (!h.boundaryHex())
			  //{
                    array<float,3> colour = morph::Tools::getJetColorF(regionNN[h.vi]);
                    nndisp.drawHex(h.position(),(h.d/2.0f),colour);
		        //cout << "just after drawHex 2"<<endl;
			   //}
			    }
                for (unsigned int i=0; i< M.vCoords[j].size();i++) {
                    std::array<float,3> pos = {M.vCoords[j][i].first, M.vCoords[j][i].second, 0};
                    nndisp.drawHex (pos,offset,4.0*hexWidth,cl_c);
                }
            }//end of loop over regions
            for (auto h : M.Hgrid->hexen) {
                if (M.Creg[h.vi] ==  1 ) {
                    nndisp.drawHex (h.position(), (h.d/2.0f), cl_b);
                    boundaryCount++;
                }
                else if (M.Creg[h.vi] > 1) {
                    nndisp.drawHex (h.position(), (h.d/2.0f), cl_c);
                }
            }
            cout << "just berore redraw display 2 i " << i << endl;
            nndisp.redrawDisplay();
		    cout << "just after redraw display 2 i " <<  i << endl;
            usleep (1000000); // one hundred seconds
		    nndisp.saveImage(logpath + "/nnField3.png");
		    cout << "just after saveImage i " << i  << endl;
            nndisp.closeDisplay();
		    cout << "just after closeDisplay i " << i  << endl;
        }//end of loop on numprint drawing fields
		cout << "after numprint loop " << i << endl;
        loopsteps = i;
    } // end of second time stepping
	cout << "after second time stepping loop " << loopsteps << endl;
       //nndisp.closeDisplay();

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData hdata(hname);
	cout <<"after creating hdata "<< hname << endl;
	for (int j=0;j<NUMPOINTS;j++) {
		std::string nstr = "n" + to_string(j);
	   	char * nst = new char[nstr.length()+1];
		std::strcpy(nst,nstr.c_str());
		std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
		std::strcpy(ccst,ccstr.c_str());
		cout <<"in third hdf5 write "<<nstr <<" , "<<ccstr<<endl;
        hdata.add_contained_vals(ccst,S[j].CC);
        hdata.add_contained_vals(nst,S[j].NN);
	}
    //data.add_contained_vals("X",M.X[0]);
    //data.add_contained_vals("Y",M.X[1]);
    // hdata.add_val ("/Dchi", Dchi);
    // hdata.add_val ("/Dn", Dn);
    // hdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;

    gfile << endl << "analysis on second morphing iteration " << endl;
    // repopulate the regions
    for (int j=0;j<NUMPOINTS;j++)
    {
        double result = M.renewRegion(j,S[j].Hgrid->hexen, S[j].NN);
        cout << "after renewRegion total NN" << result;
    }
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewBoundary(j,S[j].Hgrid->hexen);
    }
//    swap the radialAngles to the vCoords
    M.swapRadialSegments(false);
    // redissect the boundaries
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewDissect(j,2);
    }
    cout << "Edges size " << M.edges.size() << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (int j=0;j<NUMPOINTS;j++) {
	    if (M.regArea(j) != 0) {
	        int regionCount = 0;
            gfile<<"in the degree loop" << endl;
	        //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

           // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl;
            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
			for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
                radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                newdegreeRadius = L.find_zeroRadius(radiusVector,3);
                if (newdegreeRadius > degreeRadius) {
				    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT

//computing average values over all regions
    avDegreeAngle = 0;
	avDegreeRadius = 0;
	occupancy = 0;
    avAbsCorrelation = 0;
	tempPerimeter = 0;
	tempArea = 0;
	countRegions = 0;
	cout << "just before renewcorrelate_edges morph2 " << endl;
	M.random_correlate(max_comp, 2);
	cout << "just after random correlate_edges morph2 " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (int j=0;j<NUMPOINTS;j++) {
	    if (M.regArea(j) != 0){
            countRegions++;
            avAbsCorrelation += M.renewcorrelate_edges(j,2);
            occupancy += M.regNNfrac(j);
            cout << "after renNNfrac " << endl;
            tempArea = M.regArea(j);
            cout << " after regArea " << endl;
            tempPerimeter = M.regPerimeter(j);
            cout << " after regPerimeter " << endl;
            cout << "before sectorixe Dangle morph2 " << endl;
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            cout << "after sectorixe Dangle morph2 " << endl;
            degreeAngle = L.find_zeroDAngle(angleDVector);
            avDegreeAngle += degreeAngle;
            cout << "after sectorixe Dangle morph2 " << endl;
            //radial degree
            degreeRadius = 0;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            avDegreeRadius += degreeRadius;
            degfile3 << degreeAngle << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " " << tempPerimeter<<endl<<flush;
         } //end of if on non-zero regions
    } //end of loop on NUMPOINTs
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	occupancy = occupancy / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation  <<endl;
    return 0;
};
