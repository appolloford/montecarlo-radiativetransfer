#include <stdio.h>
#include <cstdlib>
#include <utility>
#include <vector>
#include <cmath>
#include <time.h>  
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#define FAIL 0
#define SUCCESS 1

#define DIM 3
#define LAYER 20
#define NPHOTON 10000

double DomainSize;
extern int    NX_TOT[DIM];

// class Patch{
//     private:
//         double 
//     public:
//         patch();
//         ~patch();
// }

using namespace std;

class Photon{
    public:
        double energy;
        double direction[DIM];
    public:
        Photon(double e, double* dir){
            energy = e;
            for (int i=0; i<DIM; i++) {
                direction[i] = dir[i];
            }
        };
        // ~Photon();
};

class Cell{
    public:
        double internalEnergy;
        double position[DIM];
        double opticalDepth;
        double absorptionCoeff;
        double forbackRatio;
        std::vector<Photon> photons;

        int initializePhotons();

        Cell() {
            internalEnergy = 0.0;
            for (int i=0; i<DIM; i++) position[i] = 0.0;
        };

        int setAbsorptionCoeff( double coeff ) {
            absorptionCoeff = coeff;
            return SUCCESS;
        };

        int setForBackRatio( double ratio ) {
            forbackRatio = ratio;
            return SUCCESS;
        };

        int scattering() {
            for (auto it = photons.begin(); it != photons.end(); it++) {
                internalEnergy += absorptionCoeff * it->energy;
                it->energy *= (1.0 - absorptionCoeff);
                if (rand()/RAND_MAX < 1.0/(1.0+forbackRatio)) {
                    for (int i=0; i<DIM; i++) {
                        it->direction[i] *= -1;
                    }
                }
            }
            return SUCCESS;
        };
    // public:
        // Cell();
        // ~Cell();
};

int Cell::initializePhotons() {
    photons.clear();
    double dir[3] = {0.0, 0.0, 1.0};
    for (int i=0; i<NPHOTON; i++) {
        Photon photon(1.0, dir);
        photons.push_back(photon);
    }
    return SUCCESS;
}

int main(int argc, char *argv[]){

    // double absorptionCoeff = 1.0;
    bool inhomogeneous = false;
    double absorptionCoeff = 0.2;
    double forbackRatio = 0.7;
    double reflect = 0.3;
    Cell exteriorSource, exteriorSink;
    Cell cells[LAYER];
    double opticalLen = 1.0;

    static struct option long_options[] = {
        { "albedo", required_argument, nullptr, 0 },
        { "reflect", required_argument, nullptr, 0 },
        { "for-back-ratio", required_argument, nullptr, 0 },
        { "inhomogeneous-medium", no_argument, nullptr, 'v'},
        { "help", no_argument, nullptr, 'h' },
        { nullptr, no_argument, nullptr, 0 }
    };

    int opt = 0;
    int longIndex = 0;
    opt = getopt_long( argc, argv, "vh", long_options, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'v':
                inhomogeneous = true;
                printf("inhomogeneous medium: on\n");
                break;
            case 'h':   /* fall-through is intentional */
            // case '?':
            //     display_usage();
            //     break;

            case 0:     /* long option without a short arg */
                if( strcmp( "albedo", long_options[longIndex].name ) == 0 ) {
                    printf("albedo=%lf \n",atof(optarg));
                    absorptionCoeff = 1.0-atof(optarg);
                }
                else if( strcmp( "reflect", long_options[longIndex].name ) == 0 ) {
                    printf("reflect=%lf \n",atof(optarg));
                    reflect = atof(optarg);
                }
                else if( strcmp( "for-back-ratio", long_options[longIndex].name ) == 0 ) {
                    printf("for-back-ratio=%lf \n",atof(optarg));
                    forbackRatio = atof(optarg);
                }
                break;

            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt_long( argc, argv, "h", long_options, &longIndex );
    }

    double opticalDepthArray[LAYER] = {
        0.025, 0.025, 0.025, 0.025, 0.025,
        0.025, 0.025, 0.025, 0.025, 0.025,
        0.075, 0.075, 0.075, 0.075, 0.075, 
        0.075, 0.075, 0.075, 0.075, 0.075
    };

    double accumDepth = 0.0;
    for (int i=0; i<LAYER; i++) {
        if (inhomogeneous) {
            accumDepth += opticalDepthArray[i];
            cells[i].opticalDepth = accumDepth;
        }
        else {
            cells[i].opticalDepth = opticalLen * (i+1)/LAYER;
        }
        // cells[i].opticalDepth = opticalLen * (i+1)/LAYER;
        cells[i].setAbsorptionCoeff(absorptionCoeff);
        cells[i].setForBackRatio(forbackRatio);
    }

    exteriorSource.initializePhotons();
    // cells[0].initializePhotons();

    // int dir[3] = {1, 0, 0};
    // cells[0].photons.clear();
    // for (int i=0; i<NPHOTON; i++) {
    //     Photon photon(1.0, dir);
    //     cells[0].photons.push_back(photon);
    // }

    /* initialize random seed: */
    srand (time(NULL));

    auto & p = exteriorSource.photons;
    int count = 0;
    for (auto it = p.begin(); it != p.end(); it++) {
        double currentOpticalDepth = 0.0;
        double initialener = it->energy;
        while ((it->energy) >= 1.0e-3*initialener )
        {
            double r = (double) (rand()/(double)RAND_MAX);
            double tau = -log(1.0-r);

            tau *= it->direction[2]; // check the direction after scattering
            // printf("%13.7e \n", tau);
            // printf("%lf %lf %lf\n", tau, it->energy, currentOpticalDepth);

            if (currentOpticalDepth + tau > opticalLen) {
                if (reflect) {
                    exteriorSink.internalEnergy += (1.0-reflect) * it->energy;
                    it->energy = reflect * it->energy;
                    currentOpticalDepth = opticalLen;
                    for (int dim=0; dim<DIM; dim++) {
                        it->direction[dim] *= -1;
                    }
                }
                else {
                    exteriorSink.internalEnergy += it->energy;
                    it->energy = 0.0;
                }
                // exteriorSink.photons.push_back(std::move(*it));
                // p.erase(it--);
            }
            else if (currentOpticalDepth + tau < 0.0) {
                exteriorSource.internalEnergy += it->energy;
                it->energy = 0.0;
                // exteriorSink.photons.push_back(std::move(*it));
                // p.erase(it--);
            }
            else if (currentOpticalDepth + tau < cells[0].opticalDepth) {
                currentOpticalDepth = cells[0].opticalDepth;
                cells[0].internalEnergy += absorptionCoeff * it->energy;
                it->energy *= (1.0-absorptionCoeff);
                if ( it->energy < 1.0e-3*initialener ) {
                    // cells[0].photons.push_back(std::move(*it));
                    // p.erase(it--);
                }
                // // else if (rand()/(double)RAND_MAX < 1.0/(1.0+forbackRatio)) {
                // // else if (rand()/(double)RAND_MAX < forbackRatio) {
                // else if (rand()/(double)RAND_MAX < (1.0-forbackRatio)) {
                //     // for (int dim=0; dim<DIM; dim++) {
                //     //     it->direction[dim] *= -1;
                //     // }
                //     it->direction[2] = -1;
                // }
                // else {
                //     it->direction[2] = 1;
                // }
                else if (rand()/(double)RAND_MAX > forbackRatio) {
                    for (int dim=0; dim<DIM; dim++) {
                        it->direction[dim] *= -1;
                    }
                }
            }
            else {
                for (int i=1; i<LAYER; i++) {
                    if ( currentOpticalDepth + tau <= cells[i].opticalDepth &&
                        currentOpticalDepth + tau > cells[i-1].opticalDepth ) {

                        int idx = tau >= 0.0 ? i : i-1;

                        // printf("%lf, %lf, %d\n", tau, currentOpticalDepth, i);
                        currentOpticalDepth = cells[idx].opticalDepth;
                        cells[idx].internalEnergy += absorptionCoeff * it->energy;
                        it->energy *= (1.0-absorptionCoeff);
                        if ( it->energy < 1.0e-3*initialener ) {
                            // cells[i].photons.push_back(std::move(*it));
                            // p.erase(it--);
                        }
                        // // else if (rand()/(double)RAND_MAX < 1.0/(1.0+forbackRatio)) {
                        // // else if (rand()/(double)RAND_MAX < forbackRatio) {
                        // else if (rand()/(double)RAND_MAX < (1.0-forbackRatio)) {
                        //     // for (int dim=0; dim<DIM; dim++) {
                        //     //     it->direction[dim] *= -1;
                        //     // }
                        //     it->direction[2] = -1;
                        // }
                        // else {
                        //     it->direction[2] = 1;
                        // }
                        else if (rand()/(double)RAND_MAX > forbackRatio) {
                            for (int dim=0; dim<DIM; dim++) {
                                it->direction[dim] *= -1;
                            }
                        }
                        break;
                    }
                }
            }
        }
        
        // double r = (double) (rand()/(double)RAND_MAX);
        // double tau = -log(1.0-r);

        // for (int i=0; i<LAYER; i++) {
        //     if ( tau < cells[i].opticalDepth ) {
        //         cells[i].photons.push_back(std::move(*it));
        //         p.erase(it--);
        //     }
        // }

        // int l = min((int) (LAYER * tau / 1.0), LAYER-1);
        // if (l>0) {
        //     cells[l].photons.push_back(std::move(*it));
        //     p.erase(it--);
        // }

        // printf("%d\n", count);
        // count ++;
    }

    count = 0;
    FILE* outfile = fopen("res.dat", "w");
    double intensity = exteriorSink.internalEnergy / (double) NPHOTON;
    double absorption = 0.0;
    fprintf(outfile, "%lf %d \n", intensity, LAYER);
    for (int i=LAYER-1; i>=0; i--) {
        // count += cells[i].photons.size();
        // double intensity = (double) count / (double) NPHOTON;

        // double intensity = cells[i].internalEnergy / (double) NPHOTON;
        intensity += cells[i].internalEnergy / (double) NPHOTON;
        absorption += cells[i].internalEnergy / (double) NPHOTON;
        // printf("%lf %d \n", intensity, i);
        fprintf(outfile, "%lf %d \n", intensity, i);
    }
    printf("I(up)/I(tot) = %lf \n", exteriorSource.internalEnergy / (double) NPHOTON);
    printf("I(down)/I(tot) = %lf \n", exteriorSink.internalEnergy / (double) NPHOTON);
    printf("I(abs)/I(tot) = %lf \n", absorption);
    fclose(outfile);

    return 0;
}
