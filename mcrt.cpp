#include <stdio.h>
#include <cstdlib>
#include <utility>
#include <vector>
#include <cmath>
#include <time.h>
#include <random>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#define FAIL 0
#define SUCCESS 1

#define DIM 3
#define ITER_MAX 200
#define ITER_PRECISE 0.001

using namespace std;

double cellSizes[DIM]= {1.0, 1.0, 1.0};

class Photon{
    public:
        bool outofDomain;
        double initEner;
        double energy;
        double optlen;
        double direction[DIM];

        static random_device rd1;
        static random_device rd2;
        static random_device rd3;
        static random_device rd4;

        static default_random_engine lengen;
        // static default_random_engine scatgen;
        static std::mt19937_64 scatgen;
        static std::mt19937_64 scatgenII;
        static default_random_engine halfplaneanglegen;
        static default_random_engine rayleighanglegen;

        static uniform_real_distribution<double> lendistr;
        static uniform_real_distribution<double> scatdistr;
        static bernoulli_distribution angledistr;

    public:
        Photon() {
            outofDomain = false;
            initEner = 0.0;
            energy = 0.0;
            for (int i=0; i<DIM; i++) direction[i] = 0.0;
        };

        Photon(double e, double* dir){
            outofDomain = false;
            initEner = e;
            energy = e;
            for (int i=0; i<DIM; i++) {
                direction[i] = dir[i];
            }
        };

        bool isActive() {
            return (energy > 0.001 * initEner);
        }

        // make the maximum projection of direction equal to cellsize
        int normalizeDirection () {
            double maxratio = 0.0;
            double maxdir = 0.0;
            for (int dim=0; dim<DIM; dim++) {
                if (abs(direction[dim]/cellSizes[dim]) > maxratio) {
                    maxratio = abs(direction[dim]/cellSizes[dim]);
                    maxdir = direction[dim];
                }
            }
            for (int dim=0; dim<DIM; dim++) {
                direction[dim] /= maxratio;
            }
            return SUCCESS;
        }

        int RayleighScatteting () {
            double cosTheta = 1.0;
            double sinTheta = 0.0;

            double r1 = scatdistr(scatgenII);
            // double r1 = rand()/(double)RAND_MAX;

            // use Newton's method to solve the cos(theta) satisfying Rayleigh phase function integration
            for (int i=0; i<ITER_MAX; ++i) {
                double tmp = cosTheta;
                cosTheta = cosTheta - rayleighPhase(cosTheta, r1) / diffRayleighPhase(cosTheta);
                if (abs(tmp - cosTheta) < ITER_PRECISE) {
                    if (cosTheta > 1.0 || cosTheta < -1.0) {
                        printf("ITERATION ERROR: cos(theta) out of range.\n");
                        exit(1);
                    }
                    double sign = 2.0 * (double)angledistr(rayleighanglegen) - 1.0;
                    sinTheta = sign * sqrt(1.0 - pow(cosTheta, 2.0));
                    break;
                }
                else if (i+1 == ITER_MAX) {
                    printf("INTERATION FAILED: over maximal iteration\n");
                    exit(1);
                }
            }

            // rotation in xz plane
            direction[0] = cosTheta*direction[0] - sinTheta*direction[2];
            direction[2] = sinTheta*direction[0] + cosTheta*direction[2];

            return SUCCESS;
        };

        int halfPlaneScattering () {
            double sign = 2.0 * (double)angledistr(halfplaneanglegen) - 1.0;
            direction[0] = sign*scatdistr(scatgen);
            direction[2] = -sqrt(1.0 - pow(direction[0], 2.0));
            return SUCCESS;
        }

        double rayleighPhase (double t, double r1) {
            return 0.5 * (1.0 - 0.75*t - 0.25*pow(t, 3.0) ) - r1;
        }

        double diffRayleighPhase (double t) {
            return 0.3755 * (-1.0 - pow(t, 2.0) );
        }

        int resetOpticalLength () {
            optlen = -log( 1.0 - lendistr(lengen));
            return SUCCESS;
        }
        // ~Photon();
};

default_random_engine Photon::lengen(time(nullptr));
// default_random_engine Photon::scatgen(time(nullptr));
std::mt19937_64 Photon::scatgen(time(nullptr));
std::mt19937_64 Photon::scatgenII(time(nullptr));
default_random_engine Photon::halfplaneanglegen(time(nullptr));
default_random_engine Photon::rayleighanglegen(time(nullptr));

// default_random_engine Photon::lengen;
// default_random_engine Photon::scatgen;
// std::mt19937_64 Photon::scatgenII;
// default_random_engine Photon::halfplaneanglegen;
// default_random_engine Photon::rayleighanglegen;

std::uniform_real_distribution<double> Photon::lendistr(0.0,1.0);
std::uniform_real_distribution<double> Photon::scatdistr(0.0,1.0);
std::bernoulli_distribution Photon::angledistr(0.5);

class Cell{
    public:

        int nFrozenPhotons = 0;
        bool isExteriorSink;
        double internalEnergy;
        double position[DIM];
        // 3-dim symmetric tensor
        // double opticalDepth[6];
        double optDepth;
        double albedo;
        std::vector<Photon> injectPhotons;
        std::vector<Photon> ejectPhotons;

        // int initializePhotons( int nPhot, double ener ) {
        //     injectPhotons.clear();
        //     ejectPhotons.clear();
        //     double dir[DIM] = {0.0, 0.0, 1.0};
        //     for (int i=0; i<nPhot; i++) {
        //         Photon photon(ener, dir);
        //         injectPhotons.push_back(photon);
        //     }
        //     return SUCCESS;
        // }

        Cell() {
            internalEnergy = 0.0;
            optDepth = 0.0;
            for (int i=0; i<DIM; i++) position[i] = 0.0;
        };

        int activePhotons() {
            int n = 0;
            for (auto it = ejectPhotons.begin(); it != ejectPhotons.end(); it++) {
                if (it->isActive()) {
                    ++n;
                }
            }
            return n;
        }

        int setAlbedo(double alb) {
            albedo = alb;
            return SUCCESS;
        };

        int setOpticalDepth(double tau) {
            optDepth = tau;
            return SUCCESS;
        };

        int setPosition(double *pos) {
            for (int i=0; i<DIM; i++) position[i] = pos[i];
            return SUCCESS;
        };

        int scattering() {
            for (auto it = ejectPhotons.begin(); it != ejectPhotons.end(); it++) {
                if (it->isActive()) {
                    if (it->outofDomain) {
                        it->halfPlaneScattering();
                        it->outofDomain = false;
                    }
                    else {
                        it->RayleighScatteting();
                        internalEnergy += (1.0-albedo) * it->energy;
                        it->energy *= albedo;
                    }
                    it->normalizeDirection();
                }
                else {
                    internalEnergy += it->energy;
                    ejectPhotons.erase(it--);
                    nFrozenPhotons += 1;
                }
            }
            return SUCCESS;
        };

        int initializeNextTransition() {
            ejectPhotons.swap(injectPhotons);
            return SUCCESS;
        }

};

class IO{

    private:
        int cycleCounter = 0;

    public:
        void cycleIncrement(){
            cycleCounter++;
        }

        void WriteOutput(Cell *cells, int nCells){
            char filename[100];
            FILE *dump;

            sprintf(filename, "Cube_Ascii_%6.6d.dat", cycleCounter);
            dump = fopen(filename, "w");

            for (int i=0; i<nCells; i++) {
                fprintf(dump, "%13.7e %13.7e %13.7e ", cells[i].position[0], cells[i].position[1], cells[i].position[2]);
                fprintf(dump, "%13.7e %7.6d ", cells[i].internalEnergy, cells[i].nFrozenPhotons);
                fprintf(dump, "\n");
            }
            fclose(dump);

            sprintf(filename, "Array_photon_%6.6d.bin", cycleCounter);
            dump = fopen(filename, "w");

            for (int i=0; i<nCells; i++) {
                int np = cells[i].nFrozenPhotons+cells[i].ejectPhotons.size()+cells[i].injectPhotons.size();
                fwrite((void*)&np, sizeof(int), 1, dump);
            }
            fclose(dump);

            sprintf(filename, "Array_energy_%6.6d.bin", cycleCounter);
            dump = fopen(filename, "w");

            for (int i=0; i<nCells; i++) {
                double ie = cells[i].internalEnergy;
                fwrite((void*)&ie, sizeof(double), 1, dump);
            }
            fclose(dump);
        }
};

IO io;

int main(int argc, char *argv[]){

    double domainMin[DIM] = {0.0, 0.0, 0.0};
    double domainMax[DIM] = {10.0, 1.0, 20.0};
    int NX_TOT[DIM] = {10, 1, 20};
    int nPhotons = 100;

    int nCells = 1;
    // TODO: simplify codes, cells must be square
    for (int i=0; i<DIM; ++i) {
        nCells *= NX_TOT[i];
        cellSizes[i] = (domainMax[i] - domainMin[i]) / (double) NX_TOT[i];
    }

    double albedo = 0.5;
    double albedo_surf = 0.8;
    double domainOpticalDepthZ = 1.0;

    static struct option long_options[] = {
        { "albedo", required_argument, nullptr, 0 },
        { "reflect", required_argument, nullptr, 0 },
        { "nPhotons", required_argument, nullptr, 0 },
        // { "for-back-ratio", required_argument, nullptr, 0 },
        // { "inhomogeneous-medium", no_argument, nullptr, 'v'},
        { "help", no_argument, nullptr, 'h' },
        { nullptr, no_argument, nullptr, 0 }
    };

    int opt = 0;
    int longIndex = 0;
    opt = getopt_long( argc, argv, "vh", long_options, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'v':
                break;
            case 'h':   /* fall-through is intentional */
                printf("usage: ./mcrt.exe [--albedo=<value>]\n");
                printf("                  [--reflect=<surface albedo>]\n");
                printf("                  [--nPhotons=<total photon number>]\n");
                exit(1);
            // case '?':
            //     display_usage();
            //     break;
 
            case 0:     /* long option without a short arg */
                if( strcmp( "albedo", long_options[longIndex].name ) == 0 ) {
                    albedo = 1.0-atof(optarg);
                }
                else if( strcmp( "reflect", long_options[longIndex].name ) == 0 ) {
                    albedo_surf = atof(optarg);
                }
                else if( strcmp( "nPhotons", long_options[longIndex].name ) == 0 ) {
                    nPhotons = atoi(optarg);
                }
                // else if( strcmp( "for-back-ratio", long_options[longIndex].name ) == 0 ) {
                //     printf("for-back-ratio=%lf \n",atof(optarg));
                //     forbackRatio = atof(optarg);
                // }
                break;
                 
            default:
                /* You won't actually get here. */
                break;
        }
         
        opt = getopt_long( argc, argv, "h", long_options, &longIndex );
    }

    printf("albedo=%lf \n",albedo);
    printf("surface albedo=%lf \n",albedo_surf);

    Cell* cells = new Cell [nCells];
    Cell sourceCell, sinkCell;
    for (int k=0; k<NX_TOT[2]; k++) {
    for (int j=0; j<NX_TOT[1]; j++) {
    for (int i=0; i<NX_TOT[0]; i++) {
        int idx = i + j*NX_TOT[0] + k*NX_TOT[1]*NX_TOT[0];
        double pos[DIM];
        pos[0] = (i + 0.5) * cellSizes[0];
        pos[1] = (j + 0.5) * cellSizes[1];
        pos[2] = (k + 0.5) * cellSizes[2];
        cells[idx].setPosition(pos);

        cells[idx].setAlbedo(albedo);
        cells[idx].setOpticalDepth(domainOpticalDepthZ/NX_TOT[2]);
    }}}
    sourceCell.setAlbedo(0.0);
    sinkCell.setAlbedo(albedo_surf);

    std::vector<Photon> initialPhotons;
    for (int i=0; i<nPhotons; i++) {
        double dir[DIM] = {0.0, 0.0, 1.0};
        Photon photonsample(1.0, dir);
        photonsample.normalizeDirection();
        initialPhotons.push_back(photonsample);
    }

    int pCount = 0;
    int avePhot = nPhotons/NX_TOT[0];
    for (auto it = initialPhotons.begin(); it != initialPhotons.end(); it++) {
        // must set a proper starting point
        double stepPos[DIM];
        int signDir[DIM];
        for (int dim=0; dim<DIM; dim++) {
            stepPos[dim] = 0.5 * (domainMin[dim]+domainMax[dim]);
            signDir[dim] = abs(it->direction[dim]) ? it->direction[dim]/abs(it->direction[dim]):0;
        }
        stepPos[0] = (pCount/avePhot + 0.5)*cellSizes[0];
        stepPos[2] = 0.0;
        pCount += 1;
        it->resetOpticalLength();
        // printf("optlen: %lf\n", it->optlen);

        bool propagating = true;
        double optPathLen = 0.0;
        int currentIdx=-1;
        int stepIdx[DIM];
        while (propagating) {
            double stepLen, stepSize[DIM];
            // limited the maximum step to one cell
            double minStep = cellSizes[2];
            int minDir;
            for (int dim=0; dim<DIM; dim++) {
                stepIdx[dim] = stepPos[dim]/cellSizes[dim];
                stepSize[dim] = (stepIdx[dim]+signDir[dim])*cellSizes[dim] - stepPos[dim];
                if (it->direction[dim] && abs(stepSize[dim]) < minStep) {
                    minStep = abs(stepSize[dim]);
                    minDir = dim;
                }
            }
            stepLen = 0.0;
            for (int dim=0; dim<DIM; dim++) {
                stepSize[dim] = it->direction[dim] * minStep/cellSizes[minDir];
                stepLen += stepSize[dim] * stepSize[dim];
            }
            // printf("stepPos: %lf, %lf, %lf\n", stepPos[0], stepPos[1], stepPos[2]);
            // printf("stepSize: %lf, %lf, %lf\n", stepSize[0], stepSize[1], stepSize[2]);
            // printf("stepIdx: %d, %d, %d\n", stepIdx[0], stepIdx[1], stepIdx[2]);
            // printf("signDir: %d, %d, %d\n", signDir[0], signDir[1], signDir[2]);
            // printf("optlen: %lf\n", it->optlen);
            stepLen = sqrt(stepLen);
            int idx = stepIdx[0] + stepIdx[1]*NX_TOT[0] + stepIdx[2]*NX_TOT[1]*NX_TOT[0];
            optPathLen += stepLen * cells[idx].optDepth;
            // printf("cell idx: %d, stepLen: %lf, optPath: %lf\n", idx, stepLen, optPathLen);
            for (int dim=0; dim<DIM; dim++) {
                stepPos[dim] += stepSize[dim];
                stepIdx[dim] = stepPos[dim]/cellSizes[dim];
            }

            idx = stepIdx[0] + stepIdx[1]*NX_TOT[0] + stepIdx[2]*NX_TOT[1]*NX_TOT[0];
            if (stepPos[2] >= domainMax[2]) {
                it->outofDomain = true;
                sinkCell.internalEnergy += it->energy * (1.0 - sinkCell.albedo);
                it->energy *= sinkCell.albedo;
                stepPos[2] -= stepSize[2];
                idx -= NX_TOT[1]*NX_TOT[0];
                cells[idx].injectPhotons.push_back(std::move(*it));
                initialPhotons.erase(it--);
                propagating = false;
            }
            else if (stepPos[2] <= domainMin[2]) {
                sourceCell.internalEnergy += it->energy;
                it->energy = 0.0;
                sourceCell.injectPhotons.push_back(std::move(*it));
                initialPhotons.erase(it--);
                propagating = false;
            }
            else if (optPathLen >= it->optlen) {
                cells[idx].injectPhotons.push_back(std::move(*it));
                initialPhotons.erase(it--);
                propagating = false;
            }

            if (stepPos[0] < 0.0) {
                stepPos[0] += domainMax[0];
            }
            else if (stepPos[0] > domainMax[0]) {
                stepPos[0] -= domainMax[0];
            }

        }

    }

    io.WriteOutput(cells, nCells);
    io.cycleIncrement();

    int count = 0;
    // for (int i=nCells-1; i>=0; i--) {
    //     count += cells[i].injectPhotons.size();
    //     printf("# of photons: %d\n", count);
    // }

    // for (int i=0; i<nCells; i++) {
    //     count = cells[i].injectPhotons.size();
    //     printf("# of photons: %d\n", count);
    // }

    FILE *scat = fopen("firstscatter.dat", "w");
    for (int i=0; i<nCells; i++) {
        fprintf(scat, "%13.7e %13.7e %13.7e ", cells[i].position[0], cells[i].position[1], cells[i].position[2]);
        fprintf(scat, "%13.7e %13.7lu ", cells[i].internalEnergy, cells[i].injectPhotons.size());
        fprintf(scat, "\n");
    }
    fclose(scat);

    bool evolving = true;
    printf("cell evolution start\n");
    while (evolving) {

        // printf("new cycle!\n");
        
        evolving = false;
        for (int i=0; i<nCells; i++) {
            cells[i].initializeNextTransition();
        }

        for (int i=0; i<nCells; i++) {
            cells[i].scattering();
            auto & plist = cells[i].ejectPhotons;
            evolving = evolving || cells[i].activePhotons();
            for (auto it = plist.begin(); it != plist.end(); it++) {
                double stepPos[DIM];
                int signDir[DIM];
                int dirIdx[DIM];
                for (int dim=0; dim<DIM; dim++) {
                    // use the center of cell as starting point
                    stepPos[dim] = cells[i].position[dim];
                    // find the rough direction of ray
                    signDir[dim] = abs(it->direction[dim]) ? it->direction[dim]/abs(it->direction[dim]):0;
                    dirIdx[dim] = (it->direction[dim]>0.0) ? it->direction[dim]/abs(it->direction[dim]):0;
                }
                it->resetOpticalLength();

                bool propagating = true;
                double optPathLen = 0.0;
                int currentIdx=-1;
                int stepIdx[DIM];

                while (propagating) {
                    double stepLen, stepSize[DIM];
                    // limited the miaximum step to one cell
                    double minStep = cellSizes[2];
                    int minDir;
                    for (int dim=0; dim<DIM; dim++) {
                        stepIdx[dim] = stepPos[dim]/cellSizes[dim];
                        // stepSize[dim] = (stepIdx[dim]+signDir[dim])*cellSizes[dim] - stepPos[dim];
                        stepSize[dim] = (stepIdx[dim]+dirIdx[dim])*cellSizes[dim] - stepPos[dim];
                        if (it->direction[dim] && abs(stepSize[dim]/it->direction[dim]) && abs(stepSize[dim]/it->direction[dim]) < minStep) {
                            minStep = abs(stepSize[dim]/it->direction[dim]);
                            minDir = dim;
                        }
                    }
                    // printf("stepSize: %lf, %lf, %lf, minStep: %lf\n", stepSize[0], stepSize[1], stepSize[2], minStep);
                    stepLen = 0.0;
                    for (int dim=0; dim<DIM; dim++) {
                        stepSize[dim] = it->direction[dim] * minStep/cellSizes[minDir];
                        stepLen += stepSize[dim] * stepSize[dim];
                        // if (abs(stepSize[dim]) > stepLen) stepLen = abs(stepSize[dim]);
                    }
                    // printf("stepSize: %lf, %lf, %lf\n", stepSize[0], stepSize[1], stepSize[2]);
                    // printf("stepIdx: %d, %d, %d\n", stepIdx[0], stepIdx[1], stepIdx[2]);
                    // printf("signDir: %d, %d, %d\n", signDir[0], signDir[1], signDir[2]);
                    // printf("photonDir: %lf, %lf, %lf\n", it->direction[0], it->direction[1], it->direction[2]);
                    // printf("optlen: %lf\n", it->optlen);
                    stepLen = sqrt(stepLen);
                    int idx = stepIdx[0] + stepIdx[1]*NX_TOT[0] + stepIdx[2]*NX_TOT[1]*NX_TOT[0];
                    optPathLen += stepLen * cells[idx].optDepth;
                    // printf("cell idx: %d, stepLen: %lf, optPath: %lf\n", idx, stepLen, optPathLen);
                    for (int dim=0; dim<DIM; dim++) {
                        stepPos[dim] += stepSize[dim];
                        stepIdx[dim] = stepPos[dim]/cellSizes[dim];
                    }

                    idx = stepIdx[0] + stepIdx[1]*NX_TOT[0] + stepIdx[2]*NX_TOT[1]*NX_TOT[0];
                    // printf("stepPos: %lf, %lf, %lf\n", stepPos[0], stepPos[1], stepPos[2]);
                    // printf("photonDir: %lf, %lf, %lf\n", it->direction[0], it->direction[1], it->direction[2]);
                    // printf("cell idx: %d, stepIdx: %d, %d, %d\n", idx, stepIdx[0], stepIdx[1], stepIdx[2]);
                    if (stepPos[2] >= domainMax[2]) {
                        it->outofDomain = true;
                        sinkCell.internalEnergy += it->energy * (1.0 - sinkCell.albedo);
                        it->energy *= sinkCell.albedo;
                        stepPos[2] -= stepSize[2];
                        idx -= NX_TOT[1]*NX_TOT[0];
                        cells[idx].injectPhotons.push_back(std::move(*it));
                        plist.erase(it--);
                        propagating = false;
                    }
                    else if (stepPos[2] <= domainMin[2]) {
                        sourceCell.internalEnergy += it->energy;
                        it->energy = 0.0;
                        sourceCell.injectPhotons.push_back(std::move(*it));
                        plist.erase(it--);
                        propagating = false;
                    }
                    else if (optPathLen >= it->optlen) {
                        cells[idx].injectPhotons.push_back(std::move(*it));
                        plist.erase(it--);
                        propagating = false;
                    }

                    if (stepPos[0] <= domainMin[0] && signDir[0] < 0) {
                        stepPos[0] += domainMax[0];
                    }
                    else if (stepPos[0] >= domainMax[0] && signDir[0] > 0) {
                        stepPos[0] -= domainMax[0];
                    }

                }
            }
        }
        // count = 0;
        // for (int i=0; i<nCells; i++) {
        //     count += cells[i].injectPhotons.size();
        // }
        // printf("check the total injections: %d\n", count);

        io.WriteOutput(cells, nCells);
        io.cycleIncrement();
    }

    double internalEner = 0.0;
    for (int i=0; i<nCells; i++) {
        internalEner += cells[i].internalEnergy;
    }

    // WriteOutput(&cells[0], nCells, 99);
    io.WriteOutput(cells, nCells);

    printf("I_up/I_tot : %lf \n", sourceCell.internalEnergy/nPhotons);
    printf("I_sfc_abs/I_tot : %lf \n", sinkCell.internalEnergy/nPhotons);
    printf("I_abs/I_tot : %lf \n", internalEner/nPhotons);
    printf("Energy Conservation rate: %lf\n", (sourceCell.internalEnergy+sinkCell.internalEnergy+internalEner)/nPhotons);

    return 0;
}
