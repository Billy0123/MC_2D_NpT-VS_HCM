#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "MTGenerator.h"

/*State of optional functions:  (remember to always update text below after changing this file)
 -analyticMethodForHCM - INACTIVE
 -analyticCheck - INACTIVE
 -analyticMinDistanceCheckForHighAnisotropy - INACTIVE    (przy highAnisotropy - dla cząstek NIEgwiaździstych metoda 'sama' analityczna nie działa (zwraca zawsze 'wewnętrzne' minDist), a 'sama moja' może mieć przeskoki do 'wnętrza cząstek' (w niskich gęstościach, gdy cząstki mogą dokonywać 'zahaczeń' przy takiej wzajemnej orientacji, że w dalszym ciągu są od siebie DALEJ niż minDist). Metoda łączona to nic innego, jak metoda 'moja', ale z następującym przed nią samą warunkiem na minDist w funkcji kąta (a nie o stałej wartości minDist dla kąta krytycznego)
 -central interaction from KVT PhD - INACTIVE (może być aktywne tylko wtedy, gdy analyticMinDistanceCheckForHighAnisotropy jest NIEaktywne)
 -adjustAngles - ACTIVE
 -adjustOddRowsTranslation - INACTIVE
 -rotationsAnalysis - ACTIVE
 -utrzymywanie kata w granicy od C/3 do C - INACTIVE
 -rectangularMovesOnly - INACTIVE
 -matrix - INACTIVE
*/

int N,gaps,activeN,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,multimerN,countCollidingPairs,ODFLength,OCFMode,skipFirstIteration,saveConfigurations,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    onlyMath[2]={0,0},initMode,neighUpdatingFrequency;
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,intervalOrientations,savedConfigurationsInt,generatorStartPoint=0;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg[2],deltaR,deltaPhi,deltaV=0.1, //*sigma(=1)
       multimerS,multimerD,randomStartStep[2],detBoxMatrix,
       neighRadius,neighRadius2,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][3],pi=M_PI; //Problem typu-Pi (program przestawał działać dla wyższej precyzji Pi) to problem precyzji double. AKURAT dla 0.500001 się on zaczyna (przy 0.500010 jeszcze wszystko działa dobrze). Chodzi o to, że dla programu w C++ d=0.500001 to tak naprawdę: 0.500001000000000023 (np.), a dla takiej średnicy, procedury zwracają już INNE absoluteMinimum i minDistance. Normalnie cyfry tak 'dalekie' mogą być pominięte, ale przy d=0.500001 zaczyna już to odgrywać rolę. Zmniejszenie precyzji PI tak naprawdę ZMNIEJSZAŁO PI (ucięcie ostatnich liczb rozwinięcia) co jakoś przekładało się na to, że program zwracał 'wyglądające na dobre' wyniki.
       //UWAGA-update: coś dziwnego dalej się dzieje dla niektórych cząstek z wysoką anizotropią. Program działa poprawnie dla d/s>=0.516000 lub 0.500100>=d/s>=0.500010. W przedziałach: 0.516>d/s>0.5001 oraz 0.500010>d/s dzieje się coś dziwnego, co sprawia, że minDist (analitycznie zwracany) jest BŁĘDNY (absoluteDistance jest zawsze OK). Na szczęście: cząstki 0.51 można badać samą moją metodą, a cząstki 0.50001 metodą łączoną.
       //WARTE ODNOTOWANIA: powyższe problemy dotyczą tylko metody analitycznej, która czasami zwraca minDist 'wewnętrzne' a czasami 'zewnętrzne' (dla cząstek NIEgwiaździstych są dwie takie wielkości), a powinna ZAWSZE zwracać 'wewnętrzne' (w Mathematice tak robi). Moja metoda sprawdza JEDYNIE przekrycia, przez co dla bardzo wysokich anizotropii może po prostu czasami pozwolić na 'zabroniony jump' do 'środka' sąsiedniej molekuły (jednak wówczas warunek dr<minDist powinien się załączyć - chyba, że w niskich gęstościach, gdy cząstki mogą dokonywać 'zahaczeń' przy takiej wzajemnej orientacji, że w dalszym ciągu są od siebie DALEJ niż minDist). Problem w tym, że te minDist obliczane jest z metody analitycznej, która, jak opisałem powyżej, dla niektórych anizotropii przestaje działać - zawsze można wbić te minDist na sztywno z Mathematici.
       //UWAGA-update#2: to jest jeszcze bardziej pojebane, dla kolejnych: 0.500001, 0.501, 0.502, 0.503, 0.504, 0.506... zupełnie losowo zwracane są poprawne/błędne wartości minDist. Zmniejszając dokładność pi, niektóre się naprawiają, a inne psują, naprzemian, zupełnie losowo... lepiej minDist z Mathematici (7sem/NajszybszaMetoda/analityczna.nb) wprowadzać po prostu dla badanych mD.
double L,C,ROkreguOpisanego,absoluteMinimum,absoluteMinimum2,minDistance,maxDistance,VcpPerParticle,dr[3];
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferMN[20],bufferMS[20],bufferMD[20],bufferFolderIndex[5],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     orientationsFileName[200]="Orientations",
     orientatCorrelFunFileName[200]="OrientatCorrelFun",
     orientationsResultsFileName[200]="OrientatRes",
     configurationsListFileName[200]="ConfigurationsList.txt",
     loadConfigurationsFileName[200]="Configurations",
     loadedJOBID[50]="j-none";
/* //matrix 1/15
int boxCellsN[200][200],boxCellsIndex[200][200][10];
double bCSize[2];
*/
/////////////////  PARTICLE functions{
typedef struct {
    double r[2], normR[2];  //x,y
    double phi;   //kąt mierzony od kierunku x
    int neighbours[15], neighCounter;   //dla N=49920 trzeba zmniejszyć (z 50 np. na 10), jeżeli chce się działać na stacku. Uwaga na przepełnienie neighbours[] - wówczas przy tworzeniu listy sąsiadów program może 'pisać' w innych zmiennych, np. w r[] zamiast informować (już tak bywało)
    //int cell[2]; //matrix 2/15
} particle;


void getParticlesDistanceSquared (particle *p1, particle *p2, double boxMatrix[2][2]) {
    double normalizedRX=p1->normR[0]-p2->normR[0],
           normalizedRY=p1->normR[1]-p2->normR[1],
           rx=p1->r[0]-p2->r[0],
           ry=p1->r[1]-p2->r[1];
    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
    dr[0]=rx; dr[1]=ry; dr[2]=rx*rx+ry*ry;
}

/*void assignParticlesToCells (particle *particles, double boxMatrix[2][2]) { //matrix 3/15
    bCSize[0]=boxMatrix[0][0]/maxDistance; bCSize[1]=boxMatrix[1][1]/maxDistance;
    for (int i=0;i<bCSize[0];i++) for (int j=0;j<bCSize[1];j++) {
        boxCellsN[i][j]=0;
        for (int k=0;k<10;k++) boxCellsIndex[i][j][k]=-1;
    }
    for (int i=0;i<N;i++) {
        int indexes[2] = {(int)(particles[i].normR[0]*bCSize[0]),(int)(particles[i].normR[1]*bCSize[1])};
        boxCellsIndex[indexes[0]][indexes[1]][boxCellsN[indexes[0]][indexes[1]]++]=i;
        particles[i].cell[0]=indexes[0]; particles[i].cell[1]=indexes[1];
    }
}*/

void adjustNeighRadius (const double & volume) {
    neighRadius=1.6*sqrt(volume)/sqrt(N);
    neighRadius2=neighRadius*neighRadius;
}

void updateNeighbourList (particle *particles, double boxMatrix[2][2], const double & volume) {
    adjustNeighRadius(volume);
    for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        if (dr[2]<neighRadius2) {
            particles[i].neighbours[particles[i].neighCounter++]=j;
            particles[j].neighbours[particles[j].neighCounter++]=i;
        }
    }
}

double normalizeAngle (double phi) {//inaczej niż w Mathematice (loopedAngle[phi_]:=Mod[phi+C,2*C]-C), bo Mod[a,b] działa inaczej. Dla a>0 i b>0 nie ma rożnicy, ale dla a<0 i b>0 JEST. Mathematica zwraca wartości dodatnie (Mod[-1,10]=9), a C++ nie (fmod(-1,10)=-1).
    phi=fmod(phi+C,2*C);
    return phi<0?phi+C:phi-C;
}

void checkSinglePeriodicBoundaryConditions (particle *p, double boxMatrix[2][2]) {
    for (int j=0;j<2;j++) {
        for (int i=0;i<2;i++) p->r[i]-=floor(p->normR[j])*boxMatrix[i][j];
        p->normR[j]=fmod(p->normR[j],1); if (p->normR[j]<0) p->normR[j]++;
    }
    //particle->phi=normalizeAngle(particle->phi);    //rotationsAnalysis(comment) 1/1
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

double get2DiscsDistance (const int & aNum, const int & bNum, const double & dr, const double & aAngle, const double & bAngle) {
    double discAngle = aAngle+aNum*2*C;
    double aNumDiscX = cos(discAngle)*ROkreguOpisanego,
           aNumDiscY = sin(discAngle)*ROkreguOpisanego;
    discAngle = bAngle+bNum*2*C;
    double bNumDiscX = dr+cos(discAngle)*ROkreguOpisanego,
           bNumDiscY = sin(discAngle)*ROkreguOpisanego;
    double xDistance = aNumDiscX-bNumDiscX, yDistance = aNumDiscY-bNumDiscY;
    return sqrt(xDistance*xDistance+yDistance*yDistance);
}

int checkOverlaps (const double & dr, const double & aAngle, const double & bAngle) {
    int overlap=0;
    double startAAngle,startBAngle;
    int a,b;

    if (aAngle>absoluteMinimum2)
        if (aAngle>bAngle) {
            startAAngle=aAngle-2*C;
            a=2;
        } else {
            startAAngle=aAngle;
            a=1;
        }
    else if (aAngle>-absoluteMinimum2) {
        startAAngle=aAngle;
        a=1;
    } else
        if (aAngle<bAngle) {
            startAAngle=aAngle;
            a=2;
        } else {
            startAAngle=aAngle;
            a=1;
        }
    if (bAngle>absoluteMinimum2)
        if (a==2) {
            startBAngle=bAngle+(double)multimerN*C;
            b=1;
        } else {
            startBAngle=bAngle+((double)multimerN-2.0)*C;
            b=2;
        }
    else if (bAngle>-absoluteMinimum2) {
        startBAngle=bAngle+(double)multimerN*C;
        b=1;
    } else
        if (a==2) {
            startBAngle=bAngle+(double)multimerN*C;
            b=1;
        } else {
            startBAngle=bAngle+(double)multimerN*C;
            b=2;
        }
    for (int j=0;j<a;j++) for (int k=0;k<b;k++)
        if (get2DiscsDistance(j,k,dr,startAAngle,startBAngle)<multimerD) {
            overlap=1;
            j=a; k=b;
            break;
        }
    return overlap;
}

int checkOverlaps2DiscsMethod (const double & dr, const double & aAngle, const double & bAngle) {  //niska efektywność - można obliczyć najpierw położenia wszystkich sprawdzanych dysków, a potem tylko porównywać odległości zamiast kilkukrotnie obliczać położenie tych samych dysków metodą 'get2DiscsDistance'
    int overlap=0;
    double startAAngle=aAngle>0?aAngle-2*C:aAngle,
           startBAngle=bAngle>0?bAngle+(multimerN-2)*C:bAngle+multimerN*C;
    for (int j=0;j<2;j++) for (int k=0;k<2;k++)
        if (get2DiscsDistance(j,k,dr,startAAngle,startBAngle)<multimerD) {
            overlap=1;
            j=2; k=2;
            break;
        }
    return overlap;
}

int checkOverlaps3DiscsMethod (const double & dr, const double & aAngle, const double & bAngle) {  //niska efektywność - można obliczyć najpierw położenia wszystkich sprawdzanych dysków, a potem tylko porównywać odległości zamiast kilkukrotnie obliczać położenie tych samych dysków metodą 'get2DiscsDistance'
    int overlap=0;
    double startAAngle=aAngle-2*C,
           startBAngle=bAngle+(multimerN-2)*C;
    for (int j=0;j<3;j++) for (int k=0;k<3;k++)
        if (get2DiscsDistance(j,k,dr,startAAngle,startBAngle)<multimerD) {
            overlap=1;
            j=3; k=3;
            break;
        }
    return overlap;
}

int checkOverlapsAllDiscsMethod (const double & dr, const double & aAngle, const double & bAngle) {  //niska efektywność - można obliczyć najpierw położenia wszystkich sprawdzanych dysków, a potem tylko porównywać odległości zamiast kilkukrotnie obliczać położenie tych samych dysków metodą 'get2DiscsDistance'
    int overlap=0;
    for (int j=0;j<multimerN;j++) for (int k=0;k<multimerN;k++)
        if (get2DiscsDistance(j,k,dr,aAngle,bAngle)<multimerD) {
            overlap=1;
            j=multimerN; k=multimerN;
            break;
        }
    return overlap;
}

double minimalDistanceAnalyticalMethodHCM (const int & indeksPowierzchni, const double & aAngle, const double & bAngle, double a[4], double b[4]) {
    double angleA=aAngle+a[indeksPowierzchni],angleB=bAngle+b[indeksPowierzchni],
           buffer=sin(angleA)+sin(angleB);
    return ROkreguOpisanego*(cos(angleA)+cos(angleB))+sqrt(multimerD*multimerD-ROkreguOpisanego*ROkreguOpisanego*buffer*buffer);
}

double getMinimalDistanceAnalyticalMethodForEvenHCM (const double & aAngle, const double & bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-C,-C,C,C};
    if (aAngle<-absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13=-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<0) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve34=-C-ARCSIN;
        double intersectionCurve14=aAngle;
        if (bAngle>intersectionCurve14) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve34) return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve12=C-ARCSIN;
        double intersectionCurve14=aAngle;
        if (bAngle>intersectionCurve12) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve14) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    } else {
        double ARCSIN=asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24=-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    }
}

double getMinimalDistanceAnalyticalMethodForOddHCM (const double & aAngle, const double & bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-2*C,0,0,2*C};
    if (aAngle<-absoluteMinimum) {
        double ARCSIN = asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13 = C-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<0) {
        double ARCSIN = asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve23 = -ARCSIN;
        double intersectionCurve12 = aAngle+C;
        if (bAngle>intersectionCurve12) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve23) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN = asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve23 = -ARCSIN;
        double intersectionCurve34 = aAngle-C;
        if (bAngle>intersectionCurve23) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve34) return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    } else {
        double ARCSIN = asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24 = -C-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    }
}

int checkOverlapsAnalyticalMethodHCM (const double & dr, double aAngle, double bAngle) {
    int overlap=0;

    aAngle=normalizeAngle(aAngle+C);
    bAngle=normalizeAngle(bAngle+C);

    if (multimerN%2==0) {
        if (dr<getMinimalDistanceAnalyticalMethodForEvenHCM(aAngle,bAngle)) overlap=1;
    } else if (dr<getMinimalDistanceAnalyticalMethodForOddHCM(aAngle,bAngle)) overlap=1;
    return overlap;
}

int createRandomGaps (particle *particles, double boxMatrix[2][2], const double & volume) {
    printf("Creating %d random gaps... ",gaps); fflush(stdout);
    adjustNeighRadius(volume);
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTRandom0to1(randomStartStep)*N);
            for (int j=0;j<i;j++) {
                getParticlesDistanceSquared(&particles[gapsIndexes[i]],&particles[gapsIndexes[j]],boxMatrix);
                if (dr[2]<neighRadius2) {
                    i--;
                    innerAttempt++;
                    break;
                }
                if (j+1==i) innerAttempt=0;
            }
            if (innerAttempt>100000) break;
            if (innerAttempt==0 && i+1==gaps) allReady=1;
        }
    } while (!allReady && attempt<10000000);

    if (attempt>=10000000) {
        printf("ERROR: Couldn't create %d random gaps in %d steps.\n",gaps,attempt);
        return 0;
    } else {
        bool change; do {
            change=false;
            for (int i=gaps-1;i>0;i--) {
                if (gapsIndexes[i-1]>gapsIndexes[i]) {
                    int buffer=gapsIndexes[i];
                    gapsIndexes[i]=gapsIndexes[i-1]; gapsIndexes[i-1]=buffer;
                    change=true;
                }
            }
        } while (change);
        int actualGapIndex=0;
        for (int i=0;i<activeN;i++) {
            if (actualGapIndex<gaps) while (i+actualGapIndex==gapsIndexes[actualGapIndex]) {
                actualGapIndex++;
                if (actualGapIndex==gaps) break;
            }
            for (int j=0;j<2;j++) {
                particles[i].r[j]=particles[i+actualGapIndex].r[j];
                particles[i].normR[j]=particles[i+actualGapIndex].normR[j];
            }
        }
        printf("done\n"); fflush(stdout);
        return 1;
    }
}

int initPositions (particle *particles, double boxMatrix[2][2], double matrixOfParticlesSize[2], int n[2], double matrixCellXY[6][6][2], double matrixCellPhi[6][6], const double & pacFrac, const double & volume) {
    if (generatorStartPoint==0) {
        printf("Setting start position of p-random number generator to actual CPU time (for INIT PROCEDURES)...\n");
        InitRandomMT();
    } else {
        printf("Setting start position of p-random number generator to %ld (for INIT PROCEDURES)...\n",generatorStartPoint);
        InitMT((unsigned int)generatorStartPoint);
    }

    double mod=sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]), interval[2][2], actualPosition[2];
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) interval[i][j]=boxMatrix[i][j]/matrixOfParticlesSize[j]/mod*n[j];
    int rowCounter=0, columnCounter=0;
    for (int i=0;i<N;i++) {
        int cellNumber[2][2]={{columnCounter/n[0],rowCounter/n[1]},{columnCounter%n[0],rowCounter%n[1]}}; //cellNumber[0/1][X/Y]: 0-numer komorki, 1-kolumna/rzad W komorce
        actualPosition[0]=cellNumber[0][0]*interval[0][0]+cellNumber[0][1]*interval[0][1]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][0]*sqrt(pacFrac);
        actualPosition[1]=cellNumber[0][1]*interval[1][1]+cellNumber[0][0]*interval[1][0]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][1]*sqrt(pacFrac);

        for (int j=0;j<2;j++) particles[i].r[j]=actualPosition[j];
        particles[i].normR[0]=(boxMatrix[1][1]*particles[i].r[0]-boxMatrix[0][1]*particles[i].r[1])/detBoxMatrix;
        particles[i].normR[1]=-(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/detBoxMatrix;
        particles[i].phi=matrixCellPhi[cellNumber[1][0]][cellNumber[1][1]];

        columnCounter++;
        if (columnCounter*1.000001>=matrixOfParticlesSize[0]*mod) {
            rowCounter++;
            columnCounter=0;
        }
    }
    if (gaps>0) return createRandomGaps(particles,boxMatrix,volume);
    else return 1;
}

int getEnergyAll (particle *particles, double boxMatrix[2][2]) {
    int energy=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        double drRoot=sqrt(dr[2]);
        if (drRoot<maxDistance) {
            if (drRoot<minDistance) energy=1;   //analyticMinDistanceCheckForHighAnisotropy 1/3
            else {                              //
            //double gamma=atan(dr[1]/dr[0]),aAngle=particles[j].phi-gamma,bAngle=particles[i].phi-gamma; if ((energy=checkOverlapsAnalyticalMethodHCM(drRoot,dr[0]>=0?normalizeAngle(aAngle):normalizeAngle(bAngle),dr[0]>=0?normalizeAngle(bAngle):normalizeAngle(aAngle)))==0) { //analyticMinDistanceCheckForHighAnisotropy
                double gamma=atan(dr[1]/dr[0]),
                       aAngle=particles[j].phi-gamma,
                       bAngle=particles[i].phi-gamma;
                if (multimerN%2!=0) {
                    if (dr[0]>0) bAngle-=C;
                    else aAngle-=C;
                }
                aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                energy=checkOverlaps(drRoot,aAngle,bAngle);
            }
            if (energy==1) {
                i=activeN; j=activeN; break;
            }
        }
    }
    return energy;
}

void adjustAngles (particle *particles, double boxMatrix[2][2]) {
    int tryNumber=-1,energy;
    printf("Angle adjusting... "); fflush(stdout);
    do {
        tryNumber++;
        energy=getEnergyAll(particles,boxMatrix);
        if (energy==1) for (int k=0;k<activeN;k++) particles[k].phi+=0.00001;
    } while (energy!=0);
    printf("Adjusted after: %d approach.\n",tryNumber); fflush(stdout);
}

void adjustOddRowsTranslation (particle *particles, double boxMatrix[2][2], const int & partInRow) {
    int tryNumber=-1,energy;
    printf("Odd rows translation adjusting... "); fflush(stdout);
    do {
        tryNumber++;
        energy=getEnergyAll(particles,boxMatrix);
        if (energy==1) for (int k=0;k<activeN;k++) if ((k/partInRow)%2==1) {
            particles[k].r[0]+=0.00000001*multimerS;
            particles[k].normR[0]=(boxMatrix[1][1]*particles[k].r[0]-boxMatrix[0][1]*particles[k].r[1])/detBoxMatrix;
            particles[k].normR[1]=-(boxMatrix[1][0]*particles[k].r[0]-boxMatrix[0][0]*particles[k].r[1])/detBoxMatrix;
            checkSinglePeriodicBoundaryConditions(&particles[k],boxMatrix);
        }
    } while (energy!=0);
    printf("Adjusted after: %d approach.\n",tryNumber); fflush(stdout);
}

int getEnergy (particle *particles, particle *dispPart, const int & index, double boxMatrix[2][2]) {
    int energy=0;

    //matrix 4/15
    /*for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) {
        int cellIndex[2]={particles[index].cell[0]+i,particles[index].cell[1]+j};
        for (int k=0;k<2;k++) {
            if (cellIndex[k]<0) cellIndex[k]+=(int)bCSize[k];
            if (cellIndex[k]>=bCSize[k]) cellIndex[k]-=(int)bCSize[k];
        }
        for (int k=0;k<boxCellsN[cellIndex[0]][cellIndex[1]];k++) if (boxCellsIndex[cellIndex[0]][cellIndex[1]][k]!=index) {
            double normalizedRX=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].normR[0]-particles[index].normR[0],
                   normalizedRY=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].normR[1]-particles[index].normR[1],
                   rx=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].r[0]-particles[index].r[0],
                   ry=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].r[1]-particles[index].r[1];
            rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
            ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
            double r2=rx*rx+ry*ry, dr=sqrt(r2);
            if (dr<maxDistance) {
                double gamma=atan(ry/rx),
                       aAngle=particles[index].phi-gamma,
                       bAngle=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][k]].phi-gamma;
                if (multimerN%2!=0) {
                    if (rx>0) bAngle-=C;
                    else aAngle-=C;
                }
                aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                energy=(double)checkOverlaps(dr,aAngle,bAngle);
                if (energy==1) {k=boxCellsN[cellIndex[0]][cellIndex[1]];i=2;j=2;}
            }
        }
    }*/

    for (int i=0;i<particles[index].neighCounter;i++) {
        getParticlesDistanceSquared(&particles[particles[index].neighbours[i]],dispPart,boxMatrix);
        double drRoot=sqrt(dr[2]);
        if (drRoot<maxDistance) {
            if (drRoot<minDistance) energy=1;   //analyticMethodForHCM 1/8    //analyticMinDistanceCheckForHighAnisotropy 2/3
            else {                              //                            //
            //double gamma=atan(dr[1]/dr[0]),aAngle=particles[index].phi-gamma,bAngle=particles[particles[index].neighbours[i]].phi-gamma; if ((energy=checkOverlapsAnalyticalMethodHCM(drRoot,dr[0]>=0?normalizeAngle(aAngle):normalizeAngle(bAngle),dr[0]>=0?normalizeAngle(bAngle):normalizeAngle(aAngle)))==0) { //analyticMinDistanceCheckForHighAnisotropy
                double gamma=atan(dr[1]/dr[0]),
                       aAngle=dispPart->phi-gamma,
                       bAngle=particles[particles[index].neighbours[i]].phi-gamma;
                if (multimerN%2!=0) {
                    //rozważanie, która molekuła jest 'po lewej', a która 'po prawej' (moja metoda-odwraca cząstkę po prawej; analityczna metoda-istotna kolejność)
                    if (dr[0]>0) bAngle-=C;     //analytiCMethodForHCM 2/8
                    else aAngle-=C;             //
                    /*if (dr[0]<0) {
                        double buffer=aAngle; aAngle=bAngle; bAngle=buffer;
                    }*/
                }
                aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);   //analyticMethodForHCM 3/8
                energy=checkOverlaps(drRoot,aAngle,bAngle);                     //
                //energy=checkOverlapsAnalyticalMethodHCM(dr,aAngle,bAngle);
            } //analyticMethodForHCM 4/8
            if (energy==1) i=particles[index].neighCounter;
        }
    }
    return energy;
}

int attemptToDisplaceAParticle (particle *particles, const int & index, double boxMatrix[2][2]) {
    int result=1;
    particle displacedParticle;
    for (int i=0;i<2;i++) displacedParticle.r[i]=particles[index].r[i]+(MTRandom0to1(randomStartStep)-0.5)*deltaR;
    displacedParticle.normR[0]=(boxMatrix[1][1]*displacedParticle.r[0]-boxMatrix[0][1]*displacedParticle.r[1])/detBoxMatrix;
    displacedParticle.normR[1]=-(boxMatrix[1][0]*displacedParticle.r[0]-boxMatrix[0][0]*displacedParticle.r[1])/detBoxMatrix;
    displacedParticle.phi=particles[index].phi+(MTRandom0to1(randomStartStep)-0.5)*deltaPhi;
    //int oldCell[2]={particles[index].cell[0],particles[index].cell[1]}; //matrix 5/15
    //while (displacedParticle.phi>C || displacedParticle.phi<C/3.0) displacedParticle.phi=particles[index].phi+(MTRandom0to1(randomStartStep)-0.5)*deltaPhi;  //utrzymywanie kata w granicy od C/3 do C, zeby czastki sie ladnie ustawialy
    //checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix); //matrix 6/15
    //for (int i=0;i<2;i++) particles[index].cell[i]=(int)(particles[index].normR[i]*bCSize[i]); //matrix 7/15
    int newEnPot=getEnergy(particles,&displacedParticle,index,boxMatrix);
    if (newEnPot==1) {
        //for (int i=0;i<2;i++) particles[index].cell[i]=oldCell[i]; //matrix 8/15
        result=0;
    } else {
        for (int i=0;i<2;i++) {
            particles[index].r[i]=displacedParticle.r[i];
            particles[index].normR[i]=displacedParticle.normR[i];
        }
        particles[index].phi=displacedParticle.phi;
        checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix); //matrix(comment) 9/15
        /*if (particles[index].cell[0]!=oldCell[0] || particles[index].cell[1]!=oldCell[1]) { //matrix 10/15
            for (int i=0;i<boxCellsN[oldCell[0]][oldCell[1]];i++) if (boxCellsIndex[oldCell[0]][oldCell[1]][i]==index) {
                boxCellsIndex[oldCell[0]][oldCell[1]][i]=boxCellsIndex[oldCell[0]][oldCell[1]][--boxCellsN[oldCell[0]][oldCell[1]]];
                i=boxCellsN[oldCell[0]][oldCell[1]];
            }
            boxCellsIndex[particles[index].cell[0]][particles[index].cell[1]][boxCellsN[particles[index].cell[0]][particles[index].cell[1]]++]=index;
        }*/
    }
    return result;
}

void cloneParticlesForSpecificBoxMatrix (particle* clonedParticles, particle *particles, double boxMatrix[2][2], double xAxisPhiFactor) {
    for (int i=0;i<activeN;i++) {
        for (int j=0;j<2;j++) clonedParticles[i].normR[j]=particles[i].normR[j];
        for (int j=0;j<2;j++) clonedParticles[i].r[j]=boxMatrix[j][0]*particles[i].normR[0]+boxMatrix[j][1]*particles[i].normR[1];
        clonedParticles[i].phi=particles[i].phi+xAxisPhiFactor;
    }
}

int attemptToChangeVolume (particle *particles, const double & pressure, double boxMatrix[2][2], double & volume, double & xAxisPhi) {
    int result=1;
    double /*lnNewBoxMatrix[2][2], */newBoxMatrix[2][2];
    if (volume/VcpPerParticle/N<pressureRealOfNotFluid) {
    //if (pressure>pressureRealOfNotFluid) {  //dozwolone zmiany ksztaltu pudla (faza stala)
        /*for (int i=0;i<2;i++) for (int j=0;j<2;j++) lnNewBoxMatrix[i][j]=log(boxMatrix[i][j]);
        lnNewBoxMatrix[0][0]+=(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        lnNewBoxMatrix[1][1]+=(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        if (boxMatrix[0][1]==0) lnNewBoxMatrix[0][1]=-20.0+(MTRandom0to1(randomStartStep)-0.5)*deltaV; //log(0) -> -\infty; E^~0 za duze jak na start
        else lnNewBoxMatrix[0][1]+=(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        lnNewBoxMatrix[1][0]=lnNewBoxMatrix[0][1];
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) newBoxMatrix[i][j]=exp(lnNewBoxMatrix[i][j]);*/

        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[1][1]=boxMatrix[1][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;  //rectangularMovesOnly(comment second component on the right)
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        /*lnNewBoxMatrix[0][0]=log(boxMatrix[0][0])+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][0]=exp(lnNewBoxMatrix[0][0]);*/
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        double modifier=newBoxMatrix[0][0]/boxMatrix[0][0];
        newBoxMatrix[1][1]=boxMatrix[1][1]*modifier;
        newBoxMatrix[0][1]=boxMatrix[0][1]*modifier;
    }
    newBoxMatrix[1][0]=newBoxMatrix[0][1];
    double newDetBoxMatrix=newBoxMatrix[0][0]*newBoxMatrix[1][1]-newBoxMatrix[1][0]*newBoxMatrix[0][1],
           newVolume=fabs(newDetBoxMatrix), newXAxisPhi=atan(newBoxMatrix[0][1]/newBoxMatrix[0][0]);

    //matrix 11/15
    /*double bufferR[activeN][2];
    for (int i=0;i<activeN;i++) {
        bufferR[i][0]=newBoxMatrix[0][0]*particles[i].normR[0]+newBoxMatrix[0][1]*particles[i].normR[1];
        bufferR[i][1]=newBoxMatrix[1][0]*particles[i].normR[0]+newBoxMatrix[1][1]*particles[i].normR[1];
    }

//UWAGA: tutaj jest to zrobione na iterowaniu po komorkach, wyszlo jakies skomplikowane skanowanie tylko
//komorek z gory i z prawej, zeby nie powtarzac, jakies rozpoznawanie komorki 'srodkowej'... pojawia sie przez to duzo if'ow, etc.
//Moze lepiej zrobic analogicznie jak w przypadku listy sasiadow? Tzn. iterowac PO CZASTKACH, kazda sprawdzac z sasiadami z przyleglych
//komorek, a zeby nie bylo powtarzania, zastosowac ten sam prosty warunek: porownywac tylko, gdy index rozpatrywanej czastki jest mniejszy od porownywanej.

    for (int i=0;i<bCSize[0];i++) for (int j=0;j<bCSize[1];j++) {
        if (boxCellsN[i][j]>0) for (int k=0;k<5;k++) {
            int cellIndex[2];
            int zakresL=boxCellsN[i][j], startM=0;
            if (k==0) {cellIndex[0]=i; cellIndex[1]=j; zakresL--;}
            else if (k==1) {cellIndex[0]=i-1; cellIndex[1]=j+1;}
            else if (k==2) {cellIndex[0]=i; cellIndex[1]=j+1;}
            else if (k==3) {cellIndex[0]=i+1; cellIndex[1]=j+1;}
            else if (k==4) {cellIndex[0]=i+1; cellIndex[1]=j;}
            for (int l=0;l<2;l++) {

//UWAGA: tu chyba jest blad (cellIndex=-1 -> cellIndex=(int)bCSize-1, tymczasem komorka O INDEKSIE (int)bcSize TEZ moze cos zawierac, bo: .cell[i]=(int)(particles[index].normR[i]*bCSize[i])
//analogicznie, przy cellIndex=(int)bCSize+1 -> cellIndex=1 (zapomina sie o cellIndex=0). TEN SAM PROBLEM JEST przy matrix 4/16.
//Nalezy: sprawdzic dzialanie tego (wyswietlac cellIndex sprawdzanej czastki a po przecinku cellIndexy identyfikowane jako sasiednie.
//prawdopodobnie dla cellIndexu 0 i (int)bCSize blednie beda identyfikowane periodyczne sasiady.
//Poprawa kodu (na szybko myslac), nalezy dodawac/odejmowac nie (int)bCSize[l] a raczej ceil(bCSize[l]), bo (int)bCSize wskazuje na INDEKS
//najwyzszej komorki (tak zostalo to zaprojektowane), natomiast dopiero ceil(bCSize) wskazuje na index+1, czyli na LICZBE komorek (w pionie/poziomie).

                if (cellIndex[l]<0) cellIndex[l]+=(int)bCSize[l];
                if (cellIndex[l]>=bCSize[l]) cellIndex[l]-=(int)bCSize[l];
            }
            for (int l=0;l<zakresL;l++) {if (k==0) startM=l;
                for (int m=startM;m<boxCellsN[cellIndex[0]][cellIndex[1]];m++) {
                    double normalizedRX=particles[boxCellsIndex[i][j][l]].normR[0]-particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]].normR[0],
                           normalizedRY=particles[boxCellsIndex[i][j][l]].normR[1]-particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]].normR[1],
                           rx=bufferR[boxCellsIndex[i][j][l]][0]-bufferR[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]][0],
                           ry=bufferR[boxCellsIndex[i][j][l]][1]-bufferR[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]][1];
                    rx-=round(normalizedRX)*newBoxMatrix[0][0]+round(normalizedRY)*newBoxMatrix[0][1];
                    ry-=round(normalizedRX)*newBoxMatrix[1][0]+round(normalizedRY)*newBoxMatrix[1][1];
                    double r2=rx*rx+ry*ry, dr=sqrt(r2);
                    if (dr<maxDistance) {
                        double gamma=atan(ry/rx),
                               aAngle=particles[boxCellsIndex[cellIndex[0]][cellIndex[1]][m]].phi-gamma,
                               bAngle=particles[boxCellsIndex[i][j][l]].phi-gamma;
                        if (multimerN%2!=0) {
                            if (rx>0) bAngle-=C;
                            else aAngle-=C;
                        }
                        aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                        double energy=(double)checkOverlaps(dr,aAngle,bAngle);
                        if (energy==1) {
                            result=0;
                            l=zakresL; k=5; j=bCSize[1]; i=bCSize[0]; break;
                        }
                    }
                }
            }
        }
    }*/

    particle *particlesInNewBox = new particle[activeN];  //dla dużych N (testowane dla N=16384 w 3D) stack nie wystarcza, trzeba użyć heapu (stack zostaje dla podstawowej tablicy czastek, jest ciut szybszy)
    cloneParticlesForSpecificBoxMatrix(particlesInNewBox,particles,newBoxMatrix,newXAxisPhi-xAxisPhi);
    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
        if (i<particles[i].neighbours[j]) {
            getParticlesDistanceSquared(&particlesInNewBox[i],&particlesInNewBox[particles[i].neighbours[j]],newBoxMatrix);
            double drRoot=sqrt(dr[2]);
            if (drRoot<maxDistance) {
                int energy;
                if (drRoot<minDistance) energy=1; //analyticMethodForHCM 5/8    //analyticMinDistanceCheckForHighAnisotropy 3/3
                else {                            //                            //
                //double gamma=atan(dr[1]/dr[0]),aAngle=particles[particles[i].neighbours[j]].phi-gamma,bAngle=particles[i].phi-gamma; if ((energy=checkOverlapsAnalyticalMethodHCM(drRoot,dr[0]>=0?normalizeAngle(aAngle):normalizeAngle(bAngle),dr[0]>=0?normalizeAngle(bAngle):normalizeAngle(aAngle)))==0) { //analyticMinDistanceCheckForHighAnisotropy
                    double gamma=atan(dr[1]/dr[0]),
                           aAngle=particlesInNewBox[particles[i].neighbours[j]].phi-gamma,
                           bAngle=particlesInNewBox[i].phi-gamma;
                    if (multimerN%2!=0) {
                        //rozważanie, która molekuła jest 'po lewej', a która 'po prawej' (moja metoda-odwraca cząstkę po prawej; analityczna metoda-istotna kolejność)
                        if (dr[0]>0) bAngle-=C;     //analytiCMethodForHCM 6/8
                        else aAngle-=C;             //
                        /*if (dr[0]<0) {
                            double buffer=aAngle; aAngle=bAngle; bAngle=buffer;
                        }*/
                    }
                    aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle); //analyticMethodForHCM 7/8
                    energy=checkOverlaps(drRoot,aAngle,bAngle);                   //
                    //energy=checkOverlapsAnalyticalMethodHCM(drRoot,aAngle,bAngle);
                } //analyticMethodForHCM 8/8
                if (energy==1) {
                    result=0;
                    i=activeN; break;
                }
            }
        }
    }

    if (result) {
        double arg=-(pressure*(newVolume-volume)-(((double)N+1.0)*log(newVolume/volume)+log((newBoxMatrix[0][0]+newBoxMatrix[1][1])/(boxMatrix[0][0]+boxMatrix[1][1]))));  //kT=1[hard] przy czynniku ciśnieniowym (pominięte)
        if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
        else {
            volume=newVolume; xAxisPhi=newXAxisPhi;
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
            for (int i=0;i<activeN;i++) {
                for (int j=0;j<2;j++) particles[i].r[j]=particlesInNewBox[i].r[j];
                particles[i].phi=particlesInNewBox[i].phi;
            }
            detBoxMatrix=newDetBoxMatrix;
        }
    }

    delete [] particlesInNewBox; particlesInNewBox=NULL;
    return result;
}
/////////////////  } PARTICLE functions

int createIterationTable () {
    char startArguments[50]; FILE *fileStartArguments = fopen("startArguments.txt","rt");
    if (fileStartArguments==NULL) {printf("Missing file: startArguments.txt\n"); return 1;}
    while (fgets(startArguments,50,fileStartArguments)!=NULL) {
        sscanf(startArguments,"%c",startArguments); char *pEnd;
        iterationTable[fileIterateIterationsNumber][0]=strtod(startArguments,&pEnd);
        iterationTable[fileIterateIterationsNumber][1]=strtod(pEnd,&pEnd);
        iterationTable[fileIterateIterationsNumber++][2]=strtod(pEnd,NULL);
    }
    fclose(fileStartArguments);
    return 0;
}

void addAppendix (char *fileName, char *JOBID, bool jobIdOn) {
    strcpy(buffer,"2D_N-"); strncat(buffer,bufferN,20);
    strncat(buffer,"_gaps-",10); strncat(buffer,bufferGaps,20);
    strncat(buffer,"_G-",5); strncat(buffer,bufferG,5);
    strncat(buffer,"_badanie-",10); strncat(buffer,bufferFolderIndex,5);
    strncat(buffer,"_mN-",5); strncat(buffer,bufferMN,20);
    strncat(buffer,"_mS-",5); strncat(buffer,bufferMS,20);
    strncat(buffer,"_mD-",5); strncat(buffer,bufferMD,20);
    mkdir(buffer,S_IRWXU);
    strncat(buffer,"/",2);
    if (jobIdOn) {
        strncat(buffer,JOBID,50);
        strncat(buffer,"_",2);
    }
    strncat(buffer,fileName,200);
    strcpy(fileName,buffer);
}

void adjustOrientationsFile (FILE *file, char *path) {
    if (file==NULL) {
        file=fopen(path,"a"); fprintf(file,"{"); fclose(file);
    } else if (!onlyMath[0]) {
        char bufferForEraseLastChar[200],linia[110]; strcpy(bufferForEraseLastChar,path); strncat(bufferForEraseLastChar,"_BUFF",6);
        FILE *bFELC = fopen(bufferForEraseLastChar,"w");
        int poziomNawiasu=0;
        while (fgets(linia,100,file)!=NULL) {
            sscanf(linia,"%c",linia); int lastIndex=100;
            for (int i=0;i<100;i++) {
                if (linia[i]=='{') poziomNawiasu++;
                else if (linia[i]=='}') poziomNawiasu--;
                if (poziomNawiasu==0) {
                    lastIndex=i;
                    break;
                }
            }
            if (lastIndex==100) fprintf(bFELC,"%s",linia);
            else {
                for (int i=0;i<lastIndex;i++) fprintf(bFELC,"%c",linia[i]);
                fprintf(bFELC,"%c",',');
            }
        }
        fclose(bFELC); fclose(file);
        remove(path); rename(bufferForEraseLastChar,path);
    }
}

void getNextArgument (double prevArg[2], bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            for (int i=0;i<2;i++) prevArg[i]=growing?iterationTable[actIteration][i+1]:iterationTable[fileIterateIterationsNumber-1-actIteration][i+1];
            if (growing) startMinPacFrac=iterationTable[actIteration][0];
            else startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg[0]*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>=round(intervalMin[i]*10000) && round(prevArg[0]*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)+round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg[0]/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>round(intervalMin[i]*10000) && round(prevArg[0]*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)-round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)<round(minArg*10000)) growing=-1;
    }
}

inline bool isLineCorrect(char linia[4096]) {
    sscanf(linia,"%c",linia);
    int actIndex=0, dataIndex=0; while (dataIndex<3) {
        char data[50]="";
        int licznik=0, dotCounter=0;
        while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
        if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;} //test of single dot in a number
        actIndex++; dataIndex++;
        if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10; //test of dot position after first digit
    } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10; //test of max 3 numbers in a row

    return (dataIndex<10);
}

inline double getAvErrorFromSumEps (const double & sum, const double & denominator) {
    return sqrt(sum/denominator);
}

void updateTableAndGetActualMean (double table[100], double & mean, int const & changeIndex, double const & changeValue) {
    mean-=table[changeIndex]*0.01; table[changeIndex]=changeValue; mean+=changeValue*0.01;
}

int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[500];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,500,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20]="";
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++]; data[licznik]=' ';
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:multimerN=strtol(data,NULL,10);break;
                case 4:initMode=strtol(data,NULL,10);break;
                case 5:multimerS=strtod(data,NULL);break;
                case 6:multimerD=strtod(data,NULL);break;
                case 7:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 8:growing=strtol(data,NULL,10);break;
                case 9:loadedConfiguration=strtol(data,NULL,10);break;
                case 10:loadedArg=strtod(data,NULL);break;
                case 11:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 12:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 13:loadedSetGenerator=strtol(data,NULL,10);break;
                case 14:iterationsNumber=strtol(data,NULL,10);break;
                case 15:countCollidingPairs=strtol(data,NULL,10);break;
                case 16:intervalSampling=strtol(data,NULL,10);break;
                case 17:intervalOutput=strtol(data,NULL,10);break;
                case 18:saveConfigurations=strtol(data,NULL,10);break;
                case 19:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 20:ODFLength=strtol(data,NULL,10);break;
                case 21:OCFMode=strtol(data,NULL,10);break;
                case 22:neighUpdatingFrequency=strtol(data,NULL,10);break;
                case 23:intervalOrientations=strtol(data,NULL,10);break;
                case 24:skipFirstIteration=strtol(data,NULL,10);break;
                case 25:useSpecificDirectory=strtol(data,NULL,10);break;
                case 26:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 27:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 28:intervalResults=strtol(data,NULL,10);break;
                case 29:maxDeltaR=strtod(data,NULL);break;
                case 30:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 31:desiredAcceptanceRatioV=strtod(data,NULL);break;
                case 32:useFileToIterate=strtol(data,NULL,10);break;
                case 33:startMinPacFrac=strtod(data,NULL);break;
                case 34:startMaxPacFrac=strtod(data,NULL);break;
                case 35:minArg=strtod(data,NULL);break;
                case 36:maxArg=strtod(data,NULL);break;
                case 37:multiplyArgument=strtol(data,NULL,10);break;
                case 38:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-39)%3) {
                        case 0: intervalMin[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 1: intervalMax[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 2: intervalDelta[intervalLicznik++/3]=strtod(data,NULL);break;
                    }
                    break;
            }
            dataIndex++;
        }
        fclose(fileConfig);
    } while (testValue!=12345);

    //zlecanie parametrow z poziomu wiersza polecen:
    char JOBID[50]="j-"; int pointNumber=0;
    if (argumentsNumber==1) {
        strncat(JOBID,"none",50);
        if (useFileToIterate) if(createIterationTable()) return 0;
    } else {
        int correctNumberOfArguments=1;
        switch (strtol(arguments[1],NULL,10)) {
            case 0: //ustaw JOBID
                if (argumentsNumber==3) {
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (useFileToIterate) if(createIterationTable()) return 0;
                } else correctNumberOfArguments=0; break;
            case 1: //ustaw JOBID, singleRun dla parametrow zadanych bezposrednio
                if (argumentsNumber==13) {
                    useFileToIterate=0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    growing=strtol(arguments[9],NULL,10);
                    if (growing) {
                        startMinPacFrac=strtod(arguments[3],NULL); minArg=strtod(arguments[4],NULL);
                    } else {
                        startMaxPacFrac=strtod(arguments[3],NULL); maxArg=strtod(arguments[4],NULL);
                    }
                    N=strtol(arguments[5],NULL,10);
                    gaps=strtol(arguments[6],NULL,10);
                    multimerS=strtod(arguments[7],NULL);
                    multimerD=strtod(arguments[8],NULL);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    multimerN=strtol(arguments[12],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    N=strtol(arguments[3],NULL,10);
                    gaps=strtol(arguments[4],NULL,10);
                    multimerS=strtod(arguments[5],NULL);
                    multimerD=strtod(arguments[6],NULL);
                    growing=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    pointNumber=strtol(arguments[11],NULL,10);
                    generatorStartPoint=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1;
                    loadedArg=strtod(arguments[9],NULL);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1; loadType=1;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                    generatorStartPoint=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==13) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    onlyMath[0]=1;
                    onlyMath[1]=strtol(arguments[3],NULL,10);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=0;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    multimerN=strtol(arguments[12],NULL,10);
                    skipFirstIteration=0;
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, multimerN\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, generatorStartPoint, multimerN\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration, multimerN\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, multimerN, generatorStartPoint\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, multimerN\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        if (growing) {
            startMinPacFrac=iterationTable[0][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[0][i+1];
        } else {
            startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[fileIterateIterationsNumber-1][i+1];
        }
    } else {
        startArg[0]=growing?minArg:maxArg;
        startArg[1]=0;
    }

    //pressureRealOfNotFluid - dostosowanie do anisotropii, tutaj używane jako graniczne v* (w przypadku używania jako p* należy pamiętać o komendzie <<pressureRealOfNotFluid/=(multimerS*multimerS)>>)
    switch (multimerN) {
        case 3: switch ((int)round(multimerD/multimerS*100000)) {
            case 50001: pressureRealOfNotFluid=1.75; break; //0.50001
            case 51000: pressureRealOfNotFluid=1.50; break; //0.51
            case 55000: pressureRealOfNotFluid=1.36; break; //0.55
            case 60000: pressureRealOfNotFluid=1.31; break; //0.6
            case 66667: pressureRealOfNotFluid=1.27; break; //0.(6)
            case 80000: pressureRealOfNotFluid=1.23; break; //0.8
            case 100000: pressureRealOfNotFluid=1.19; break; //1
            case 120000: pressureRealOfNotFluid=1.16; break; //1.2
            case 150000: pressureRealOfNotFluid=1.25; break; //1.5
            case 200000: pressureRealOfNotFluid=1.25; break; //2
            case 300000: pressureRealOfNotFluid=1.24; break; //3
            case 500000: pressureRealOfNotFluid=1.24; break; //5
            case 140000: pressureRealOfNotFluid=1.25; break; //1.4
            case 135000: pressureRealOfNotFluid=1.25; break; //1.35
            case 130000: pressureRealOfNotFluid=1.25; break; //1.3
            case 96000: pressureRealOfNotFluid=1.19; break; //0.96
            case 92000: pressureRealOfNotFluid=1.19; break; //0.92
            case 88000: pressureRealOfNotFluid=1.19; break; //0.88
            case 84000: pressureRealOfNotFluid=1.19; break; //0.84
        } break;
    }

    deltaR=maxDeltaR*multimerS; deltaV*=multimerS;
    for (int i=0;i<pointNumber;i++) getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg[0];
    activeN=N-gaps;

    if (N%56!=0 && N%780!=0 && fabs(sqrt(N)-floor(sqrt(N)))>0.000001) {
        printf("ERROR: Not supported N: %d.\n",N);
        return 0;
    }

    //stale wynikajace z zadanych parametrow multimerow
    L=multimerS/multimerD;
    C=pi/(double)multimerN; deltaPhi=deltaR*2.0*sin(C)/multimerS;
    ROkreguOpisanego=multimerS/(2.0*sin(C));
    absoluteMinimum=atan(L/(2.0*L/tan(C)+sqrt(4.0-L*L)));
    absoluteMinimum2=C-absoluteMinimum;
    switch (multimerN) {
        case 6: switch ((int)(multimerD*1000000)) {
            case 500001: minDistance=1.8037364284584916; break;
            case 501000: minDistance=1.8331940077554698; break;
            case 502000: minDistance=1.845827233464671; break;
            case 503000: minDistance=1.8555403664338204; break;
            case 504000: minDistance=1.8637442860393154; break;
            case 505000: minDistance=1.8709851905229433; break;
            case 506000: minDistance=1.8775430588681767; break;
            case 507000: minDistance=1.8835841266820283; break;
            case 508000: minDistance=1.8892166507639028; break;
            case 509000: minDistance=1.8945157885020365; break;
            case 510000: minDistance=1.899536233850406; break;
            case 511000: minDistance=1.9043192554184105; break;
            case 512000: minDistance=1.9088968985500718; break;
            case 513000: minDistance=1.913294634275949; break;
            case 514000: minDistance=1.9175331038931438; break;
            case 515000: minDistance=1.9216293101072457; break;
            case 516000: minDistance=1.9255974552691306; break;
            case 517000: minDistance=1.9294495466426982; break;
            case 518000: minDistance=1.9331958432559069; break;
            case 519000: minDistance=1.9368451922349774; break;
            case 520000: minDistance=1.9404052862930774; break;
            default: minDistance=getMinimalDistanceAnalyticalMethodForEvenHCM(absoluteMinimum,absoluteMinimum); break;
        } break;
        default: if (multimerN%2==0) minDistance=getMinimalDistanceAnalyticalMethodForEvenHCM(absoluteMinimum,absoluteMinimum);
                 else minDistance=getMinimalDistanceAnalyticalMethodForOddHCM(absoluteMinimum,absoluteMinimum-C); break;
    }
    //minDistance=0.5*multimerS*sqrt(3)+sqrt(multimerD*multimerD-multimerS*multimerS*0.25);  //central interaction from KVT PhD 1/1
    maxDistance=ROkreguOpisanego*2+multimerD;

    /*printf("minDistance (n=%d, d=%.6f): %.17f\n",multimerN,multimerD,minDistance);    //testy porównawcze z Mathematicą (7 semestr.../1 Najszybsza metoda obliczania potencjału 2 multimerów/analityczna.nb). UWAGA: w przypadku rozbieżności, sprawdzić też absoluteMinimum! (bo z niego liczony jest absoluteMinimum2, który wykorzystywany jest przez moją metodę)
    return 0;*/

    //nazwy folderow na podstawie parametrow programu
    sprintf(bufferG,"%d",growing); sprintf(bufferN,"%d",N); sprintf(bufferGaps,"%d",gaps);
    sprintf(bufferMN,"%d",multimerN); sprintf(bufferMS,"%.2f",multimerS); sprintf(bufferMD,"%.6f",multimerD);

    int folderIndex=useSpecificDirectory, checkNext;
    char bufferCheckFolderExisting[200];
    FILE *checkFolderExisting;
    if (!folderIndex) do {
        sprintf(bufferFolderIndex,"%d",++folderIndex);
        strcpy(bufferCheckFolderExisting,resultsFileName);
        addAppendix(bufferCheckFolderExisting,JOBID,false);
        checkFolderExisting = fopen(bufferCheckFolderExisting,"rt");
        if (checkFolderExisting!=NULL) {
            fclose(checkFolderExisting);
            checkNext=1;
        } else checkNext=0;
    } while (checkNext);
    sprintf(bufferFolderIndex,"%d",folderIndex);
    addAppendix(resultsFileName,JOBID,false);
    addAppendix(excelResultsFileName,JOBID,false);
    addAppendix(configurationsFileName,JOBID,true);
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.4E",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(orientationsFileName,JOBID,true);
    if (OCFMode) addAppendix(orientatCorrelFunFileName,JOBID,true);
    addAppendix(orientationsResultsFileName,JOBID,true);
    addAppendix(configurationsListFileName,JOBID,false);

    particle particles[N];
    double args[10];

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileOrientations, *fileOrientatCorrelFun, *fileConfigurationsList, *fileAllResults, *fileAllOrientations, *fileOrientationsResults, *fileAllOrientationsResults;
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        if (saveConfigurations) fprintf(fileResults,"Cycles\tPressure*\tVolume\tdVolume\tBoxMatrix[0][0]\tdBoxMatrix[0][0]\tBoxMatrix[1][1]\tdBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tdBoxMatrix[1][0]([0][1])\tRho\tdRho\tV/V_cp\tdV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavNu2\tdAvNu2\tavB*\tdAvB*\tavMy*\tdAvMy*\tavE*\tdAvE*\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileResults,"Cycles\tPressure*\tVolume\tdVolume\tBoxMatrix[0][0]\tdBoxMatrix[0][0]\tBoxMatrix[1][1]\tdBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tdBoxMatrix[1][0]([0][1])\tRho\tdRho\tV/V_cp\tdV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavNu2\tdAvNu2\tavB*\tdAvB*\tavMy*\tdAvMy*\tavE*\tdAvE*\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        if (saveConfigurations) fprintf(fileExcelResults,"Pressure*\tV/V_cp\tavNu\tavNu2\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tavB*\tavMy*\tavE*\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileExcelResults,"Pressure*\tV/V_cp\tavNu\tavNu2\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tavB*\tavMy*\tavE*\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg[2]={startArg[0],startArg[1]}, oldBoxMatrix[2][2];
    while (growing>=0) {
        double pressureReduced=arg[0], pressure=pressureReduced/multimerS/multimerS, //pressureReduced=\tau*\sigma^2/kT, kT=1[hard] - Dwa punkty widzenia: 1) zmniejszenie sigma ZWIĘKSZA JEDNOSTKE p^*, zatem ten sam STAN FIZYCZNY jak przy sigma=1 bedzie przy mniejszym p^*. pReal to tak naprawde pReducedAdjusted. Inny, równoważny punkt widzenia, to 2) pReal redukuje objetosc, ktora NIE jest wyrazana w jednostkach sigma. Objetosc jest obliczana z boxMatrix, ktory jest inicjowany z czynnikiem *sigma, zatem MA jednostkę, a NIE jest zredukowany. Przeciez gdyby sigma=2, to boxMatrix bylby 2x wiekszy, a 'w jednostkach sigma' (zredukowany) powinien pozostac identyczny
               boxMatrix[2][2],matrixOfParticlesSize[2],unitCellAtCP[2],                                                                //obydwa sprowadzają się do tego, że przy liczeniu prawdopodobieństwa ma być jednostka zredukowana (bezwymiarowa): 1) zakłada, że volume jest zredukowane, więc dostosowuje pReduced do stanu fizycznego; 2) zakłada, że pReduced już jest OK (w końcu jest zredukowane), tylko po prostu objętość NIE jest zredukowana, i trzeba ją zredukować dzieląc przez sigma^2
               matrixCellXY[6][6][2],matrixCellPhi[6][6];
        int n[2]; //n[X/Y], matrixCell[n[XMax]][n[YMax]][x/y], zatem: n[X/Y](max)=6
        switch (multimerN) {
            case 6: {//dla heksamerow o dowolnym d/\sigma
                unitCellAtCP[0]=minDistance; unitCellAtCP[1]=sqrt(3)*minDistance;
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=absoluteMinimum2;
                matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=absoluteMinimum2;
            } break;
            case 5: {//dla pentamerow o d/\sigma=1
                unitCellAtCP[0]=2.4048671732*multimerS; unitCellAtCP[1]=4.2360679772*multimerS;
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                matrixCellXY[0][1][0]=1.0131106571*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
            } break;
            /*case 7: switch (initMode) {//dla heptamerow o d/\sigma=1, struktury jak w WojTreKow2003PRE
                case 0: {//struktura A
                    unitCellAtCP[0]=3.0566685376*multimerS; unitCellAtCP[1]=5.5150210832*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=1.9144263193*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//struktura B
                    unitCellAtCP[0]=3.0566685376*multimerS; unitCellAtCP[1]=5.5382990167*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=1.3656999121*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 2: {//struktura C
                    unitCellAtCP[0]=6.0309063912*multimerS; unitCellAtCP[1]=5.5371391576*multimerS;
                    n[0]=2; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.1644029030;
                    matrixCellXY[1][0][0]=unitCellAtCP[0]/2.0; matrixCellXY[1][0][1]=0.1230590461*multimerS; matrixCellPhi[1][0]=0.1644029030;
                    matrixCellXY[0][1][0]=1.2832327269*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.2843960475;
                    matrixCellXY[1][1][0]=4.2986859225*multimerS; matrixCellXY[1][1][1]=2.8916286250*multimerS; matrixCellPhi[1][1]=0.2843960475;
                } break;
            } break;*/
            case 3: switch (initMode) {//struktury trimerow jak z pracy KVT
                case 0: {//struktura INC (Isotropic Nonchiral Crystal)
                    unitCellAtCP[0]=0.5*multimerS*sqrt(3)+sqrt(multimerD*multimerD-multimerS*multimerS*0.25); unitCellAtCP[1]=unitCellAtCP[0]*sqrt(3);
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
                case 1: {//struktura R2 (Rectangular Crystal with 2 molecules in the unit cell)
                    unitCellAtCP[0]=0.5*multimerS*sqrt(3)+sqrt(multimerD*multimerD-multimerS*multimerS*0.25); unitCellAtCP[1]=2.0*cos(C/2.0-absoluteMinimum2)*minDistance;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=0;
                    matrixCellXY[0][1][0]=-sin(C/2.0-absoluteMinimum2)*minDistance; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 4: switch (initMode) {
                case 0: {//struktura kwadratowa skrecona
                    unitCellAtCP[0]=minDistance; unitCellAtCP[1]=minDistance;
                    n[0]=1; n[1]=1;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=absoluteMinimum2;
                } break;
            } break;
            /*case 95: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=31.1028*multimerS; unitCellAtCP[1]=54.084*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.5404*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie (przechodzi w DISORDERED nawet przy p*=4000)
                    unitCellAtCP[0]=31.1028*multimerS; unitCellAtCP[1]=54.1536*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.5423*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 97: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=31.7394*multimerS; unitCellAtCP[1]=55.1314*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.8679*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie
                    unitCellAtCP[0]=31.7394*multimerS; unitCellAtCP[1]=55.2852*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=15.8679*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 98: switch (initMode) {
                case 0: {//naprzemienne katy gamma (jednakowe gamma w rzedzie)
                    unitCellAtCP[0]=32.0695*multimerS; unitCellAtCP[1]=55.7368*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=0.0155995;
                    matrixCellXY[0][1][0]=16.0347*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.0485146;
                } break;
                case 1: {//jednakowe gamma wszedzie
                    unitCellAtCP[0]=32.0695*multimerS; unitCellAtCP[1]=55.7804*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=0.01559950076973243;
                    matrixCellXY[0][1][0]=16*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.01559950076973243;
                } break;
            } break;
            case 100: switch (initMode) {
                case 0: {//naprzemienne katy gamma (jednakowe gamma w rzedzie)
                    unitCellAtCP[0]=32.6904*multimerS; unitCellAtCP[1]=56.7792*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0161203;
                    matrixCellXY[0][1][0]=16.3452*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.0161203;
                } break;
                case 1: {//jednakowe gamma wszedzie
                    unitCellAtCP[0]=32.6904*multimerS; unitCellAtCP[1]=56.9516*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0161203;
                    matrixCellXY[0][1][0]=16.3452*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.0161203;
                } break;
            } break;
            case 101: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=33.0128*multimerS; unitCellAtCP[1]=57.4142*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.5064*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie
                    unitCellAtCP[0]=33.0128*multimerS; unitCellAtCP[1]=57.435*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.5064*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 103: switch (initMode) {
                case 0: {//naprzemienne katy o C (jednakowe katy w rzedzie)
                    unitCellAtCP[0]=33.6495*multimerS; unitCellAtCP[1]=58.4394*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.8248*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//jednakowe katy C wszedzie
                    unitCellAtCP[0]=33.6495*multimerS; unitCellAtCP[1]=58.5904*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=16.8248*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            case 104: switch (initMode) {
                case 0: {//naprzemienne katy gamma (jednakowe gamma w rzedzie)
                    unitCellAtCP[0]=33.963916981946*multimerS; unitCellAtCP[1]=59.062029530994*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0154856;
                    matrixCellXY[0][1][0]=16.982008490998*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0.0154856;
                } break;
                case 1: {//jednakowe gamma wszedzie
                    unitCellAtCP[0]=33.963916981946*multimerS; unitCellAtCP[1]=59.104*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.0154856;
                    matrixCellXY[0][1][0]=16.982008490998*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.0154856;
                } break;
            } break;*/
            default: {
                unitCellAtCP[0]=maxDistance; unitCellAtCP[1]=maxDistance*sqrt(3);
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
            } break;
        } VcpPerParticle=unitCellAtCP[0]*unitCellAtCP[1]/(double)n[0]/(double)n[1];
        if (N%56==0) {matrixOfParticlesSize[0]=7; matrixOfParticlesSize[1]=8;}
        else if (N%780==0) {matrixOfParticlesSize[0]=26; matrixOfParticlesSize[1]=30;}
        else if (floor(sqrt(N))==sqrt(N)) {matrixOfParticlesSize[0]=matrixOfParticlesSize[1]=sqrt(N);}
        double NLinearMod = sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]);
        if (startArg[0]==arg[0]) {
            for (int i=0;i<2;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]/(double)n[i]*NLinearMod*(growing?sqrt(startMinPacFrac):sqrt(startMaxPacFrac));
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0;
        } else for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
        detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
        double volume=fabs(detBoxMatrix), rho=N/volume, pacFrac=1.0/VcpPerParticle/rho, xAxisPhi=atan(boxMatrix[0][1]/boxMatrix[0][0]);

        if (!onlyMath[0]) {
            if (arg[0]==startArg[0] && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (StartDen: %.7E, startPacFrac: %.7E), mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,startArg[0],rho,pacFrac,multimerN,multimerS,multimerD);
                if (!initPositions(particles,boxMatrix,matrixOfParticlesSize,n,matrixCellXY,matrixCellPhi,pacFrac,volume)) return 0;
                adjustAngles(particles,boxMatrix);   //dla układów, w których obracanie cząstek pomaga uzyskać initowy układ (np. HCH)
                //adjustOddRowsTranslation(particles,boxMatrix,(int)(matrixOfParticlesSize[0]*NLinearMod));   //dla układów, w których translacja co drugiego rzędu pomaga uzyskać initowy układ (np. HCT)
                updateNeighbourList(particles,boxMatrix,volume); //matrix(comment) 12/15
            } else if (loadedConfiguration) {
                char configurations[4096];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                for (int i=0;i<3;i++) fgets(configurations,4096,fileCTL);
                int character,dataType=0,pIndex=-1; char data[50]=""; int actIndex=0;
                while (dataType<11) {
                    character=fgetc(fileCTL); //character is in int, but it can work as char
                    if (dataType<10) { //stage #1 configuration parameters
                        if (character==' ') {
                            data[actIndex++]=' '; //end of data without clearing the entire array
                            args[dataType++]=strtod(data,NULL);
                            if (dataType==10) {
                                boxMatrix[0][0]=args[5]; boxMatrix[1][1]=args[6]; boxMatrix[1][0]=args[7]; boxMatrix[0][1]=args[7];
                                deltaR=args[8]; deltaPhi=deltaR*2.0*sin(C)/multimerS; deltaV=args[9];
                                arg[0]=args[3]; pressureReduced=arg[0]; pressure=pressureReduced/multimerS/multimerS;

                                detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
                                volume=fabs(detBoxMatrix); rho=N/volume; pacFrac=1.0/VcpPerParticle/rho; xAxisPhi=atan(boxMatrix[0][1]/boxMatrix[0][0]);
                            }
                            actIndex=0;
                            continue;
                        }
                        data[actIndex++]=character;
                    } else { //stage #2 configuration parameters (coordinates of particles)
                        if (character=='m') {pIndex++; continue;}
                        if (pIndex>=0) {
                            for (int i=0;i<3;i++) {
                                character=fgetc(fileCTL);
                                while (character!=',') {
                                    data[actIndex++]=character;
                                    character=fgetc(fileCTL);
                                } data[actIndex++]=' '; actIndex=0;
                                if (i<2) particles[pIndex].r[i]=strtod(data,NULL);
                                else particles[pIndex].phi=strtod(data,NULL);
                            }
                            particles[pIndex].normR[0]=(boxMatrix[1][1]*particles[pIndex].r[0]-boxMatrix[0][1]*particles[pIndex].r[1])/detBoxMatrix;
                            particles[pIndex].normR[1]=-(boxMatrix[1][0]*particles[pIndex].r[0]-boxMatrix[0][0]*particles[pIndex].r[1])/detBoxMatrix;
                            while (character!=']') character=fgetc(fileCTL); fgetc(fileCTL); //next to read: 'm'
                            if (pIndex>=activeN-1) dataType++;
                        }
                    }
                }
                fclose(fileCTL);
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (startDen: %.7E, startPacFrac: %.7E), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4E, DeltaV: %.4E\n",N,gaps,growing,args[3],args[2],pacFrac,args[0],args[1],(long)args[4],args[8],args[9]);
                //for (int i=0;i<2;i++) for (int j=0;j<2;j++) printf("boxMatrix[%d][%d]=%.17E\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<activeN;i++) printf("%d: %.17E,  %.17E,  %.17E\n",i,particles[i].r[0],particles[i].r[1],particles[i].phi);return 0;
                updateNeighbourList(particles,boxMatrix,volume); //matrix(comment) 13/15
                printf("Checking overlaps in loaded file... "); fflush(stdout);
                if (getEnergyAll(particles,boxMatrix)==1) {
                    printf("Configuration from loaded file contains overlap(s) [energy=1].\n");
                    return 0;
                } else  {printf("done\n"); fflush(stdout);}
            }
        }

        if (skipFirstIteration) {
            printf("Skipping first iteration...\n");
            skipFirstIteration=0;
        } else {
            if (!onlyMath[0]) {
                if (loadedConfiguration) {
                    if (loadedSetStartGenerator) {
                        printf("Setting start position of p-random number generator to position from file...\n");
                        InitMT((unsigned int)args[0]);
                        randomStartStep[0]=args[0];
                    } else {
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                        if (generatorStartPoint==0) {
                            generatorStartPoint=time(0);
                            printf("Setting start position of p-random number generator to actual CPU time...\n");
                        } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                        InitMT((unsigned int)generatorStartPoint);
                        randomStartStep[0]=generatorStartPoint;
                    }
                    randomStartStep[1]=0;
                    if (loadedSetGenerator) {
                        printf("Setting p-random number generator to last position from file...\n");
                        for (double i=0;i<args[1];i++) MTGenerate(randomStartStep);
                    } else printf("Setting p-random number generator to last position from file - DISABLED\n");
                } else {
                    if (generatorStartPoint==0) {
                        generatorStartPoint=time(0);
                        printf("Setting start position of p-random number generator to actual CPU time...\n");
                    } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                    InitMT((unsigned int)generatorStartPoint);
                    randomStartStep[0]=generatorStartPoint;
                    randomStartStep[1]=0;
                }
                printf("Start of equilibration at reduced pressure: %.7E (startDen: %.7E, startPacFrac: %.7E)... (%ld cycles)\n",pressureReduced,rho,pacFrac,cyclesOfEquilibration);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.7E, mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,pressureReduced,multimerN,multimerS,multimerD);
            fflush(stdout);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                fullCycle=activeN+volumeMoveChance,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberV=0, displacedNumberV=0,
                cyclesOfMeasurementBuffer=arg[1]==0?cyclesOfMeasurement:0;
            double deltaRTable[100], deltaRMean=deltaR, deltaVTable[100], deltaVMean=deltaV;
            for (int i=0;i<100;i++) {deltaRTable[i]=deltaRMean; deltaVTable[i]=deltaVMean;}
            int simulationStage=cyclesOfEquilibration>0?0:cyclesOfMeasurementBuffer>0?1:2;  //0-equilibration, 1-measurement, 2-end
            int volumeMove=0, cycleCounter=0, indexScanned=(matrixOfParticlesSize[0]*round(matrixOfParticlesSize[1]*NLinearMod/2.0)-round(matrixOfParticlesSize[0]/2.0))*NLinearMod;
            //assignParticlesToCells(particles,boxMatrix); //matrix 14/15

            char allResultsFileName[200],bufferConfigurations[200],bufferSavedConfigurations[200],bufferOrientations[200],bufferOrientatCorrelFun[200],allOrientationsFileName[200],bufferOrientationsResults[200],allOrientationsResultsFileName[200],bufferPressure[100];
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferOrientations,orientationsFileName); strcpy(allOrientationsFileName,orientationsFileName);
            if (OCFMode) strcpy(bufferOrientatCorrelFun,orientatCorrelFunFileName);
            strcpy(bufferOrientationsResults,orientationsResultsFileName); strcpy(allOrientationsResultsFileName,orientationsResultsFileName);
            sprintf(bufferPressure,"%.4E",pressureReduced);
            strncat(allResultsFileName,"_arg-",6); strncat(allResultsFileName,bufferPressure,100); strncat(allResultsFileName,"_Results.txt",13);
            strncat(bufferConfigurations,"_arg-",6); strncat(bufferConfigurations,bufferPressure,100); strcpy(bufferSavedConfigurations,bufferConfigurations); strncat(bufferConfigurations,".txt",5); strncat(bufferSavedConfigurations,"_transient.txt",15);
            strncat(bufferOrientations,"_arg-",6); strncat(bufferOrientations,bufferPressure,100); strncat(bufferOrientations,".txt",5);
            strncat(allOrientationsFileName,"_arg-",6); strncat(allOrientationsFileName,bufferPressure,100); strncat(allOrientationsFileName,"_allOnt.txt",12);
            if (OCFMode) {
                strncat(bufferOrientatCorrelFun,"_arg-",6); strncat(bufferOrientatCorrelFun,bufferPressure,100); strncat(bufferOrientatCorrelFun,".txt",5);
            }
            strncat(bufferOrientationsResults,"_arg-",6); strncat(bufferOrientationsResults,bufferPressure,100); strncat(bufferOrientationsResults,".txt",5);
            strncat(allOrientationsResultsFileName,"_arg-",6); strncat(allOrientationsResultsFileName,bufferPressure,100); strncat(allOrientationsResultsFileName,"_allOnt.txt",12);

            fileAllResults = fopen(allResultsFileName,"a");
            adjustOrientationsFile(fileOrientations=fopen(bufferOrientations,"rt"),bufferOrientations);
            fileOrientations = fopen(bufferOrientations,"a");
            fileAllOrientations = fopen(allOrientationsFileName,"a");
            if (saveConfigurations) fileSavedConfigurations = fopen(bufferSavedConfigurations,"a");
            if (OCFMode) {
                adjustOrientationsFile(fileOrientatCorrelFun=fopen(bufferOrientatCorrelFun,"rt"),bufferOrientatCorrelFun);
                fileOrientatCorrelFun = fopen(bufferOrientatCorrelFun,"a");
            }
            if (onlyMath[0]) simulationStage=2;

            timeStart=time(0);
            while (simulationStage<2) {
                int randIndex;
                if (volumeMove) {
                    randIndex = (int)(MTRandom0to1(randomStartStep)*activeN);
                    volumeMove=0;
                } else randIndex = (int)(MTRandom0to1(randomStartStep)*fullCycle);
                if (randIndex<activeN) {
                    attemptedNumberR++;
                    if (attemptToDisplaceAParticle(particles,randIndex,boxMatrix))
                        displacedNumberR++;
                } else {
                    volumeMove=1;
                    attemptedNumberV++;
                    if (attemptToChangeVolume(particles,pressure,boxMatrix,volume,xAxisPhi))
                        displacedNumberV++;
                }

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++;

                    if (cycle%intervalSampling==0) {
                        if (simulationStage==0 && cycle>cyclesOfEquilibration) {
                            simulationStage=1;
                            printf("Equilibration finished after: %ld cycles (%ldsec).\n",cyclesOfEquilibration,time(0)-timeStart);
                            fflush(stdout);
                        }
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR,
                               acceptanceRatioV = displacedNumberV/(double)attemptedNumberV;
                        if (cycle%neighUpdatingFrequency==0) updateNeighbourList(particles,boxMatrix,volume);  //matrix(comment) 15/15

                        /////wypisywanie danych czesciej niz normalnie i PRZED zrownowagowaniem
                        /*if (cycle%50==0) {
                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    getParticlesDistanceSquared(particles[i],particles[j],boxMatrix);
                                    double drRoot=sqrt(dr[2]);
                                    int energy;
                                    if (drRoot<maxDistance) {
                                        if (drRoot<minDistance) energy=1;
                                        else {
                                            double gamma=atan(dr[1]/dr[0]),
                                                    aAngle=particles[j].phi-gamma,
                                                    bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (dr[0]>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                                            energy=checkOverlaps(drRoot,aAngle,bAngle);
                                            //energy=checkOverlapsAnalyticalMethodHCM(drRoot,aAngle,bAngle);
                                            if (energy==1) printf("colliding: distance- %.12E, alpha: %.12E, beta: %.12E, analytic: %.12E, i: %d, j: %d\n",drRoot,aAngle,bAngle,getMinimalDistanceAnalyticalMethodForEvenHCM(normalizeAngle(aAngle+C),normalizeAngle(bAngle+C)),i,j);
                                        }
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+(long)args[4]),collidingPairs);
                            else printf("Cycle: %ld\n",(cycle+(long)args[4]));
                            printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                            printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                            printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                            fflush(stdout);
                        }*/
                        /////

                        if (simulationStage==1) {
                            if (timeEquilibration==0) {
                                timeEquilibration=time(0);

                                printf("Checking overlaps in equilibrated configuration... "); fflush(stdout);
                                if (getEnergyAll(particles,boxMatrix)==1) {
                                    printf("Equilibrated configuration contains overlap(s) [energy=1]. ");
                                    char allResultsErrorFileName[200];
                                    strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                                    if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err).\n");
                                    else printf("Error renaming results file (.err).\n");
                                    return 0;
                                } else  {printf("done\n"); fflush(stdout);}
                            }

                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
                                    double drRoot=sqrt(dr[2]);
                                    int energy;
                                    if (drRoot<maxDistance) {
                                        if (drRoot<minDistance) energy=1;   //analyticCheck 1/4
                                        else {                              //
                                            double gamma=atan(dr[1]/dr[0]),
                                                   aAngle=particles[j].phi-gamma,
                                                   bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                //rozważanie, która molekuła jest 'po lewej', a która 'po prawej' (moja metoda-odwraca cząstkę po prawej; analityczna metoda-istotna kolejność)
                                                if (dr[0]>0) bAngle-=C;     //analyticCheck 2/4
                                                else aAngle-=C;             //
                                                /*if (dr[0]<0) {
                                                    double buffer=aAngle; aAngle=bAngle; bAngle=buffer;
                                                }*/
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);   //analyticCheck 3/4
                                            energy=checkOverlaps(drRoot,aAngle,bAngle);                     //
                                            //energy=checkOverlapsAnalyticalMethodHCM(drRoot,aAngle,bAngle);
                                        }  //analyticCheck 4/4
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0)
                                fprintf(fileAllResults,"%.17E\t%.17E\t%.17E\t\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);

                            if (cycle%intervalOrientations==0) {
                                /////skan po konkretnej cząstce (np. dla 224: 14*(16/2)-(14/2)=105, etc.) - w srodku by nie skakala na granicy pudla periodycznego; bezposrednie uzycie w Mathematice (format tablicy)
                                if (cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer)
                                    fprintf(fileOrientations,"{%.12E,%.12E,%.12E}}",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                else fprintf(fileOrientations,"{%.12E,%.12E,%.12E},",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                /////skan po wszystkich cząstkach
                                for (int i=0;i<activeN-1;i++) fprintf(fileAllOrientations,"%.12E,",particles[i].phi);
                                fprintf(fileAllOrientations,"%.12E\n",particles[activeN-1].phi);
                                /////OCF
                                if (OCFMode) {
                                    char bufferText[4096]="", bufferAngle[20];
                                    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                        if (i<particles[i].neighbours[j]) {
                                            strcpy(bufferText,"");
                                            getParticlesDistanceSquared(&particles[i],&particles[particles[i].neighbours[j]],boxMatrix);
                                            double gamma=atan(dr[1]/dr[0]),
                                                   aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                                                   bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (dr[0]>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle+C); bAngle=normalizeAngle(bAngle+C);

                                            strncat(bufferText,"{",2); sprintf(bufferAngle,"%.12E",aAngle); strncat(bufferText,bufferAngle,20);
                                            strncat(bufferText,",",2); sprintf(bufferAngle,"%.12E",bAngle); strncat(bufferText,bufferAngle,20);
                                            if (i>=activeN-2 && cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer) strncat(bufferText,"}}",3); // -2 because it's last particle which has neighbour UNtested (the last one)
                                            else strncat(bufferText,"},",3);
                                            fprintf(fileOrientatCorrelFun,"%s",bufferText);
                                        }
                                    }
                                }
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.17E\t%.17E\t%.17E\t{",(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++)
                                    fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E],",particles[i].r[0],particles[i].r[1],particles[i].phi);
                                fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E]}",particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].phi);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d, ",(cycle+(long)args[4]),collidingPairs);
                                else printf("Cycle: %ld, ",(cycle+(long)args[4]));
                                printf("simulation time: full-%ldsec, measurement-%ldsec\n",time(0)-timeStart,time(0)-timeEquilibration);
                                printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                                printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                fflush(stdout);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                                for (int i=0;i<activeN;i++)
                                    fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                fclose(fileConfigurations);
                            }
                        } else {//dostosowywanie delt w trakcie pomiaru
                            if (acceptanceRatioR>desiredAcceptanceRatioR) deltaR*=1.05; else deltaR*=0.95;
                            if (deltaR>maxDeltaR*multimerS) deltaR=maxDeltaR*multimerS;
                            deltaPhi=deltaR*2.0*sin(C)/multimerS;
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; else deltaV*=0.95;

                            int sampleNumberMod100=(cycle/intervalSampling)%100;
                            updateTableAndGetActualMean(deltaRTable,deltaRMean,sampleNumberMod100,deltaR); deltaR=deltaRMean;
                            updateTableAndGetActualMean(deltaVTable,deltaVMean,sampleNumberMod100,deltaV); deltaV=deltaVMean;
                        }
                        //printf("c: %ld, v*: %.17E, rho: %.17E\n\trR: %.17E -d/aR: %ld/%ld, dR: %.17E\n\trV: %.17E -d/aV: %ld/%ld, dV: %.17E\n\tb00: %.17E, b11: %.17E, b01: %.17E\n",cycle,volume/VcpPerParticle/N,N/volume,acceptanceRatioR,displacedNumberR,attemptedNumberR,deltaR,acceptanceRatioV,displacedNumberV,attemptedNumberV,deltaV,boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);

                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberV=0; displacedNumberV=0;
                    }
                    if (simulationStage==1 && cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer) simulationStage=2;
                }
            }
            fclose(fileAllResults);
            fclose(fileOrientations); fclose(fileAllOrientations);
            if (saveConfigurations) fclose(fileSavedConfigurations); if (OCFMode) fclose(fileOrientatCorrelFun);
            if (timeEquilibration==0) timeEquilibration=time(0);
            printf("Checking overlaps in final configuration... "); fflush(stdout);
            if (!onlyMath[0] && getEnergyAll(particles,boxMatrix)==1) {
                printf("Final configuration contains overlap(s) [energy=1]. ");
                char allResultsErrorFileName[200],configurationErrorFileName[200];
                strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                strcpy(configurationErrorFileName,bufferConfigurations); strncat(configurationErrorFileName,".err",5);
                if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err). ");
                else printf("Error renaming results file (.err). ");
                if (rename(bufferConfigurations,configurationErrorFileName)==0) printf("Configuration file successfully renamed (.err).\n");
                else printf("Error renaming configuration file (.err).\n");
                return 0;
            } else {printf("done\n"); fflush(stdout);}
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            printf("Start of calculation of results...\n");

            //obliczenie liczby linii danych (potrzebne do podziału na zespoły i obliczenia średnich błędów)
            printf("Calculation of data lines... "); fflush(stdout);
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[4096]; double dataLicznik=0; int faultyLines=0, onlyMathLinesBuffer=onlyMath[1];
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) {
                if (fgets(linia,300,fileAllResults)!=NULL && !isLineCorrect(linia)) onlyMathLinesBuffer++;
            }
            while(fgets(linia,300,fileAllResults)!=NULL) {
                if (isLineCorrect(linia)) dataLicznik++; else faultyLines++;
            }
            fclose(fileAllResults);
            printf("done (Found %ld data lines [%d faulty lines occurred].",(long)dataLicznik,faultyLines);
            if ((long)dataLicznik%10>0) printf(" Last %ld lines won't be considered, due to calculations of averages in 10 sets.)\n",(long)dataLicznik%10); else printf(")\n");
            dataLicznik-=(long)dataLicznik%10;

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... "); fflush(stdout);
            double avVolumeSet[10], avBoxMatrixSet[10][3], avRhoSet[10], avPacFracSet[10]; //wyniki dzielone są na 10 zespołów (obliczane są nieskorelowane wzajemnie średnie "lokalne", do obliczenia błędu średniej "globalnej")
            for (int i=0;i<10;i++) {
                avVolumeSet[i]=0; for (int j=0;j<3;j++) avBoxMatrixSet[i][j]=0;
                avRhoSet[i]=0; avPacFracSet[i]=0;
            }
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            double lineCounter=0;
            /*while(fgets(linia,300,fileAllResults)!=NULL) {//OLD-VERSION - large files including: cycles,volume,3xBoxMatrix,rho,pacFrac  1/2
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++; if (actIndex>=300) continue;
                int dataIndex=0; double dataD[6]; while (dataIndex<6) {
                    char data[50]="";
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    dataD[dataIndex++]=strtod(data,NULL);
                }
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    avVolumeSet[setIndex]+=dataD[0]; avBoxMatrixSet[setIndex][0]+=dataD[1];
                    avBoxMatrixSet[setIndex][1]+=dataD[2]; avBoxMatrixSet[setIndex][2]+=dataD[3];
                    avRhoSet[setIndex]+=dataD[4]; avPacFracSet[setIndex]+=dataD[5];
                }
                lineCounter++;
            }
            fclose(fileAllResults);
            double avVolume=0, avBoxMatrix[3]={0,0,0}, avRho=0, avPacFrac=0;
            for (int i=0;i<10;i++) {
                avVolumeSet[i]/=dataLicznik*0.1; avVolume+=avVolumeSet[i];
                avRhoSet[i]/=dataLicznik*0.1; avRho+=avRhoSet[i];
                avPacFracSet[i]/=dataLicznik*0.1; avPacFrac+=avPacFracSet[i];
                for (int j=0;j<3;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
            }*/
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {//NEW-VERSION - small files including only: 3xBoxMatrix  1/2
                sscanf(linia,"%c",linia);
                int actIndex=0;
                int dataIndex=0; double dataD[3]; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10;
                    else dataD[dataIndex++]=strtod(data,NULL);
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    for (int i=0;i<3;i++) avBoxMatrixSet[setIndex][i]+=dataD[i];
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            double avVolume=0, avBoxMatrix[3]={0,0,0}, avRho=0, avPacFrac=0;
            for (int i=0;i<10;i++) {
                for (int j=0;j<3;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
                avVolumeSet[i]=fabs(avBoxMatrixSet[i][0]*avBoxMatrixSet[i][1]-avBoxMatrixSet[i][2]*avBoxMatrixSet[i][2]); avVolume+=avVolumeSet[i];
                avRhoSet[i]=N/avVolumeSet[i]; avRho+=avRhoSet[i];
                avPacFracSet[i]=1.0/VcpPerParticle/avRhoSet[i]; avPacFrac+=avPacFracSet[i];
            }
            avVolume*=0.1; avRho*=0.1; avPacFrac*=0.1; for (int i=0;i<3;i++) avBoxMatrix[i]*=0.1;
            //obliczenie bledow mierzonych wielkosci
            double dAvVolume=0, dAvBoxMatrix[3]={0,0,0}, dAvRho=0, dAvPacFrac=0;
            for (int i=0;i<10;i++) {double epsilon=avVolume-avVolumeSet[i]; dAvVolume+=epsilon*epsilon;} dAvVolume=getAvErrorFromSumEps(dAvVolume,90.0); //10*9 (n(n-1))
            for (int j=0;j<3;j++) {for (int i=0;i<10;i++) {double epsilon=avBoxMatrix[j]-avBoxMatrixSet[i][j]; dAvBoxMatrix[j]+=epsilon*epsilon;} dAvBoxMatrix[j]=getAvErrorFromSumEps(dAvBoxMatrix[j],90.0);}
            for (int i=0;i<10;i++) {double epsilon=avRho-avRhoSet[i]; dAvRho+=epsilon*epsilon;} dAvRho=getAvErrorFromSumEps(dAvRho,90.0);
            for (int i=0;i<10;i++) {double epsilon=avPacFrac-avPacFracSet[i]; dAvPacFrac+=epsilon*epsilon;} dAvPacFrac=getAvErrorFromSumEps(dAvPacFrac,90.0);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... "); fflush(stdout);
            double e1111Set[10], e1122Set[10], e1212Set[10], e2222Set[10], e1112Set[10], e1222Set[10],
                   HxyHyx,HxxHyy,HxxHxy,Hxx2,Hyy2,mod0,mod1;
            for (int i=0;i<10;i++) {e1111Set[i]=0; e1122Set[i]=0; e1212Set[i]=0; e2222Set[i]=0; e1112Set[i]=0; e1222Set[i]=0;}
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            lineCounter=0; int setIndex, oldSetIndex=-1;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                setIndex=(int)(lineCounter/dataLicznik*10.0);
                if (setIndex!=oldSetIndex) {
                    HxyHyx=avBoxMatrixSet[setIndex][2]*avBoxMatrixSet[setIndex][2];
                    HxxHyy=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][1];
                    HxxHxy=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][2];
                    Hxx2=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][0];
                    Hyy2=avBoxMatrixSet[setIndex][1]*avBoxMatrixSet[setIndex][1];
                    mod0=HxyHyx-HxxHyy;
                    mod1=1.0/(2.0*mod0*mod0);
                    oldSetIndex=setIndex;
                }
                sscanf(linia,"%c",linia);
                int actIndex=0;
                /*while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++;  //OLD-VERSION - large files including: cycles,volume,3xBoxMatrix,rho,pacFrac  2/2
                if (actIndex>=300) continue;
                double h11,h22,h12;
                int dataIndex=0; while (dataIndex<4) {
                    char data[50]="";
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    switch (dataIndex++) {
                        case 1: h11=strtod(data,NULL);break;
                        case 2: h22=strtod(data,NULL);break;
                        case 3: h12=strtod(data,NULL);break;
                    }
                }*/
                double h[3],h11,h22,h12;  //NEW-VERSION - small files including only: 3xBoxMatrix  2/2
                int dataIndex=0; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10) {
                        if ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.')) dataIndex=10;
                        else h[dataIndex++]=strtod(data,NULL);
                    }
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    h11=h[0]; h22=h[1]; h12=h[2];
                    double hxyhyx=h12*h12,
                           hxxPhyy=h11+h22,
                           hxx2=h11*h11,
                           hyy2=h22*h22,

                           e11=mod1*(HxyHyx*(hxyhyx-HxyHyx+hyy2)-2.0*(h11*avBoxMatrixSet[setIndex][2]*h12-avBoxMatrixSet[setIndex][0]*HxyHyx+avBoxMatrixSet[setIndex][2]*h12*h22)*avBoxMatrixSet[setIndex][1]+(hxx2-Hxx2+hxyhyx)*Hyy2),
                           e22=mod1*(HxyHyx*(hxx2+hxyhyx-HxyHyx)-2.0*(HxxHxy*h12*hxxPhyy-HxxHyy*HxyHyx)+Hxx2*(hxyhyx+hyy2-Hyy2)),
                           e12=mod1*(-HxxHxy*(hxyhyx+hyy2)+HxxHyy*h12*hxxPhyy+avBoxMatrixSet[setIndex][2]*(h12*avBoxMatrixSet[setIndex][2]*hxxPhyy-(hxx2+hxyhyx)*avBoxMatrixSet[setIndex][1]));

                    e1111Set[setIndex]+=e11*e11;
                    e1122Set[setIndex]+=e11*e22;
                    e1212Set[setIndex]+=e12*e12;
                    e2222Set[setIndex]+=e22*e22;
                    e1112Set[setIndex]+=e11*e12;
                    e1222Set[setIndex]+=e12*e22;
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            double e1111=0, e1122=0, e1212=0, e2222=0, e1112=0, e1222=0;
            for (int i=0;i<10;i++) {
                e1111Set[i]/=dataLicznik*0.1; e1111+=e1111Set[i];
                e1122Set[i]/=dataLicznik*0.1; e1122+=e1122Set[i];
                e1212Set[i]/=dataLicznik*0.1; e1212+=e1212Set[i];
                e2222Set[i]/=dataLicznik*0.1; e2222+=e2222Set[i];
                e1112Set[i]/=dataLicznik*0.1; e1112+=e1112Set[i];
                e1222Set[i]/=dataLicznik*0.1; e1222+=e1222Set[i];
            }
            e1111*=0.1; e1122*=0.1; e1212*=0.1; e2222*=0.1; e1112*=0.1; e1222*=0.1;
            //obliczenie bledow iloczynow elementow tensora odkształceń
            double dE1111=0, dE1122=0, dE1212=0, dE2222=0, dE1112=0, dE1222=0;
            for (int i=0;i<10;i++) {double epsilon=e1111-e1111Set[i]; dE1111+=epsilon*epsilon;} dE1111=getAvErrorFromSumEps(dE1111,90.0); //10*9 (n(n-1))
            for (int i=0;i<10;i++) {double epsilon=e1122-e1122Set[i]; dE1122+=epsilon*epsilon;} dE1122=getAvErrorFromSumEps(dE1122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1212-e1212Set[i]; dE1212+=epsilon*epsilon;} dE1212=getAvErrorFromSumEps(dE1212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e2222-e2222Set[i]; dE2222+=epsilon*epsilon;} dE2222=getAvErrorFromSumEps(dE2222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1112-e1112Set[i]; dE1112+=epsilon*epsilon;} dE1112=getAvErrorFromSumEps(dE1112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1222-e1222Set[i]; dE1222+=epsilon*epsilon;} dE1222=getAvErrorFromSumEps(dE1222,90.0);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... "); fflush(stdout);
            //Eijkl - bezwymiarowe [strain], volume - [sigma^2], kT=1[hard]
            double s1111=e1111*avVolume, dS1111=fabs(e1111*dAvVolume)+fabs(dE1111*avVolume),
                   s1122=e1122*avVolume, dS1122=fabs(e1122*dAvVolume)+fabs(dE1122*avVolume),
                   s1212=e1212*avVolume, dS1212=fabs(e1212*dAvVolume)+fabs(dE1212*avVolume),
                   s2222=e2222*avVolume, dS2222=fabs(e2222*dAvVolume)+fabs(dE2222*avVolume),
                   s1112=e1112*avVolume, dS1112=fabs(e1112*dAvVolume)+fabs(dE1112*avVolume),
                   s1222=e1222*avVolume, dS1222=fabs(e1222*dAvVolume)+fabs(dE1222*avVolume),

                   S11=(s1111+s2222)*0.5,
                   S66=4.0*s1212,
                   avNu=-s1122/S11, dAvNu=fabs(dS1122/S11)+fabs((dS1111+dS2222)*0.5*s1122/S11/S11), //nu obliczane z Sxxyy
                   avNu2=S66/S11*0.5-1, dAvNu2=fabs(0.5/S11*4.0*dS1212)+fabs(0.5*S66/S11/S11*(dS1111+dS2222)*0.5), //nu obliczane z Sxyxy (inny rodzaj scinania, przy izotropowych ukladach powinno byc tyle samo co nu1)

                   //stary sposob obliczania nu (nu_x, nu_y i srednia)
                   /*nu2211_1111=-s1122/s1111, dNu2211_1111=fabs(dS1122/s1111)+fabs(dS1111*s1122/s1111/s1111),
                   nu1122_2222=-s1122/s2222, dNu1122_2222=fabs(dS1122/s2222)+fabs(dS2222*s1122/s2222/s2222),
                   avNu=(nu2211_1111+nu1122_2222)/2.0, dAvNu=(dNu2211_1111+dNu1122_2222)/2.0,*/

                   //\lambdaReduced=\lambda*\sigma^2/kT, kT=1[hard]
                   l11=multimerS*multimerS/(8.0*(s1111+s1122)), dL11=multimerS*multimerS*(fabs(dS1111)+fabs(dS1122))/(8.0*fabs(s1111+s1122)*fabs(s1111+s1122)),
                   l12=multimerS*multimerS/(8.0*(s2222+s1122)), dL12=multimerS*multimerS*(fabs(dS2222)+fabs(dS1122))/(8.0*fabs(s2222+s1122)*fabs(s2222+s1122)),

                   B1=4.0*l11, dB1=4.0*dL11,
                   B2=4.0*l12, dB2=4.0*dL12,
                   avB=(B1+B2)/2.0, dAvB=(dB1+dB2)/2.0,

                   my1=(B1-B1*avNu)/(1.0+avNu), dMy1=fabs(dB1*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-B1/(1.0+avNu)-(B1-B1*avNu)/(1.0+avNu)/(1.0+avNu))),
                   my2=(B2-B2*avNu)/(1.0+avNu), dMy2=fabs(dB2*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-B2/(1.0+avNu)-(B2-B2*avNu)/(1.0+avNu)/(1.0+avNu))),
                   avMy=(my1+my2)/2.0, dAvMy=(dMy1+dMy2)/2.0,

                   avE=4.0*avB*avMy/(avB+avMy), dAvE=4.0*(fabs(avB*avB*dAvMy)+fabs(dAvB*avMy*avMy))/fabs(avB+avMy)/fabs(avB+avMy);
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla 1 czastki
            printf("Creation of 1-particle orientation file for Origin... "); fflush(stdout);
            double componentCounter=0, averageCos6PhiOne=0, ODFMaxOne=0;
            fileOrientations=fopen(bufferOrientations,"rt");
            double ODF_1P[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_1P[i]=0; //Orientational Distribution Function (1 Particle)
            int licznik=1,character;
            while ((character=fgetc(fileOrientations))!=EOF) {
                if (character==',') licznik++;
                if (licznik==3) {
                    licznik=0;
                    char data[50]=""; int actIndex=0;
                    while (true) {
                        character=fgetc(fileOrientations);
                        if (character!='}' && actIndex<50) data[actIndex++]=character;
                        else {
                            double angle = normalizeAngle(strtod(data,NULL)+C);
                            averageCos6PhiOne+=cos(6.0*angle); componentCounter++;
                            int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                            ODF_1P[index]++;
                            break;
                        }
                    }
                }
            }
            fclose(fileOrientations);
            double dPhi=2.0*C/((double)(ODFLength-1.0));
            double suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_1P[i]; for (int i=0;i<ODFLength;i++) ODF_1P[i]/=suma*dPhi;
            averageCos6PhiOne/=componentCounter; for (int i=0;i<ODFLength;i++) if (ODFMaxOne<ODF_1P[i]) ODFMaxOne=ODF_1P[i];
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla wszystkich czastek
            printf("Creation of ALL-particle orientation file for Origin... "); fflush(stdout);
            componentCounter=0; double averageCos6PhiAll=0, ODFMaxAll=0, averagePhiAll=0;
            fileAllOrientations = fopen(allOrientationsFileName,"rt");
            double ODF_AllP[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]=0; //Orientational Distribution Function (All Particles)
            char data[50]=""; int actIndex=0;
            while ((character=fgetc(fileAllOrientations))!=EOF) {
                if (character!=',' && character!='\n') data[actIndex++]=character;
                else {
                    data[actIndex++]=' '; actIndex=0;
                    double angle = normalizeAngle(strtod(data,NULL)+C);
                    averagePhiAll+=angle;
                    averageCos6PhiAll+=cos(6.0*angle); componentCounter++;
                    int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                    ODF_AllP[index]++;
                }
            }
            fclose(fileAllOrientations);
            suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_AllP[i]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]/=suma*dPhi;
            averagePhiAll/=componentCounter; averageCos6PhiAll/=componentCounter;
            int maxODFAllIndex;
            for (int i=0;i<ODFLength;i++) if (ODFMaxAll<ODF_AllP[i]) {
                ODFMaxAll=ODF_AllP[i];
                maxODFAllIndex=i;
            }
            printf("done\n");

            //analiza konfiguracji przejsciowych (dPhi z konfiguracji na konfiguracje)
            double avAbsDPhi=0;
            if (saveConfigurations) {
                printf("Transient configurations analysis... "); fflush(stdout);
                fileSavedConfigurations = fopen(bufferSavedConfigurations,"rt");
                double prevCfg[activeN][3]; bool undefinedPrevCfg=true;
                int CfgQuantity=0,CfgFileQuantity=0;
                while ((character=fgetc(fileSavedConfigurations))!=EOF) {
                    if (character=='n') {//text: newCfgFile (after merging) [it's NOT \n -> it's n, first char of 'new']
                        undefinedPrevCfg=true;
                        while (fgetc(fileSavedConfigurations)!='\n');
                        continue;
                    } else CfgQuantity++;
                    while (fgetc(fileSavedConfigurations)!='[');
                    for (int i=0;i<activeN;i++) {
                        for (int j=0;j<3;j++) {
                            strcpy(data,""); actIndex=0;
                            while ((character=fgetc(fileSavedConfigurations))!=',' && character!=']') data[actIndex++]=character;
                            data[actIndex++]=' ';
                            if (j==2) {
                                for (int k=0;k<3;k++) fgetc(fileSavedConfigurations); //skip ",m["
                                if (!undefinedPrevCfg) {
                                    avAbsDPhi+=fabs(strtod(data,NULL)-prevCfg[i][j]);
                                }
                            }
                            prevCfg[i][j]=strtod(data,NULL);
                        }
                    }
                    if (undefinedPrevCfg) {undefinedPrevCfg=false; CfgFileQuantity++;}
                }
                fclose(fileSavedConfigurations);
                avAbsDPhi/=(double)activeN*(CfgQuantity-CfgFileQuantity);
                printf("done\n");
            }

            long timeEndMath=time(0);




/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... "); fflush(stdout);
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) {
                fileConfigurations = fopen(bufferConfigurations,"w");
                fileConfigurationsList = fopen(configurationsListFileName,"a");
            }
            fileOrientationsResults = fopen(bufferOrientationsResults,"w");
            fileAllOrientationsResults = fopen(allOrientationsResultsFileName,"w");

            if (saveConfigurations) {
                fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%ld\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avVolume,dAvVolume,avBoxMatrix[0],dAvBoxMatrix[0],avBoxMatrix[1],dAvBoxMatrix[1],avBoxMatrix[2],dAvBoxMatrix[2],avRho,dAvRho,avPacFrac,dAvPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,s1112,dS1112,s1222,dS1222,avNu,dAvNu,avNu2,dAvNu2,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,savedConfigurationsInt,avAbsDPhi);
                fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%ld\t%.17E\n",pressureReduced,avPacFrac,avNu,avNu2,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,avB,avMy,avE,savedConfigurationsInt,avAbsDPhi);
            } else {
                fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avVolume,dAvVolume,avBoxMatrix[0],dAvBoxMatrix[0],avBoxMatrix[1],dAvBoxMatrix[1],avBoxMatrix[2],dAvBoxMatrix[2],avRho,dAvRho,avPacFrac,dAvPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,s1112,dS1112,s1222,dS1222,avNu,dAvNu,avNu2,dAvNu2,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll);
                fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",pressureReduced,avPacFrac,avNu,avNu2,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,avB,avMy,avE);
            }

            if (!onlyMath[0]) {
                rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                fprintf(fileConfigurationsList,"multimers[x_,y_,kI_]:={");
                for (int i=0;i<activeN;i++) {
                    fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                    fprintf(fileConfigurationsList,"m[%.12E+x,%.12E+y,%.12E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                }
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                fprintf(fileConfigurationsList,"{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[If[x==0 && y==0,0.4,0]],Green,Disk[{%.12E,%.12E},%.12E]}};\nconfigurationsList=Append[configurationsList,g[%.12E,multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                        boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,pacFrac,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
            }

            for (int i=0;i<ODFLength;i++) {
                fprintf(fileOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_1P[i]);
                fprintf(fileAllOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_AllP[i]);
            }

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) {
                fclose(fileConfigurations); fclose(fileConfigurationsList);
            }
            fclose(fileOrientationsResults); fclose(fileAllOrientationsResults);
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        getNextArgument(arg,true);
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        args[4]=0;
        loadedConfiguration=0;
        generatorStartPoint=0;
    }
    printf("\nTime for equilibrations: %ldsec, time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n"); fflush(stdout);
}
