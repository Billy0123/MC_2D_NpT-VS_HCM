#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "MTGenerator.h"

int N,gaps,activeN,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,multimerN,countCollidingPairs,ODFLength,skipFirstIteration,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    lambdaSetIndex,onlyMath[2]={0,0};
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,intervalOrientations;
double maxDeltaR,desiredAcceptanceRatioR,
       minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg,deltaR, //*sigma(=1)
       multimerS,multimerD,randomStartStep[2],normalizingDenominator,
       neighRadius,neighRadius2,neighSafeDistance,multiplyFactor,
       iterationTable[1000][2],pi=3.1415926535898,
       lambda[2], displacementsSum[2], avPhi, massCenter[2][2], massCenterLattice[2][2];  //massCenter[0]:r, massCenter[1]:normR
double L,C,ROkreguOpisanego,absoluteMinimum,absoluteMinimum2,minDistance,maxDistance,VcpPerParticle;
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferMN[20],bufferMS[20],bufferMD[20],bufferFolderIndex[5],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     orientationsFileName[200]="Orientations",
     orientationsResultsFileName[200]="OrientatRes",
     configurationsListFileName[200]="ConfigurationsList.txt",
     loadConfigurationsFileName[200]="Configurations",
     loadedJOBID[50]="j-none";

/////////////////  PARTICLE functions{
typedef struct particle {
    double r[2], normR[2];  //x,y
    double phi;   //kąt mierzony od kierunku x
    int neighbours[50], neighCounter;
} particle;

void updateNeighbourList (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
               normalizedRY=particles[i].normR[1]-particles[j].normR[1],
               rx=particles[i].r[0]-particles[j].r[0],
               ry=particles[i].r[1]-particles[j].r[1];
        rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
        ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
        double r2=rx*rx+ry*ry;
        if (r2<neighRadius2) {
            particles[i].neighbours[particles[i].neighCounter++]=j;
            particles[j].neighbours[particles[j].neighCounter++]=i;
        }
    }
}

double normalizeAngle (double phi) {//utrzymywanie phi w zakresie (-30st,30st) = (-Pi/n,Pi/n)
    int change;
    do {
        change=0;
        if (phi<-C) {
            phi+=2*C;
            change=1;
        } else if (phi>C) {
            phi-=2*C;
            change=1;
        }
    } while (change==1);
    return phi;
}

double get2DiscsDistance (int aNum, int bNum, double dr, double aAngle, double bAngle) {
    double discAngle = aAngle+aNum*2*C;
    double aNumDiscX = cos(discAngle)*ROkreguOpisanego,
           aNumDiscY = sin(discAngle)*ROkreguOpisanego;
    discAngle = bAngle+bNum*2*C;
    double bNumDiscX = dr+cos(discAngle)*ROkreguOpisanego,
           bNumDiscY = sin(discAngle)*ROkreguOpisanego;
    double xDistance = aNumDiscX-bNumDiscX, yDistance = aNumDiscY-bNumDiscY;
    return sqrt(xDistance*xDistance+yDistance*yDistance);
}

int checkOverlaps (double dr, double aAngle, double bAngle) {
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

double minimalDistanceForEvenHCM (int indeksPowierzchni, double aAngle, double bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-C,-C,C,C};
    double angleA=aAngle+a[indeksPowierzchni],angleB=bAngle+b[indeksPowierzchni],
           buffer=sin(angleA)+sin(angleB);
    return (ROkreguOpisanego*(cos(angleA)+cos(angleB))+sqrt(multimerD*multimerD-ROkreguOpisanego*ROkreguOpisanego*buffer*buffer));
}

double getMinimalDistanceAnalyticalMethodForEvenHCM (double aAngle, double bAngle) {
    if (aAngle<-absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13=-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceForEvenHCM(0,aAngle,bAngle);
        else return minimalDistanceForEvenHCM(2,aAngle,bAngle);
    } else if (aAngle<0) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve14=aAngle,
               intersectionCurve34=-C-ARCSIN;
        if (bAngle>intersectionCurve14) return minimalDistanceForEvenHCM(0,aAngle,bAngle);
        else if (bAngle>intersectionCurve34) return minimalDistanceForEvenHCM(3,aAngle,bAngle);
        else return minimalDistanceForEvenHCM(2,aAngle,bAngle);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve14=aAngle,
               intersectionCurve12=C-ARCSIN;
        if (bAngle>intersectionCurve12) return minimalDistanceForEvenHCM(1,aAngle,bAngle);
        else if (bAngle>intersectionCurve14) return minimalDistanceForEvenHCM(0,aAngle,bAngle);
        else return minimalDistanceForEvenHCM(3,aAngle,bAngle);
    } else {
        double ARCSIN=asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24=-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceForEvenHCM(1,aAngle,bAngle);
        else return minimalDistanceForEvenHCM(3,aAngle,bAngle);
    }
}

int checkOverlapsAnalyticalMethodForEvenHCM (double dr, double aAngle, double bAngle) {
    int overlap=0;

    aAngle=normalizeAngle(aAngle+C);
    bAngle=normalizeAngle(bAngle+C);

    if (dr<getMinimalDistanceAnalyticalMethodForEvenHCM(aAngle,bAngle)) overlap=1;
    return overlap;
}

int createRandomGaps (particle *particles, double boxMatrix[2][2]) {
    printf("Creating %d random gaps... ",gaps);
    InitRandomMT();
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTGenerate(randomStartStep)%1000000/1000000.0*N);
            for (int j=0;j<i;j++) {
                double normalizedRX=particles[gapsIndexes[i]].normR[0]-particles[gapsIndexes[j]].normR[0],
                       normalizedRY=particles[gapsIndexes[i]].normR[1]-particles[gapsIndexes[j]].normR[1],
                       rx=particles[gapsIndexes[i]].r[0]-particles[gapsIndexes[j]].r[0],
                       ry=particles[gapsIndexes[i]].r[1]-particles[gapsIndexes[j]].r[1];
                rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                double r2=rx*rx+ry*ry, dr=sqrt(r2);
                if (dr<neighRadius) {
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

        //massCenter change due gaps existence
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) massCenter[i][j]=0;
        for (int i=0;i<activeN;i++) {
            massCenter[0][0]+=particles[i].r[0]; massCenter[0][1]+=particles[i].r[1];
            massCenter[1][0]+=particles[i].normR[0]; massCenter[1][1]+=particles[i].normR[1];
        }
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) massCenter[i][j]/=(double)activeN;

        printf("done\n");
        return 1;
    }
}

int initPositions (particle *particles, double boxMatrix[2][2], bool lattice) {//dziala dla heksamerow (dla penta juz nie)
    double modX, modY;
    switch (N) {
        case 56:{modX=7.0;modY=8.0;}break;
        case 224:{modX=14.0;modY=16.0;}break;
        case 504:{modX=21.0;modY=24.0;}break;
        case 780:{modX=26.0;modY=30.0;}break;
        case 896:{modX=28.0;modY=32.0;}break;
        case 3120:{modX=52.0;modY=60.0;}break;
    }
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) if (lattice) massCenterLattice[i][j]=0; else massCenter[i][j]=0;
    double xInterval=boxMatrix[0][0]/modX, yInterval=boxMatrix[1][1]/modY,
           xPosition=0, yPosition=0;
    int licznik=1;
    for (int i=0;i<N;i++) {
        particles[i].r[0]=xPosition; particles[i].r[1]=yPosition;
        particles[i].normR[0]=(-boxMatrix[1][1]*particles[i].r[0]+boxMatrix[1][0]*particles[i].r[1])/normalizingDenominator;
        particles[i].normR[1]=(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/normalizingDenominator;
        xPosition+=xInterval;
        if (xPosition+minDistance/2.0>boxMatrix[0][0]) {
            xPosition=(licznik++%2)*xInterval/2.0;
            yPosition+=yInterval;
        }
        particles[i].phi=absoluteMinimum2;
        if (lattice) {
            massCenterLattice[0][0]+=particles[i].r[0]; massCenterLattice[0][1]+=particles[i].r[1];
            massCenterLattice[1][0]+=particles[i].normR[0]; massCenterLattice[1][1]+=particles[i].normR[1];
        } else {
            massCenter[0][0]+=particles[i].r[0]; massCenter[0][1]+=particles[i].r[1];
            massCenter[1][0]+=particles[i].normR[0]; massCenter[1][1]+=particles[i].normR[1];
        }
    }
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) if (lattice) massCenterLattice[i][j]/=(double)N; else massCenter[i][j]/=(double)N;
    if (!lattice && gaps>0) return createRandomGaps(particles,boxMatrix);
    else return 1;
}

void adjustAngles (particle *particles, double boxMatrix[2][2]) {
    int tryNumber=-1,collidingPairs=0;
    printf("Angle adjusting... ");
    do {
        tryNumber++;

        collidingPairs=0;
        for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
            double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
                   normalizedRY=particles[i].normR[1]-particles[j].normR[1],
                   rx=particles[i].r[0]-particles[j].r[0],
                   ry=particles[i].r[1]-particles[j].r[1];
            rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
            ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
            double r2=rx*rx+ry*ry, dr=sqrt(r2);
            int energy;
            if (dr<maxDistance) {
                if (dr<minDistance) energy=1;
                else {
                    double gamma=atan(ry/rx),
                            aAngle=particles[j].phi-gamma,
                            bAngle=particles[i].phi-gamma;
                    if (multimerN%2!=0) {
                        if (rx>0) bAngle-=C;
                        else aAngle-=C;
                    }
                    aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                    energy=checkOverlaps(dr,aAngle,bAngle);
                }
                if (energy==1) {
                    collidingPairs++;
                    for (int k=0;k<activeN;k++) {
                        particles[k].phi+=0.000001;
                    }
                    i=activeN; j=activeN; break;
                }
            }
        }
    } while (collidingPairs>0);
    printf("Adjusted after: %d approach.\n",tryNumber);
}

void checkSinglePeriodicBoundaryConditions (particle *particle, double boxMatrix[2][2]) {
    for (int j=0;j<2;j++) {
        int change;
        do {
            change=0;
            if (particle->normR[j]<0) {
                particle->normR[j]++;
                particle->r[0]+=boxMatrix[0][j];
                particle->r[1]+=boxMatrix[1][j];
                change=1;
            } else if (particle->normR[j]>=1) {
                particle->normR[j]--;
                particle->r[0]-=boxMatrix[0][j];
                particle->r[1]-=boxMatrix[1][j];
                change=1;
            }
        } while (change==1);
    }
    particle->phi=normalizeAngle(particle->phi);
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

void computeDisplacementsSum (particle *particles, particle *latticeNodes, double boxMatrix[2][2]) { //zalozenie - brak overlapow
    double massCenterDelta[2][2];
    for (int i=0;i<2;i++) {
        displacementsSum[i]=0;
        for (int j=0;j<2;j++) massCenterDelta[i][j]=massCenter[i][j]-massCenterLattice[i][j];
    }
    for (int i=0;i<activeN;i++) {
        double normalizedRX=particles[i].normR[0]-latticeNodes[i].normR[0]-massCenterDelta[1][0],
               normalizedRY=particles[i].normR[1]-latticeNodes[i].normR[1]-massCenterDelta[1][1],
               rx=particles[i].r[0]-latticeNodes[i].r[0]-massCenterDelta[0][0],
               ry=particles[i].r[1]-latticeNodes[i].r[1]-massCenterDelta[0][1];
        rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
        ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
        displacementsSum[0]+=rx*rx+ry*ry;
        double sinDeltaPhi=sin(particles[i].phi-avPhi);
        displacementsSum[1]+=sinDeltaPhi*sinDeltaPhi;
    }
}

int getOverlapsAndDisplacements (particle *particles, particle *latticeNodes, int index, double boxMatrix[2][2]) {
    int overlap=0;
    for (int i=0;i<particles[index].neighCounter;i++) {
        double normalizedRX=particles[particles[index].neighbours[i]].normR[0]-particles[index].normR[0],
               normalizedRY=particles[particles[index].neighbours[i]].normR[1]-particles[index].normR[1],
               rx=particles[particles[index].neighbours[i]].r[0]-particles[index].r[0],
               ry=particles[particles[index].neighbours[i]].r[1]-particles[index].r[1];
        rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
        ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
        double r2=rx*rx+ry*ry, dr=sqrt(r2);
        if (dr<maxDistance) {
            if (dr<minDistance) overlap=1; //analyticMethodForEvenHCM 1/6
            else {
                double gamma=atan(ry/rx),
                       aAngle=particles[index].phi-gamma,
                       bAngle=particles[particles[index].neighbours[i]].phi-gamma;
                if (multimerN%2!=0) {
                    if (rx>0) bAngle-=C;
                    else aAngle-=C;
                }
                aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle); //analyticMethodForEvenHCM 2/6
                overlap=checkOverlaps(dr,aAngle,bAngle);
                //energies.overlap=checkOverlapsAnalyticalMethodForEvenHCM(dr,aAngle,bAngle);
            } //analyticMethodForEvenHCM 3/6
            if (overlap==1) i=particles[index].neighCounter;
        }
    }
    if (overlap==0) computeDisplacementsSum(particles,latticeNodes,boxMatrix);
    return overlap;
}

int attemptToDisplaceAParticle (particle *particles, particle *latticeNodes, int index, double boxMatrix[2][2]) {
    int result=1;
    double oldR[2]={particles[index].r[0],particles[index].r[1]},
           oldNormR[2]={particles[index].normR[0],particles[index].normR[1]},
           oldPhi=particles[index].phi,
           oldMassCenter[2][2]={{massCenter[0][0],massCenter[0][1]},{massCenter[1][0],massCenter[1][1]}},
           oldDisplacementsSum[2]={displacementsSum[0],displacementsSum[1]};
    for (int i=0;i<2;i++) particles[index].r[i]+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaR;
    particles[index].normR[0]=(-boxMatrix[1][1]*particles[index].r[0]+boxMatrix[1][0]*particles[index].r[1])/normalizingDenominator;
    particles[index].normR[1]=(boxMatrix[1][0]*particles[index].r[0]-boxMatrix[0][0]*particles[index].r[1])/normalizingDenominator;
    particles[index].phi+=(MTGenerate(randomStartStep)%1000000/1000000.0-0.5)*deltaR;
    for (int i=0;i<2;i++) {
        massCenter[0][i]+=(particles[index].r[i]-oldR[i])/(double)activeN;
        massCenter[1][i]+=(particles[index].normR[i]-oldNormR[i])/(double)activeN;
    }
    int overlap=getOverlapsAndDisplacements(particles,latticeNodes,index,boxMatrix);
    if (overlap==0 && MTGenerate(randomStartStep)%1000000/1000000.0<exp((oldDisplacementsSum[0]-displacementsSum[0])*lambda[0]+(oldDisplacementsSum[1]-displacementsSum[1])*lambda[1])) checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix);
    else {
        for (int i=0;i<2;i++) {
            particles[index].r[i]=oldR[i];
            particles[index].normR[i]=oldNormR[i];
            for (int j=0;j<2;j++) massCenter[i][j]=oldMassCenter[i][j];
            displacementsSum[i]=oldDisplacementsSum[i];
        }
        particles[index].phi=oldPhi;
        result=0;
    }
    return result;
}

/////////////////  } PARTICLE functions

int createIterationTable () {
    char startArguments[50]; FILE *fileStartArguments = fopen("startArguments.txt","rt");
    if (fileStartArguments==NULL) {printf("Missing file: startArguments.txt\n"); return 1;}
    while (fgets(startArguments,50,fileStartArguments)!=NULL) {
        sscanf(startArguments,"%c",startArguments); char *pEnd;
        iterationTable[fileIterateIterationsNumber][0]=strtod(startArguments,&pEnd);
        iterationTable[fileIterateIterationsNumber++][1]=strtod(pEnd,NULL);
    }
    fclose(fileStartArguments);
    return 0;
}

int getLambdaSet () {
    char lambdaSet[50]; FILE *fileLambdaSets = fopen("lambdaSets.txt","rt");
    if (fileLambdaSets==NULL) {printf("Missing file: lambdaSets.txt\n"); return 1;}
    int licznik=0;
    while (fgets(lambdaSet,50,fileLambdaSets)!=NULL) {
        if (licznik++==lambdaSetIndex) {
            sscanf(lambdaSet,"%c",lambdaSet); char *pEnd;
            lambda[0]=strtod(lambdaSet,&pEnd);
            lambda[1]=strtod(pEnd,NULL);
            break;
        }
    }
    fclose(fileLambdaSets);
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

double getNextArgument (double prevArg, bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            prevArg=growing?iterationTable[actIteration][0]:iterationTable[fileIterateIterationsNumber-1-actIteration][0];
            avPhi=growing?iterationTable[actIteration][1]:iterationTable[fileIterateIterationsNumber-1-actIteration][1];
            if (growing) minArg=iterationTable[actIteration][0];
            else maxArg=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg*10000)>=round(intervalMin[i]*10000) && round(prevArg*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg*10000)+round(intervalDelta[i]*10000);
                prevArg=newArg/10000.0;
                break;
            }
        if (round(prevArg*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg*10000)>round(intervalMin[i]*10000) && round(prevArg*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg*10000)-round(intervalDelta[i]*10000);
                prevArg=newArg/10000.0;
                break;
            }
        if (round(prevArg*10000)<round(minArg*10000)) growing=-1;
    }
    return prevArg;
}

double getAvErrorFromSumEps (double sum, double denominator) {
    return sqrt(sum/denominator);
}

bool A (char txt[10], bool result) {
    printf("%s\n",txt);
    return result;
}


int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[300];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,300,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20];
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++];
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:multimerN=strtol(data,NULL,10);break;
                case 4:multimerS=strtod(data,NULL);break;
                case 5:multimerD=strtod(data,NULL);break;
                case 6:growing=strtol(data,NULL,10);break;
                case 7:loadedConfiguration=strtol(data,NULL,10);break;
                case 8:loadedArg=strtod(data,NULL);break;
                case 9:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 10:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 11:loadedSetGenerator=strtol(data,NULL,10);break;
                case 12:iterationsNumber=strtol(data,NULL,10);break;
                case 13:countCollidingPairs=strtol(data,NULL,10);break;
                case 14:intervalSampling=strtol(data,NULL,10);break;
                case 15:intervalOutput=strtol(data,NULL,10);break;
                case 16:ODFLength=strtol(data,NULL,10);break;
                case 17:intervalOrientations=strtol(data,NULL,10);break;
                case 18:skipFirstIteration=strtol(data,NULL,10);break;
                case 19:useSpecificDirectory=strtol(data,NULL,10);break;
                case 20:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 21:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 22:intervalResults=strtol(data,NULL,10);break;
                case 23:maxDeltaR=strtod(data,NULL);break;
                case 24:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 25:lambdaSetIndex=strtol(data,NULL,10);break;
                case 26:useFileToIterate=strtol(data,NULL,10);break;
                case 27:avPhi=strtod(data,NULL);break;
                case 28:minArg=strtod(data,NULL);break;
                case 29:maxArg=strtod(data,NULL);break;
                case 30:multiplyArgument=strtol(data,NULL,10);break;
                case 31:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-32)%3) {
                        case 0: intervalMin[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 1: intervalMax[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 2: intervalDelta[intervalLicznik++/3]=strtod(data,NULL);break;
                    }
                    break;
            }
            dataIndex++;
            for (int i=0;i<20;i++) data[i]=' ';
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
                if (argumentsNumber==11) {
                    useFileToIterate=0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (growing) minArg=strtod(arguments[3],NULL);
                    else maxArg=strtod(arguments[3],NULL);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    iterationsNumber=strtol(arguments[9],NULL,10);
                    useSpecificDirectory=strtol(arguments[10],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==13) {
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
                    lambdaSetIndex=strtol(arguments[12],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==13) {
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
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==14) {
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
                    lambdaSetIndex=strtol(arguments[13],NULL,10);
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
                    skipFirstIteration=0;
                    lambdaSetIndex=strtol(arguments[12],NULL,10);
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startArg (minArg/maxArg depending on growing), N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, lambdaSetIndex\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, lambdaSetIndex\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, lambdaSetIndex\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        minArg=iterationTable[0][0];
        maxArg=iterationTable[fileIterateIterationsNumber-1][0];
        avPhi=growing?iterationTable[0][1]:iterationTable[fileIterateIterationsNumber-1][1];
    }
    startArg=growing?minArg:maxArg;
    deltaR=maxDeltaR*multimerS;
    for (int i=0;i<pointNumber;i++) startArg=getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg;
    activeN=N-gaps;
    if (getLambdaSet()) return 0;


    //stale wynikajace z zadanych parametrow multimerow
    L=multimerS/multimerD;
    C=pi/(double)multimerN;
    ROkreguOpisanego=multimerS/(2.0*sin(C));
    absoluteMinimum=atan(L/(2.0*L/tan(C)+sqrt(4.0-L*L)));
    absoluteMinimum2=C-absoluteMinimum;
    minDistance=getMinimalDistanceAnalyticalMethodForEvenHCM(absoluteMinimum,absoluteMinimum); //tylko dla parzystych HCM
    maxDistance=ROkreguOpisanego*2+multimerD;
    neighRadius=1.4*/*dla 0.51 i fazy chiralnej 1.3 powinno wystarczyc*1.3*/maxDistance; neighRadius2=neighRadius*neighRadius; neighSafeDistance=neighRadius-maxDistance;
    VcpPerParticle=minDistance*minDistance*sqrt(3)/2.0;  //przynajmniej dla heksamerow
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
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.3f",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,"_l-",4); sprintf(buffer,"%d",lambdaSetIndex); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(orientationsFileName,JOBID,true);
    addAppendix(orientationsResultsFileName,JOBID,true);
    addAppendix(configurationsListFileName,JOBID,false);

    particle particles[N], latticeNodes[N];
    long arg5=0; double arg1, arg2, arg3, arg4, arg6, arg7, arg8, arg9;

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileOrientations, *fileConfigurationsList, *fileAllResults, *fileAllOrientations, *fileOrientationsResults, *fileAllOrientationsResults;
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        fprintf(fileResults,"Cycles\tVolume\tBoxMatrix[0][0]\tBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tRho\tV/V_cp\tlambda[0]\tavDisplacements[0]\tdAvDisplacements[0]\tlambda[1]\tavDisplacements[1]\tdAvDisplacements[1]\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        fprintf(fileExcelResults,"Volume\tRho\tV/V_cp\tlambda[0]\tavDisplacements[0]\tlambda[1]\tavDisplacements[1]\tODFMax_All\t<cos(6Phi)>_All\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg=startArg, oldVolume, oldBoxMatrix[2][2];
    while (growing>=0) { //OBECNA WERSJA NIE NADAJE SIĘ DO SEKWENCJI (NIE MA PROCEDUR ODPOWIADAJĄCYCH ZA SPRĘŻANIE/ROZPRĘŻANIE UKŁADU)
        double volume=(startArg==arg)?(growing?N/(1.0/VcpPerParticle/minArg):N/(1.0/VcpPerParticle/maxArg)):oldVolume, rho=N/volume, pacFrac=1.0/VcpPerParticle/rho, //(pacFrac==arg)
                boxMatrix[2][2],boxSizeCoefficient;
        if (startArg==arg) {
            if (N==56 || N==224 || N==504 || N==896) {
                boxSizeCoefficient=sqrt(volume/(28*sqrt(3)));
                boxMatrix[0][0]=7*boxSizeCoefficient;
                boxMatrix[1][1]=8*sqrt(3/4.0)*boxSizeCoefficient;
            } else if (N==780 || N==3120) {
                boxSizeCoefficient=sqrt(volume/(97.5*sqrt(3)));
                boxMatrix[0][0]=13*boxSizeCoefficient;
                boxMatrix[1][1]=15*sqrt(3/4.0)*boxSizeCoefficient;
            }
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0;
        } else for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
        normalizingDenominator=pow(boxMatrix[1][0],2)-boxMatrix[0][0]*boxMatrix[1][1];

        if (!onlyMath[0]) {
            if (arg==startArg && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartDen: %.4f, startPacFrac: %.4f, lambda[0]: %.4f, lambda[1]: %.4f, avPhi: %.4f, mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,rho,pacFrac,lambda[0],lambda[1],avPhi,multimerN,multimerS,multimerD);
                if (!initPositions(particles,boxMatrix,false)) return 0;
                adjustAngles(particles,boxMatrix);
                updateNeighbourList(particles,boxMatrix);
            } else if (loadedConfiguration) {
                char configurations[400+110*activeN];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                int numerLinii=1;
                while(fgets(configurations,400+110*activeN,fileCTL)!=NULL && numerLinii++<=3) sscanf(configurations,"%c",configurations);
                fclose(fileCTL);

                char *pEnd;
                arg1=strtod(configurations,&pEnd);  //RandStart
                arg2=strtod(pEnd,&pEnd);            //RandStep
                arg3=strtod(pEnd,&pEnd);            //Density
                arg4=strtod(pEnd,&pEnd);            //PacFrac
                arg5=strtol(pEnd,&pEnd,10);         //Cycles
                arg6=strtod(pEnd,&pEnd);            //boxMatrix[0][0]
                arg7=strtod(pEnd,&pEnd);            //boxMatrix[1][1]
                arg8=strtod(pEnd,&pEnd);            //boxMatrix[0][1]=boxMatrix[1][0] (symetryczna macierz)
                arg9=strtod(pEnd,&pEnd);            //deltaR

                int actIndex=1;
                for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
                    int licznik=0;
                    char massCenterCoordinate[50]="";
                    while (pEnd[actIndex]!=',' && pEnd[actIndex]!=' ') massCenterCoordinate[licznik++]=pEnd[actIndex++];
                    actIndex++;
                    massCenter[i][j]=strtod(massCenterCoordinate,NULL);
                }

                boxMatrix[0][0]=arg6; boxMatrix[1][1]=arg7; boxMatrix[1][0]=arg8; boxMatrix[0][1]=arg8;
                deltaR=arg9;
                arg=arg4; pacFrac=arg;
                rho=arg3; volume=N/rho;

                normalizingDenominator=pow(boxMatrix[1][0],2)-boxMatrix[0][0]*boxMatrix[1][1];
                actIndex=0;
                while (configurations[actIndex]!='{') actIndex++;
                actIndex+=3;
                for (int i=0;i<activeN;i++) {
                    for (int j=0;j<3;j++) {
                        char coordinate[50];
                        int licznik=0;
                        while (configurations[actIndex]!=',') coordinate[licznik++]=configurations[actIndex++];
                        if (j<2) {
                            actIndex++;
                            particles[i].r[j]=strtod(coordinate,NULL);
                        } else {
                            //actIndex+=22;         //dla mniejszej dokladnosci (%.5f) - dla kompatybilnosci ze starymi plikami
                            actIndex+=36;
                            particles[i].phi=strtod(coordinate,NULL);
                        }
                    }
                    particles[i].normR[0]=(-boxMatrix[1][1]*particles[i].r[0]+boxMatrix[1][0]*particles[i].r[1])/normalizingDenominator;
                    particles[i].normR[1]=(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/normalizingDenominator;
                }
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, startDen: %.4f, startPacFrac: %.4f, lambda[0]: %.4f, lambda[1]: %.4f, avPhi: %.4f, RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4f, massCenter: [%.4f,%.4f,%.4f,%.4f]\n",N,gaps,growing,arg3,pacFrac,lambda[0],lambda[1],avPhi,arg1,arg2,arg5,arg9,massCenter[0][0],massCenter[0][1],massCenter[1][0],massCenter[1][1]);
                //for (int i=0;i<2;i++) for (int j=0;j<2;j++) printf("boxMatrix[%d][%d]=%.12f\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<N;i++) printf("%.12f,  %.12f,  %.12f\n",particles[i].r[0],particles[i].r[1],particles[i].phi);return 0;
                updateNeighbourList(particles,boxMatrix);
            }
        }

        if (skipFirstIteration) {
            printf("Skipping first iteration...\n");
            skipFirstIteration=0;
        } else {
            if (loadedConfiguration) {
                if (loadedSetStartGenerator) {
                    printf("Setting start position of p-random number generator to position from file...\n");
                    InitMT((unsigned int)arg1);
                    randomStartStep[0]=arg1;
                } else {
                    randomStartStep[0]=InitRandomMT();
                    printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                }
                randomStartStep[1]=0;
                if (loadedSetGenerator) {
                    printf("Setting p-random number generator to last position from file...\n");
                    for (double i=0;i<arg2;i++) MTGenerate(randomStartStep);
                } else printf("Setting p-random number generator to last position from file - DISABLED\n");
            } else {
                printf("Setting start position of p-random number generator to actual CPU time...\n");
                randomStartStep[0]=InitRandomMT();
                randomStartStep[1]=0;
            }
            printf("Start of equilibration at pacFrac: %.4f (startDen: %.4f)... (%ld cycles)\n",pacFrac,rho,cyclesOfEquilibration);




/////////////////////////////////////////////// RDZEN MC

            long fullCycle=activeN,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0;
            double nStep=fullCycle*((double)cyclesOfEquilibration+(double)cyclesOfMeasurement+10.0),
                   possibleDistance=0;
            int cycleCounter=0, indexScanned=N==56?25:(N==224?105:(N==504?242:(N==780?377:(N==896?434:1534))));

            char allResultsFileName[200],bufferConfigurations[200],bufferOrientations[200],allOrientationsFileName[200],bufferOrientationsResults[200],allOrientationsResultsFileName[200],bufferPacFrac[100],bufferLambdaSetIndex[100];
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferOrientations,orientationsFileName); strcpy(allOrientationsFileName,orientationsFileName);
            strcpy(bufferOrientationsResults,orientationsResultsFileName); strcpy(allOrientationsResultsFileName,orientationsResultsFileName);
            sprintf(bufferPacFrac,"%.3f",pacFrac); sprintf(bufferLambdaSetIndex,"%d",lambdaSetIndex);
            strncat(allResultsFileName,"_arg-",6); strncat(allResultsFileName,bufferPacFrac,100); strncat(allResultsFileName,"_l-",4); strncat(allResultsFileName,bufferLambdaSetIndex,100); strncat(allResultsFileName,"_Results.txt",13);
            strncat(bufferConfigurations,"_arg-",6); strncat(bufferConfigurations,bufferPacFrac,100); strncat(bufferConfigurations,"_l-",4); strncat(bufferConfigurations,bufferLambdaSetIndex,100); strncat(bufferConfigurations,".txt",5);
            strncat(bufferOrientations,"_arg-",6); strncat(bufferOrientations,bufferPacFrac,100); strncat(bufferOrientations,"_l-",4); strncat(bufferOrientations,bufferLambdaSetIndex,100); strncat(bufferOrientations,".txt",5);
            strncat(allOrientationsFileName,"_arg-",6); strncat(allOrientationsFileName,bufferPacFrac,100); strncat(allOrientationsFileName,"_l-",4); strncat(allOrientationsFileName,bufferLambdaSetIndex,100); strncat(allOrientationsFileName,"_allOnt.txt",12);
            strncat(bufferOrientationsResults,"_arg-",6); strncat(bufferOrientationsResults,bufferPacFrac,100); strncat(bufferOrientationsResults,"_l-",4); strncat(bufferOrientationsResults,bufferLambdaSetIndex,100); strncat(bufferOrientationsResults,".txt",5);
            strncat(allOrientationsResultsFileName,"_arg-",6); strncat(allOrientationsResultsFileName,bufferPacFrac,100); strncat(allOrientationsResultsFileName,"_l-",4); strncat(allOrientationsResultsFileName,bufferLambdaSetIndex,100); strncat(allOrientationsResultsFileName,"_allOnt.txt",12);

            fileOrientations=fopen(bufferOrientations,"rt");
            if (fileOrientations==NULL) {
                fileOrientations=fopen(bufferOrientations,"a"); fprintf(fileOrientations,"{"); fclose(fileOrientations);
            } else if (!onlyMath[0]) {
                char bufferForEraseLastChar[200],linia[110]; strcpy(bufferForEraseLastChar,bufferOrientations); strncat(bufferForEraseLastChar,"_BUFF",6);
                FILE *bFELC = fopen(bufferForEraseLastChar,"w");
                int poziomNawiasu=0;
                while (fgets(linia,100,fileOrientations)!=NULL) {
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
                fclose(bFELC); fclose(fileOrientations);
                remove(bufferOrientations); rename(bufferForEraseLastChar,bufferOrientations);
            }
            fileAllResults = fopen(allResultsFileName,"a");
            fileOrientations = fopen(bufferOrientations,"a");
            fileAllOrientations = fopen(allOrientationsFileName,"a");
            if (onlyMath[0]) nStep=0;
            else {
                initPositions(latticeNodes,boxMatrix,true);
                computeDisplacementsSum(particles,latticeNodes,boxMatrix);
            }

            timeStart=time(0);
            for (double iStep=1;iStep<nStep;iStep++) {
                int randIndex = (int)(MTGenerate(randomStartStep)%1000000/1000000.0*fullCycle);
                attemptedNumberR++;
                if (attemptToDisplaceAParticle(particles,latticeNodes,randIndex,boxMatrix)) displacedNumberR++;

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++;

                    if (cycle%intervalSampling==0) {
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR;
                        possibleDistance+=sqrt(0.5*deltaR*deltaR)*((double)intervalSampling)*acceptanceRatioR*5.0;  //ostatnie *5.0 - dla bezpieczenstwa; sqrt(0.5*deltaR*0.5*deltaR+0.5*deltaR*0.5*deltaR)=sqrt(0.5*deltaR*deltaR)
                        if (possibleDistance>=neighSafeDistance) {
                            updateNeighbourList(particles,boxMatrix);
                            possibleDistance=0;
                        }

                        /////wypisywanie danych czesciej niz normalnie i PRZED zrownowagowaniem
                        /*if (cycle%50==0) {
                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
                                           normalizedRY=particles[i].normR[1]-particles[j].normR[1],
                                           rx=particles[i].r[0]-particles[j].r[0],
                                           ry=particles[i].r[1]-particles[j].r[1];
                                    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                                    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                                    double r2=rx*rx+ry*ry, dr=sqrt(r2);
                                    int energy;
                                    if (dr<maxDistance) {
                                        if (dr<minDistance) energy=1;
                                        else {
                                            double gamma=atan(ry/rx),
                                                    aAngle=particles[j].phi-gamma,
                                                    bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (rx>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);
                                            energy=checkOverlaps(dr,aAngle,bAngle);
                                            //energy=checkOverlapsAnalyticalMethodForEvenHCM(dr,aAngle,bAngle);
                                            if (energy==1) printf("colliding: distance- %.12f, alpha: %.12f, beta: %.12f, analytic: %.12f, i: %d, j: %d\n",dr,aAngle,bAngle,getMinimalDistanceAnalyticalMethodForEvenHCM(normalizeAngle(aAngle+C),normalizeAngle(bAngle+C)),i,j);
                                        }
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+arg5),collidingPairs);
                            else printf("Cycle: %ld\n",(cycle+arg5));
                            printf("   AccRatR: %.4f, dR: %.4f\n",acceptanceRatioR,deltaR);
                            printf("   Dens: %.4f, V/V_cp: %.4f\n",rho,pacFrac);
                            printf("   box00: %.8f, box11: %.8f, box10(01): %.8f\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                            printf("   displacementsSum[0]: %.8f, displacementsSum[1]: %.8f\n",displacementsSum[0],displacementsSum[1]);
                        }*/
                        /////

                        if (cycle>cyclesOfEquilibration) {
                            if (timeEquilibration==0) timeEquilibration=time(0);

                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    double normalizedRX=particles[i].normR[0]-particles[j].normR[0],
                                           normalizedRY=particles[i].normR[1]-particles[j].normR[1],
                                           rx=particles[i].r[0]-particles[j].r[0],
                                           ry=particles[i].r[1]-particles[j].r[1];
                                    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
                                    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
                                    double r2=rx*rx+ry*ry, dr=sqrt(r2);
                                    int energy;
                                    if (dr<maxDistance) {
                                        if (dr<minDistance) energy=1;  //analyticCheck 1/3
                                        else {
                                            double gamma=atan(ry/rx),
                                                    aAngle=particles[j].phi-gamma,
                                                    bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (rx>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle); bAngle=normalizeAngle(bAngle);  //analyticCheck 2/3
                                            energy=checkOverlaps(dr,aAngle,bAngle);
                                            //energy=checkOverlapsAnalyticalMethodForEvenHCM(dr,aAngle,bAngle);
                                        }  //analyticCheck 3/3
                                        if (energy==1) collidingPairs++;
                                    }
                                }
                            }

                            if (cycle%intervalResults==0)
                                fprintf(fileAllResults,"%ld\t%.17f\t%.17f\t\n",(cycle+arg5),displacementsSum[0]/(double)activeN,displacementsSum[1]/(double)activeN);

                            if (cycle%intervalOrientations==0) {
                                /////skan po konkretnej cząstce (np. dla 224: 14*(16/2)-(14/2)=105, etc.) - w srodku by nie skakala na granicy pudla periodycznego; bezposrednie uzycie w Mathematice (format tablicy)
                                if (cycle-cyclesOfEquilibration>=cyclesOfMeasurement)
                                    fprintf(fileOrientations,"{%.12f,%.12f,%.12f}}",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                else fprintf(fileOrientations,"{%.12f,%.12f,%.12f},",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                /////skan po wszystkich cząstkach
                                for (int i=0;i<activeN-1;i++) fprintf(fileAllOrientations,"%.12f,",particles[i].phi);
                                fprintf(fileAllOrientations,"%.12f\n",particles[activeN-1].phi);
                            }

                            if (cycle%intervalOutput==0) {
                                if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+arg5),collidingPairs);
                                else printf("Cycle: %ld\n",(cycle+arg5));
                                printf("   AccRatR: %.4f, dR: %.4f\n",acceptanceRatioR,deltaR);
                                printf("   Dens: %.4f, V/V_cp: %.4f\n",rho,pacFrac);
                                printf("   box00: %.8f, box11: %.8f, box10(01): %.8f\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                                printf("   displacementsSum[0]: %.8f, displacementsSum[1]: %.8f\n",displacementsSum[0],displacementsSum[1]);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12f\tV/V_cp: %.12f\tPressureReduced: pressureReduced\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        randomStartStep[0],randomStartStep[1],(cycle+arg5),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17f %.17f %ld %.17f %.17f %.17f %.17f %.17f,%.17f,%.17f,%.17f {",randomStartStep[0],randomStartStep[1],rho,pacFrac,(cycle+arg5),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,massCenter[0][0],massCenter[0][1],massCenter[1][0],massCenter[1][1]);
                                for (int i=0;i<activeN;i++)
                                    fprintf(fileConfigurations,"m[%.17f,%.17f,%.17f,%.12f,%.12f,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                                fprintf(fileConfigurations,"{Opacity->0.2,Red,Polygon[{{0,0},{%.12f,%.12f},{%.12f,%.12f},{%.12f,%.12f}}]},{Opacity->0.2,Green,Disk[{%.12f,%.12f},%.12f]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12f, boxMatrix[1][1]=%.12f, boxMatrix[1][0]=boxMatrix[0][1]=%.12f",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                                fclose(fileConfigurations);
                            }
                        } else {
                            if (acceptanceRatioR>desiredAcceptanceRatioR) deltaR*=1.05; else deltaR*=0.95;
                            if (deltaR>maxDeltaR) deltaR=maxDeltaR;
                        }
                        attemptedNumberR=0; displacedNumberR=0;
                    }
                    if (cycle-cyclesOfEquilibration>=cyclesOfMeasurement) iStep=nStep;
                }
            }
            fclose(fileAllResults);
            fclose(fileOrientations); fclose(fileAllOrientations);
            if (timeEquilibration==0) timeEquilibration=time(0);
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            printf("Start of calculation of results...\n");

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... ");
            double avDisplacements[2]={0,0};
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[activeN*17];
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);

            long dataLicznik=0;
            while(fgets(linia,300,fileAllResults)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++; if (actIndex>=300) continue;
                int dataIndex=0; double dataD[2]; while (dataIndex<2) {
                    char data[50];
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    dataD[dataIndex++]=strtod(data,NULL);
                }
                if (dataIndex<10) {
                    for (int i=0;i<2;i++) avDisplacements[i]+=dataD[i];
                    dataLicznik++;
                }
            }
            fclose(fileAllResults);
            for (int i=0;i<2;i++) avDisplacements[i]/=(double)dataLicznik;
            printf("done\n");

            //obliczenie bledow mierzonych wartosci
            printf("Calculation of errors of averages... ");
            double dAvDisplacements[2]={0,0};
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMath[1];i++) fgets(linia,300,fileAllResults);
            while(fgets(linia,300,fileAllResults)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                while (linia[actIndex]!='\t' && actIndex<300) actIndex++; actIndex++; if (actIndex>=300) continue;
                int dataIndex=0; double dataD[2]; while (dataIndex<2) {
                    char data[50];
                    int licznik=0;
                    while (linia[actIndex]!='\t' && licznik<50) data[licznik++]=linia[actIndex++]; if (licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    dataD[dataIndex++]=strtod(data,NULL);
                }
                if (dataIndex<10) for (int i=0;i<2;i++) {
                    double epsilon=avDisplacements[i]-dataD[i];
                    dAvDisplacements[i]+=epsilon*epsilon;
                }

            }
            fclose(fileAllResults);
            double denominator = (double)dataLicznik*((double)dataLicznik-1.0);
            for (int i=0;i<2;i++) dAvDisplacements[i]=getAvErrorFromSumEps(dAvDisplacements[i],denominator);
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla 1 czastki
            printf("Creation of 1-particle orientation file for Origin... ");
            double componentCounter=0, averageCos6PhiOne=0, ODFMaxOne=0;
            fileOrientations=fopen(bufferOrientations,"rt");
            double ODF_1P[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_1P[i]=0; //Orientational Distribution Function (1 Particle)
            int licznik=1;
            while (fgets(linia,2,fileOrientations)!=NULL) {
                sscanf(linia,"%c",linia);
                if (linia[0]==',') licznik++;
                if (licznik==3) {
                    licznik=0;
                    char data[50]; int actIndex=0;
                    while (true) {
                        fgets(linia,2,fileOrientations); sscanf(linia,"%c",linia);
                        if (linia[0]!='}') data[actIndex++]=linia[0];
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
            printf("Creation of ALL-particle orientation file for Origin... ");
            componentCounter=0; double averageCos6PhiAll=0, ODFMaxAll=0;
            fileAllOrientations = fopen(allOrientationsFileName,"rt");
            double ODF_AllP[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]=0; //Orientational Distribution Function (All Particles)
            while(fgets(linia,activeN*50,fileAllOrientations)!=NULL) {
                sscanf(linia,"%c",linia);
                int actIndex=0, licznikN=0;
                while (actIndex<activeN*30+10 && licznikN++<activeN) {
                    char data[50]; int licznik=0;
                    while (linia[actIndex]!=',' && licznik<30) data[licznik++]=linia[actIndex++]; actIndex++;
                    double angle = normalizeAngle(strtod(data,NULL)+C);
                    averageCos6PhiAll+=cos(6.0*angle); componentCounter++;
                    int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                    ODF_AllP[index]++;
                }
            }
            fclose(fileAllOrientations);
            suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_AllP[i]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]/=suma*dPhi;
            averageCos6PhiAll/=componentCounter; for (int i=0;i<ODFLength;i++) if (ODFMaxAll<ODF_AllP[i]) ODFMaxAll=ODF_AllP[i];
            printf("done\n");

            long timeEndMath=time(0);




/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... ");
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) {
                fileConfigurations = fopen(bufferConfigurations,"w");
                fileConfigurationsList = fopen(configurationsListFileName,"a");
            }
            fileOrientationsResults = fopen(bufferOrientationsResults,"w");
            fileAllOrientationsResults = fopen(allOrientationsResultsFileName,"w");

            fprintf(fileResults,"%ld\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",(cycle+arg5),volume,boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],rho,pacFrac,lambda[0],avDisplacements[0],dAvDisplacements[0],lambda[1],avDisplacements[1],dAvDisplacements[1],ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll);
            fprintf(fileExcelResults,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",volume,rho,pacFrac,lambda[0],avDisplacements[0],lambda[1],avDisplacements[1],ODFMaxAll,averageCos6PhiAll);

            if (!onlyMath[0]) {
                fprintf(fileConfigurations,"Rho: %.12f\tV/V_cp: %.12f\tPressureRed: pressureReduced\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        randomStartStep[0],randomStartStep[1],(cycle+arg5),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17f %.17f %ld %.17f %.17f %.17f %.17f %.17f,%.17f,%.17f,%.17f {",randomStartStep[0],randomStartStep[1],rho,pacFrac,(cycle+arg5),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,massCenter[0][0],massCenter[0][1],massCenter[1][0],massCenter[1][1]);
                fprintf(fileConfigurationsList,"multimers[x_,y_,kI_]:={");
                for (int i=0;i<activeN;i++) {
                    fprintf(fileConfigurations,"m[%.17f,%.17f,%.17f,%.12f,%.12f,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                    fprintf(fileConfigurationsList,"m[%.12f+x,%.12f+y,%.12f,%.12f,%.12f,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                }
                fprintf(fileConfigurations,"{Opacity->0.2,Red,Polygon[{{0,0},{%.12f,%.12f},{%.12f,%.12f},{%.12f,%.12f}}]},{Opacity->0.2,Green,Disk[{%.12f,%.12f},%.12f]},{Opacity->0.7,Black,Disk[{%.12f,%.12f},%.12f]},{Opacity->0.4,Black,Disk[{%.12f,%.12f},%.12f]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,massCenterLattice[0][0],massCenterLattice[0][1],multimerD,massCenter[0][0],massCenter[0][1],multimerD);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12f, boxMatrix[1][1]=%.12f, boxMatrix[1][0]=boxMatrix[0][1]=%.12f",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[1][0]);
                fprintf(fileConfigurationsList,"{Opacity->If[x==0 && y==0,0.4,0],Red,Polygon[{{0,0},{%.12f,%.12f},{%.12f,%.12f},{%.12f,%.12f}}]},{Opacity->If[x==0 && y==0,0.4,0],Green,Disk[{%.12f,%.12f},%.12f]}};\nconfigurationsList=Append[configurationsList,g[%.12f,multimers[%.12f,%.12f,kolorIndex=1],multimers[%.12f,%.12f,kolorIndex=1],multimers[%.12f,%.12f,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                        boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,pacFrac,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
            }

            for (int i=0;i<ODFLength;i++) {
                fprintf(fileOrientationsResults,"%.12f\t%.12f\n",-C+i*dPhi,ODF_1P[i]);
                fprintf(fileAllOrientationsResults,"%.12f\t%.12f\n",-C+i*dPhi,ODF_AllP[i]);
            }

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) {
                fclose(fileConfigurations); fclose(fileConfigurationsList);
            }
            fclose(fileOrientationsResults); fclose(fileAllOrientationsResults);
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        arg=getNextArgument(arg,true);
        oldVolume=volume;
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        arg5=0;
        loadedConfiguration=0;
    }
    printf("\nTime for equilibrations: %ldsec ,  time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n");
}
