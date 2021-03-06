Lambda <- 4;

Bonferroni <- FALSE;

#QueryStartID <-  "SPy0002";  
#QueryEndID <-  "SPy0040"

QueryStartID <-  "SPy0002";  
QueryEndID <-  "SPy2217"

Genome <- FALSE;

MintStat <- 0;


            ################# INPUT / OUTPUT  #################

StatDataDirectory <- "/home/kirkb/Playarea/Bayesian Clustering/Distribution Depend/Data/Input/CyberT/";

StatDataFile1 <- "G8bA8P2tStatNOL4.GeCyTBH";
StatDataFile2 <- "G8bA8P2tStatNOL4.CyTBH";
StatHeaderLength <- 9;  #DO NOT EDIT


OutputDirectory <- "/home/kirkb/Playarea/Bayesian Clustering/Distribution Depend/Data/Output/";

OutputFile <- "G8bA8P2tStatNO";

# From PgKMicroarrayWindow.R
###############################################################################
#                                FUNCTIONS                                    #
###############################################################################


###############################################################################
#                           CONSTANT DECLARATION                              #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x, y, or z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANENCY


                 ############################################
                 #######  DEFINE PERMANENT CONSTANTS  #######
                 ############################################

#Define Input Directory and File
StatDataIn1 <- paste(StatDataDirectory, StatDataFile1, sep="");
StatDataIn2 <- paste(StatDataDirectory, StatDataFile2, sep="");

#Define Output Dirctory with Prefix
OutputPrefix <- paste(OutputDirectory, OutputFile, sep="");


             #####################################################
             ##########  DEFINE TEMPORARY CONSTANTS  #############
             #####################################################

#Data Positions in CyberT .GeCyTBH output file
xStatIDPos1 <- 1;
xStattStatPos1 <- 2;
xStatpValuePos1  <- 3;
xStatBHpValuePos1  <- 4;
xStatBonpValuePos1 <- 5;
xStatGenomePositionPos1 <- 6;

#Data Positions in CyberT .CyTBH output file
xStatIDPos2 <- 2;
xStatMPos2 <- 15;
xStatAdjStndErrPos2 <- 16;
xStatArrayNumberPos2 <- 17;
xStatIntensityPos2 <- 18;

###############################################################################
#                             VARIABLE DECLARATION                            #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x", "y", or "z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANCY


                        #####  Input Tables  ######

#Input StatSigData 
xStatDataTable1  <- data.frame(read.delim(StatDataIn1,
                                       skip=StatHeaderLength));

xStatDataTable2  <- data.frame(read.delim(StatDataIn2,
                                       skip=StatHeaderLength));



                 #####  Create Vectors from Input Tables  #####

#Data from StatSigTable
ID <- as.character(as.vector(xStatDataTable1[,xStatIDPos1]));
tStat <- as.numeric(as.vector(xStatDataTable1[, xStattStatPos1]));
pValue  <- as.numeric(as.vector(xStatDataTable1[, xStatpValuePos1]));
BHpValue  <- as.numeric(as.vector(xStatDataTable1[, xStatBHpValuePos1]));
BonpValue <- as.numeric(as.vector(xStatDataTable1[, xStatBonpValuePos1]));
GenomePos <- as.numeric(as.vector(
                                 xStatDataTable1[, xStatGenomePositionPos1]));


#Find Stat Lambda
xLambdaScan <- scan(file=StatDataIn1, what=character(0), sep="\t",
                        nlines=StatHeaderLength);
xLambdaString <- grep("Lambda = ", xLambdaScan,
                          value=TRUE);

if (length(xLambdaString) != 0) {
  StatLambdaFound <- TRUE;
  xStringSplit <- unlist(strsplit(xLambdaString, "= "));
  StatLambda <- as.character(xStringSplit[2]); 
} else {
  StatLambdaFound <- FALSE;
  StatLambda <- "NotFound";
}

## Oragnize select data from second Stat file relative to the first
StatID2 <- as.character(as.vector(xStatDataTable2[,xStatIDPos2]));
StatM <- as.numeric(as.vector(xStatDataTable2[, xStatMPos2]));
StatAdjStndErr <- as.numeric(as.vector(xStatDataTable2[, xStatAdjStndErrPos2]));
StatArrayNumber <- as.numeric(as.vector(xStatDataTable2[, xStatArrayNumberPos2]));
StatIntensity <- as.numeric(as.vector(xStatDataTable2[, xStatIntensityPos2]));

LogFoldChange <- rep(0, length(ID));
Intensity <- rep(0, length(ID));
AdjStndErr <- rep(0, length(ID));
ArrayNumber <- rep(0, length(ID));

for (i in 1:length(StatID2)) {
  GeneFilter <- as.logical(match(ID, StatID2[i], nomatch=0));

  LogFoldChange[GeneFilter] <- StatM[i];
  Intensity[GeneFilter] <- StatIntensity[i];
  AdjStndErr[GeneFilter] <- StatAdjStndErr[i];
  ArrayNumber[GeneFilter] <- StatArrayNumber[i];
}

#Index GenomeGenes and Define Query Start and Query End

QueryStartPos <- GenomePos[as.logical(match(ID, QueryStartID,
                                             nomatch=0))];
QueryEndPos <- GenomePos[as.logical(match(ID, QueryEndID,
                                             nomatch=0))];

rm(list=ls(pat="^x"));

###############################################################################
#                               MAIN BODY                                     #
###############################################################################
#Search for Statistically Significant Clusters Using a Dynamic Windowing
#Algorithm

ClStartID <- "blank";
ClEndID <- "blank";
ClStartPos <- 0;
ClEndPos <- 0;

ClLength <- 0;

DegreesOfFreedom <- as.numeric(StatLambda)-1;

if (!Genome) {
for (i in QueryStartPos:QueryEndPos) {
  if (tStat[i] > MintStat) {
    Search <- 1;
    
    while ((tStat[(i+Search)] > MintStat) && ((Search+i) <= QueryEndPos)) {
            
      ClStartPos <- c(ClStartPos, i);
      ClEndPos <- c(ClEndPos, i + Search);

      ClStartID <- c(ClStartID, ID[i]);
      ClEndID <- c(ClEndID, ID[i+Search]);

      ClLength <- c(ClLength, Search + 1);

      Search <- Search + 1;
    }
  }
}

ClStartID <- ClStartID[2:length(ClStartID)];
ClEndID <- ClEndID[2:length(ClEndID)];
ClStartPos <- ClStartPos[2:length(ClStartPos)];
ClEndPos <- ClEndPos[2:length(ClEndPos)];
ClLength <- ClLength[2:length(ClLength)];
} else {
ClStartID <- ID[QueryStartPos];
ClEndID <- ID[QueryEndPos];;
ClStartPos <- QueryStartPos;
ClEndPos <- QueryEndPos;
ClLength <- QueryEndPos-QueryStartPos+1;

}

GeName <- "blank";
GeClLength <- 0;

ReftStat <- 0;

Ku <- 0;
Ksub <- 0;

Pg <- 0;
PKsub <- 0;
PKu <- 0;
PKmi <- 0;
PKet <- 0;
 
PgKsub <- 0;
PgKu <- 0;
PgKmi <- 0;
PgKet <- 0;

PgEKu <- 0;
PgIKu <- 0;

PKmiKu <- 0;

rPg <- 0;
rPKsub <- 0;
rPKu <- 0;
rPKmi <- 0;
rPKet <- 0;

GetStat <- 0
GeLogFoldChng <- 0;
GeStndErr <- 0;
GeIntensity <- 0;
GeArrayNumber <- 0;

ClBHpValue <- 0;

for (i in 1:length(ClStartID)) {
  
  for (j in ClStartPos[i]:ClEndPos[i]) {
    if (tStat[j] > 0) {

      ClPos <- 1:ClLength[i];
    
      ReftStat <- c(ReftStat, tStat[j]);

      tKsubSet <- tStat[ClStartPos[i]:ClEndPos[i]];

      tKsubSet <- tKsubSet[which(ClPos != (j-ClStartPos[i] + 1))];
    
      Ku <- c(Ku, sum(tStat[ClStartPos[i]:ClEndPos[i]]));

      Ksub <- c(Ksub, sum(tKsubSet));

  
      Pg <- c(Pg, 1-2*(1-pt(ReftStat[length(ReftStat)], DegreesOfFreedom)));

      PKsub <- c(PKsub, 1-2*(1-pt(Ksub[length(Ksub)], DegreesOfFreedom)));

      PKu <- c(PKu, 1-2*(1-pt(Ku[length(Ku)], DegreesOfFreedom)));

      tPKet <- PKu[length(PKu)] - PKsub[length(PKsub)] ;
      for (tg in tKsubSet) {
        tKsub <- Ku[length(Ku)] - tg;
        tPgE <- PKu[length(PKu)] - (1-2*(1-pt(tKsub, DegreesOfFreedom)));
        tPKet <- tPKet + tPgE;

      }


      PKet <- c(PKet, tPKet);

      PKmi <- c(PKmi, PKu[length(PKu)] - PKet[length(PKet)]);
  
      PgKsub <-c(PgKsub,(Pg[length(Pg)]+PKsub[length(PKsub)]-PKu[length(PKu)])/
                  PKsub[length(PKsub)]);
  
      PgKu <- c(PgKu, Pg[length(Pg)] / PKu[length(PKu)]);

      PgKmi <- c(PgKmi,(Pg[length(Pg)]+PKsub[length(PKsub)]-PKu[length(PKu)]) /
                 PKmi[length(PKmi)]);

      PgKet <- c(PgKet, (PKu[length(PKu)] - PKsub[length(PKsub)]) /
                 PKet[length(PKet)]);

      PgEKu <- c(PgEKu, (PKu[length(PKu)] - PKsub[length(PKsub)]) /
                 PKu[length(PKu)]);

      PgIKu <- c(PgIKu,(Pg[length(Pg)]+PKsub[length(PKsub)]-PKu[length(PKu)]) /
                 PKu[length(PKu)]);
      

      GeName <- c(GeName, ID[j]);
      GeClLength <- c(GeClLength, ClLength[i]);
      GetStat <- c(GetStat, tStat[j]);
      GeLogFoldChng <- c(GeLogFoldChng, LogFoldChange[j]);
      GeStndErr <- c(GeStndErr, AdjStndErr[j]);
      GeIntensity <- c(GeIntensity, Intensity[j]);
      GeArrayNumber <- c(GeArrayNumber, ArrayNumber[j]);
            
    } else {  ## tStat = 0

    
      ReftStat <- c(ReftStat, tStat[j]);

      tKsubSet <- tStat[ClStartPos[i]:ClEndPos[i]];

      tKsubSet <- tKsubSet[which(ClPos != (j-ClStartPos[i] + 1))];
    
      Ku <- c(Ku, sum(tStat[ClStartPos[i]:ClEndPos[i]]));

      Ksub <- c(Ksub, sum(tKsubSet));


      Pg <- c(Pg, 0);

      PKsub <- c(PKsub, 1-2*(1-pt(Ksub[length(Ksub)], DegreesOfFreedom)));

      PKu <- c(PKu, 1-2*(1-pt(Ku[length(Ku)], DegreesOfFreedom)));


      tPKet <- PKu[length(PKu)] - PKsub[length(PKsub)] ;
      for (tg in tKsubSet) {
        tKsub <- Ku[length(Ku)] - tg;
        tPgE <- PKu[length(PKu)] - (1-2*(1-pt(tKsub, DegreesOfFreedom)));
        tPKet <- tPKet + tPgE;

      }


      PKet <- c(PKet, tPKet);

      PKmi <- c(PKmi, PKu[length(PKu)] - PKet[length(PKet)]);
      
      PgKsub <- c(PgKsub, 0);
  
      PgKu <- c(PgKu, 0);

      PgKmi <- c(PgKmi,(Pg[length(Pg)]+PKsub[length(PKsub)]-PKu[length(PKu)]) /
                 PKmi[length(PKmi)]);

      PgKet <- c(PgKet, (PKu[length(PKu)] - PKsub[length(PKsub)]) /
                 PKet[length(PKet)]);

      PgEKu <- c(PgEKu, 0);

      PgIKu <- c(PgIKu, 0);

      PKmiKu <- c(PKmiKu, PKmi/PKu); 

      GeName <- c(GeName, ID[j]);
      GeClLength <- c(GeClLength, ClLength[i]);
      GetStat <- c(GetStat, tStat[j]);
      GeLogFoldChng <- c(GeLogFoldChng, LogFoldChange[j]);
      GeStndErr <- c(GeStndErr, AdjStndErr[j]);
      GeIntensity <- c(GeIntensity, Intensity[j]);
      GeArrayNumber <- c(GeArrayNumber, ArrayNumber[j]);

    }
    
   
    
  }

  PKmiKu <- c(PKmiKu, PKmi[length(PKmi)]/PKu[length(PKmi)]);

  tClPos <- ClStartPos[i]:ClEndPos[i];

  tpValue <- pValue[tClPos];
  tpValueFilter <- order(tpValue, decreasing=FALSE);
  tpValue <- tpValue[tpValueFilter];

  RankedClPos <- tClPos[tpValueFilter];
  tBHpValue <- rep(1, ClLength[i]);
  for (tRank in 1:ClLength[i]) {
    tBHpValue[tRank] <- (tpValue[tRank] * ClLength[i]) / tRank;
  }
  tBHpValue[which(tBHpValue>1)] <- 1;
  for (tRank in 1:ClLength[i]) {
    tBHpValue[tRank] <- min(tBHpValue[tRank:ClLength[i]]);   
  }

  ReorderFilter <- order(RankedClPos, decreasing=FALSE);
  ClBHpValue <- c(ClBHpValue, tBHpValue[ReorderFilter]);
  
  
}


GeName <- GeName[2:length(GeName)];
GetStat <- GetStat[2:length(GeName)];
GeLogFoldChng <- GeLogFoldChng[2:length(GeLogFoldChng)];
GeStndErr <- GeStndErr[2:length(GeStndErr)];
GeIntensity <- GeIntensity[2:length(GeIntensity)];
GeArrayNumber <- GeArrayNumber[2:length(GeArrayNumber)];
GeClLength <- GeClLength[2:length(GeClLength)];


ReftStat <- ReftStat[2:length(ReftStat)];
Ku <- Ku[2:length(Ku)];
Ksub <- Ksub[2:length(Ksub)];


Pg <- Pg[2:length(Pg)];
PKsub <- PKsub[2:length(PKsub)];
PKu <- PKu[2:length(PKu)];
PKmi <- PKmi[2:length(PKmi)];
PKet <- PKet[2:length(PKet)];

PgKu <- PgKu[2:length(PgKu)];
PgIKu <- PgIKu[2:length(PgIKu)];
PgEKu <- PgEKu[2:length(PgEKu)];

PgKsub <- PgKsub[2:length(PgKsub)];

PgKmi <- PgKmi[2:length(PgKmi)];
PgKet <- PgKet[2:length(PgKet)];

PKmiKu <- PKmiKu[2:length(PKmiKu)];

ClBHpValue <- ClBHpValue[2:length(ClBHpValue)];

tPKmiKu <- PKmiKu[which(PKmiKu > 0.95)];
tClStartPos <- ClStartPos[which(PKmiKu > 0.95)];
tClEndPos <- ClEndPos[which(PKmiKu > 0.95)];
tClStartID <- ClStartID[which(PKmiKu > 0.95)];
tClEndID <- ClEndID[which(PKmiKu > 0.95)];
tClLength <- ClLength[which(PKmiKu > 0.95)];

UberClStartPos <- 0;
UberClEndPos <- 0;
UberClStartID <- "None";
UberClEndID <- "None";
UberClLength <- 0;

for (k in (1:length(tPKmiKu))) {
  if (tClStartPos[k] > UberClEndPos[length(UberClEndPos)]) {
    
    UberClStartPos <- c(UberClStartPos, tClStartPos[k]);
    UberClEndPos <- c(UberClEndPos, tClEndPos[k]);

    UberClStartID <- c(UberClStartID, tClStartID[k]);
    UberClEndID <- c(UberClEndID, tClEndID[k]);;

    
  } else {
       
    if (tClEndPos[k] > UberClEndPos[length(UberClEndPos)]) {
      UberClEndPos[length(UberClEndPos)] <- tClEndPos[k];
      UberClEndID[length(UberClEndID)] <- tClEndID[k];
    }
    
  }
  
  
UberClLength <- UberClEndPos - UberClStartPos + 1;

}

UberClStartPos <- UberClStartPos[2:length(UberClStartPos)];
UberClEndPos <- UberClEndPos[2:length(UberClEndPos)];

UberClStartID <- UberClStartID[2:length(UberClStartID)];
UberClEndID <- UberClEndID[2:length(UberClEndID)];
         
UberClLength <- UberClLength[2:length(UberClLength)];


pg <- 1-Pg;

pgKu <- 1-PgKu;
pgIKu <- 1-PgIKu;
pgEKu <- 1-PgEKu;

pgKsub <- 1-PgKsub;

pgKmi <- 1-PgKmi;
pgKet <- 1-PgKet;


Bon <- GeClLength * (1-Pg);
Sidak <- 1-Pg**GeClLength;

p <- PKu - Pg;



                       ############################
                       #####  OUTPUT RESULTS  #####
                       ############################

postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/PhenomBHDiffMicroArrayDoF", as.character(DegreesOfFreedom), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);

library(grDevices);

rgb.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"),
                                space = "rgb");

Diff <- p-ClBHpValue;

Prob.palette <- rgb.palette(length(Diff));
DiffFilter <- order(Diff, decreasing=FALSE);

tDiff <- Diff[DiffFilter];
 
tPK <- PKsub[DiffFilter];
tPg <- Pg[DiffFilter];
tPK2 <- GeClLength[DiffFilter];

Min <- 0;
Max <- 1;

#x11(width=11, height=8);


library(lattice);

x11(width=11, height=8);
cloud(tPK2~tPg*tPK,
      screen=list(z=-135, x=-90),
      pch=16,
      col=Prob.palette);


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/PhenomBHFracMicroArray1DoF", as.character(DegreesOfFreedom), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);

library(grDevices);

rgb.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"),
                                space = "rgb");

DiffPercent <- (pgKsub-ClBHpValue)/pgKsub;

Prob.palette <- rgb.palette(length(DiffPercent));
DiffPercentFilter <- order(DiffPercent, decreasing=FALSE);

tDiffPercent <- DiffPercent[DiffPercentFilter];
 
tPK <- PKsub[DiffPercentFilter];
tPg <- Pg[DiffPercentFilter];
tPK2 <- GeClLength[DiffPercentFilter];

Min <- 0;
Max <- 1;

#x11(width=11, height=8);


library(lattice);

x11(width=11, height=8);
cloud(tPK2~tPg*tPK,
      screen=list(z=-135, x=-90),
      pch=16,
      col=Prob.palette);

####
plot(tPK~tPg, type="p",
     xlab="", ylab="", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=2, pch=16, col=Prob.palette);


points(tPg[which(abs(tDiff) < 0.20)],
       tPK[which(abs(tDiff) < 0.20)],
       cex=0.5, pch=16, col="white");

points(tPg[which(abs(tDiff) < 0.10)],
       tPK[which(abs(tDiff) < 0.10)],
       cex=0.5, pch=16, col="grey60");

points(tPg[which(abs(tDiff) < 0.05)],
       tPK[which(abs(tDiff) < 0.05)],
       cex=0.5, pch=16, col="black");
dev.off()


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/PhenomBHDiffPerMicroArrayDoF", as.character(DegreesOfFreedom), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);

DiffPercent <- p/ClBHpValue;

Prob.palette <- rgb.palette(length(DiffPercent));
DiffPercentFilter <- order(DiffPercent, decreasing=FALSE);

tDiffPercent <- DiffPercent[DiffPercentFilter];
 
tPK <- PKsub[ProbFilter];
tPg <- Pg[ProbFilter];

Min <- 0;
Max <- 1;

#x11(width=11, height=8);

plot(tPK~tPg, type="p",
     xlab="", ylab="", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=2, pch=16, col=Prob.palette);


points(tPg[which(tDiffPercent < 0.20)],
       tPK[which(tDiffPercent < 0.20)],
       cex=0.5, pch=16, col="white");

points(tPg[which(tDiffPercent < 0.10)],
       tPK[which(tDiffPercent < 0.10)],
       cex=0.5, pch=16, col="grey60");

points(tPg[which(tDiffPercent < 0.05)],
       tPK[which(tDiffPercent < 0.05)],
       cex=0.5, pch=16, col="black");
dev.off()


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/PgKmi3DDynamicMicroArrayDoF", as.character(DegreesOfFreedom), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);


library(grDevices);

rgb.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"),
                                space = "rgb");

Prob.palette <- rgb.palette(length(PgKmi));
ProbFilter <- order(PgKmi, decreasing=FALSE);

tProb <- PgKmi[ProbFilter];
 
tPK <- PKsub[ProbFilter];
tPK2 <- PKmi[ProbFilter];
tPg <- Pg[ProbFilter];

library(lattice);

#x11(width=11, height=8);
cloud(tPK2~tPg*tPK,
      screen=list(z=-120, x=-45),
      pch=16,
      col=Prob.palette);

dev.off();
#tPK <- tPK[which(GeClLength > 2)];
#tPg <- tPg[which(GeClLength > 2)];


Min <- 0;
Max <- 1;

#x11(width=11, height=8);

plot(tPK~tPg, type="p",
     xlab="", ylab="", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=2, pch=16, col=Prob.palette);

points(tPg[which(tProb > 0.80)],
       tPK[which(tProb > 0.80)],
       cex=0.5, pch=16, col="white");

points(tPg[which(tProb > 0.90)],
       tPK[which(tProb > 0.90)],
       cex=0.5, pch=16, col="grey60");

points(tPg[which(tProb > 0.95)],
       tPK[which(tProb > 0.95)],
       cex=0.5, pch=16, col="black");
dev.off()


points(tPg[which(tProb < 0.20)],
       tPK[which(tProb < 0.20)],
       cex=0.5, pch=16, col="white");

points(tPg[which(tProb < 0.10)],
       tPK[which(tProb < 0.10)],
       cex=0.5, pch=16, col="grey60");

points(tPg[which(tProb < 0.05)],
       tPK[which(tProb < 0.05)],
       cex=0.5, pch=16, col="black");
dev.off()



postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/PgEKuScatterDoF", as.character(DegreesOfFreedom), "DotThree.eps", sep=""),  horizontal=TRUE, onefile=TRUE);


library(grDevices);

rgb.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"),
                                space = "rgb");

Prob.palette <- rgb.palette(length(Bon));
ProbFilter <- order((1-Bon), decreasing=FALSE);
tProb <- 1-Bon[ProbFilter];
 
tPK <- PKsub[ProbFilter];
tPg <- Pg[ProbFilter];

#tPK <- tPK[which(GeClLength > 2)];
#tPg <- tPg[which(GeClLength > 2)];

Min <- 0;
Max <- 1;

#x11(width=11, height=8);

plot(tPK~tPg, type="p",  main=paste("Bon dof =",
     as.character(DegreesOfFreedom), sep=" "),
     xlab="Pg", ylab = "PKsub", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=1, pch=16, col=Prob.palette);

points(tPg[which(tProb > 0.80)],
       tPK[which(tProb > 0.80)],
       cex=0.5, pch=16, col="white");

points(tPg[which(tProb > 0.90)],
       tPK[which(tProb > 0.90)],
       cex=0.5, pch=16, col="grey60");

points(tPg[which(tProb > 0.95)],
       tPK[which(tProb > 0.95)],
       cex=0.5, pch=16, col="black");


x11(width=11, height=8);


Min <- 0;
Max <- 1;

plot(PKsub~Pg, type="p",  main=paste("dof =",
     as.character(DegreesOfFreedom), sep=" "),
     xlab="Pg", ylab = "PKsub", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=0.5, pch=16, col="black");

x11(width=11, height=8);


Min <- 0;
Max <- 1;


plot(rPKsub~rPg, type="p",  main=paste("dof =",
     as.character(DegreesOfFreedom), sep=" "),
     xlab="Pg", ylab = "PKsub", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=0.5, pch=16, col="red");

x11(width=11, height=8);


Min <- 0;
Max <- 1;

plot(PKmi~Pg, type="p",  main=paste("dof =",
     as.character(DegreesOfFreedom), sep=" "),
     xlab="Pg", ylab = "PKmi", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=0.5, pch=16, col="black");

x11(width=11, height=8);


Min <- 0;
Max <- 1;


plot(rPKmi~rPg, type="p",  main=paste("dof =",
     as.character(DegreesOfFreedom), sep=" "),
     xlab="Pg", ylab = "PKmi", xlim=c(Min, Max), ylim=c(Min, Max),
     cex=0.5, pch=16, col="red");


##############################################

BonFilter <-  order(Bon, decreasing=FALSE);
SidakFilter <-  order(Sidak, decreasing=FALSE);
PgKuFilter <-  order(PgKu, decreasing=FALSE);
PgKsubFilter <-  order(PgKsub, decreasing=FALSE);
BHpValueFilter <-  order(BHpValue, decreasing=FALSE);


Bon[BonFilter];
BHpValue[BHpValueFilter];
Sidak[SidakFilter];
PgKu[PgKuFilter];
PgKsub[PgKsubFilter];

tempGeneData <- cbind(GeName, Bon, Sidak, PgKsub, PgKu,
                      PGeKIntersect, PGe, PKu,
                      PKsub, ClStartName, ClEndName, ClStartPos, ClEndPos);

######################  PLOTS  ####################
####??????? Debug b4 use

n <- 5;
Finish <- trunc(1000/n);
ntPKu <- 0;

for (i in 1:Finish) {
  RandomFilter <- sample(1:length(trtStat), n);

  tKutStat <- sum(trtStat[RandomFilter]);
  ntPKu <- c(ntPKu, 1 - (2 * (1-pt(tKutStat, DegreesOfFreedom))));

}

ntPKu <- ntPKu[2:length(ntPKu)];

PlotLimit <- 0.5;
x11(width=11, height=8);
plot(density(ntPKu), main="PKu  n = ",xlim=c(0, 1.2));










#############################################################################
PosteriorPos <- 2;
ClLengthPos <- 4;

Posterior <- as.numeric(as.vector(tempGeneData[,PosteriorPos])) * ClLength;

Fill <- c(1, 1, 1, 0, 0, 0, 0, NA, NA, 0, 0);


PlotData <- tempGeneData[1,]; #To be removed later only used to prime



for (i in ID) {
  GeneFilter <- as.logical(match(GeName, i, nomatch=0));

  if (sum(GeneFilter) > 0) {

    GeGeneData <- tempGeneData[GeneFilter,];

    if (sum(GeneFilter) == 1) {
      PlotData <- rbind(PlotData, GeGeneData);     
    } else {
      subPosterior <- Posterior[GeneFilter];
  
      ProbFilter <- order(subPosterior);
      OrderedGeGeneData <- GeGeneData[ProbFilter,];
      PlotData <- rbind(PlotData, OrderedGeGeneData[1,]);
    }
  } else {

    tempFill <- c(i, Fill)
    PlotData <- rbind(PlotData, tempFill);     

  }
}

#rm(tempGeneData);

PlotData <- PlotData[(2:length(PlotData[,1])),];


PostBonpValue <- as.numeric(as.vector(PlotData[,PosteriorPos])) *
                                 as.numeric(as.vector(PlotData[,ClLengthPos]));

PostBonpValue[which(PostBonpValue > 1)] <- 1;


BonPairpValue <- 2 * (1-PairPost);
BonPairpValue[which(BonPairpValue > 1)] <- 1;


OutPlotData <- cbind(PlotData[,1], tStat, LogFoldChange, Intensity,
                     AdjStndErr, ArrayNumber, BHpValue, BonpValue,
                     PostBonpValue, BonPairpValue,
                     PlotData[, 2:length(PlotData[1,])]);


# Plot Data
FilterHeader1 <- paste("StatData File1 = ", StatDataIn1, sep="");
FilterHeader2 <- paste("StatData File2 = ", StatDataIn2, sep="");
FilterHeader3 <- c(paste("Query Start = ", as.character(QueryStartID),
                                           sep=""),
                 paste("Query End = ", as.character(QueryEndID), sep=""));
FilterHeader4 <- paste("Window Limit = ", as.character(WindowLimit),
                                           sep="");   
FilterHeader5 <- c("Stat Param:",
                   paste("Lambda = ", as.character(StatLambda), sep=""));
FilterHeader6 <- c(paste("Bonferroni = ",
                         as.character(Bonferroni), sep=""));
FilterHeader7 <- "Data Generated by GCKu.R";
FilterHeader8 <- date();

FilterHeaderTable <-  c("Gene", "tStat", "LogFoldChng", "RMSInt",
                        "AdjStndErr", "ArrayNumber", "BHpValue", "BonpValue",
                        "Bon Post p", "Bon Pair Post p", "Post p",
                        "Ge p", "Cl Size", "Ge I Cl P",
                        "Ge P", "P Ku", "P Ksub", "Cluster Start",
                        "Cluster End", "Cl Start Pos", "Cl End Pos");

cat(FilterHeader1, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n");
cat(FilterHeader2, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeader3, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeader4, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeader5, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeader6, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeader7, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeader8, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);
cat(FilterHeaderTable, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeKuPD", sep=""), sep="\n", append=TRUE);

write.table(OutPlotData, file=paste(OutputPrefix, 
                               "L", StatLambda, ".GeKuPD", sep=""),
                               quote=FALSE, sep="\t", col.names=FALSE,
                               row.names=FALSE, append=TRUE);



rm(PlotData);

