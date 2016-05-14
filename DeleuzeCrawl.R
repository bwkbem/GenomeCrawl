
#From PgKMicroarrayWhole.R
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
StatDataGenome <- paste(StatDataDirectory, StatDataGenomeFile, sep="");

#Define Output Dirctory with Prefix
OutputPrefix <- paste(OutputDirectory, OutputFile, sep="");


             #####################################################
             ##########  DEFINE TEMPORARY CONSTANTS  #############
             #####################################################

#Data Positions in CyberT .GeCyTDS output file
xIDPos <- 1;
xGenomePositionPos <- 2;

if (CyberT) { 
    xtStatPos <- 3;
    xpValuePos  <- 5;
} else {
    xtStatPos <- 4;
    xpValuePos  <- 6;
}

xMPos <- 7;
xIntensityPos <- 8;
xBioSamplePos <- 9;
xArrayNumberPos <- 10;


###############################################################################
#                             VARIABLE DECLARATION                            #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x", "y", or "z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANCY


                        #####  Input Table  ######

#Input StatSigData 
xStatDataTable  <- data.frame(read.delim(StatDataGenome,
                                       skip=StatHeaderLength));


                 #####  Create Vectors from Input Table  #####

#Data from StatSigTable
ID <- as.character(as.vector(xStatDataTable[,xIDPos]));
GenomePos <- as.numeric(as.vector(
                                 xStatDataTable[, xGenomePositionPos ]));
tStat <- as.numeric(as.vector(xStatDataTable[, xtStatPos]));
pValue  <- as.numeric(as.vector(xStatDataTable[, xpValuePos]));
LogFoldChange <- as.numeric(as.vector(xStatDataTable[, xMPos]));
Intensity <- as.numeric(as.vector(xStatDataTable[, xIntensityPos]));
BioSampleNumber <- as.numeric(as.vector(xStatDataTable[, xBioSamplePos]));
ArrayNumber <- as.numeric(as.vector(xStatDataTable[, xArrayNumberPos]));

#Find Stat Lambda
xLambdaScan <- scan(file=StatDataGenome, what=character(0), sep="\t",
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


#Index GenomeGenes and Define Query Start and Query End

QueryStartPos <- GenomePos[as.logical(match(ID, QueryStartID,
                                             nomatch=0))];
QueryEndPos <- GenomePos[as.logical(match(ID, QueryEndID,
                                             nomatch=0))];

rm(list=ls(pat="^x"));

###############################################################################
#                               MAIN BODY                                     #
###############################################################################
#Make Clusters Using a Dynamic Windowing Algorithm

ClStartID <- "blank";
ClEndID <- "blank";
ClStartPos <- 0;
ClEndPos <- 0;

ClPermProb <- 0;

ClLength <- 0;

DegreesOfFreedom <- as.numeric(StatLambda)-2;

if (Crawl) {
    for (i in QueryStartPos:(QueryEndPos-1)) {
        if (tStat[i] > MintStat) {
            Search <- 1;
            SearchSize <- GapLimit;
            while ((Search <= SearchSize) && ((Search+i) <= QueryEndPos)) { 
                if (tStat[(i+Search)] > MintStat) {
                    tStatSet <- sum(tStat[i:(i+Search)]);

                    b <- 0;
                    NumberOfProbes <- Search + 1;
                    for (j in 1:TotalPermutations
                         ) {
                        PertStatSet <- sum(sample(tStat, NumberOfProbes));
                        if (PertStatSet > tStatSet) {b <- b+1}  
                    }

                    PermProb <- b/TotalPermutations;
                    SearchSize <- SearchSize + 1;

        
                    ClStartPos <- c(ClStartPos, i);
                    ClEndPos <- c(ClEndPos, i + Search);
    
                    ClStartID <- c(ClStartID, ID[i]);
                    ClEndID <- c(ClEndID, ID[i+Search]);
            
                    ClPermProb <- c(ClPermProb, PermProb)
                    ClLength <- c(ClLength, Search + 1);
                }
                Search <- Search + 1;
            }
        }
    }

    ClStartID <- ClStartID[2:length(ClStartID)];
    ClEndID <- ClEndID[2:length(ClEndID)];
    ClStartPos <- ClStartPos[2:length(ClStartPos)];
    ClEndPos <- ClEndPos[2:length(ClEndPos)];
    ClPermProb <- ClPermProb[2:length(ClPermProb)];
    ClLength <- ClLength[2:length(ClLength)];
} else {
    ClStartID <- ID[QueryStartPos];
    ClEndID <- ID[QueryEndPos];;
    ClStartPos <- QueryStartPos;
    ClEndPos <- QueryEndPos;
    ClPermProb <- 1;
    ClLength <- QueryEndPos-QueryStartPos+1;
}


#cbind(ClStartID, ClEndID, ClPermProb, ClLength)
########### Edited up to here
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
    }
  }
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

rPg <- rPg[2:length(rPg)];
rPKsub <- rPKsub[2:length(rPKsub)];
rPKu <- rPKu[2:length(rPKu)];
rPKmi <- rPKmi[2:length(rPKmi)];
rPKet <- rPKet[2:length(rPKet)];

PgKu <- PgKu[2:length(PgKu)];
PgIKu <- PgIKu[2:length(PgIKu)];
PgEKu <- PgEKu[2:length(PgEKu)];

PgKsub <- PgKsub[2:length(PgKsub)];

PgKmi <- PgKmi[2:length(PgKmi)];
PgKet <- PgKet[2:length(PgKet)];


pg <- 1-Pg;

pgKu <- 1-PgKu;
pgIKu <- 1-PgIKu;
pgEKu <- 1-PgEKu;

pgKsub <- 1-PgKsub;

pgKmi <- 1-PgKmi;
pgKet <- 1-PgKet;

ClusterSize <- length(Pg[which(Pg>0)]); 
Bon <- ClusterSize * (1-Pg);
Sidak <- (1-(Pg**(ClusterSize)));


p <- PKu -Pg;

ttStat <- tStat[which(tStat>0)];
tBHp <- BHpValue[which(tStat>0)];
tStatFilter <- order(ttStat);

ttStat <- ttStat[tStatFilter];
tBHp <- tBHp[tStatFilter];
tSidak <- Sidak[tStatFilter];
tp <- p[tStatFilter];
tpg <- pg[tStatFilter];

tpgKu <- pgKu[tStatFilter];
tpgKmi <- pgKmi[tStatFilter];
tpgIKu <- pgIKu[tStatFilter];
tpgKsub <- pgKsub[tStatFilter];

tPKu <- PKu[tStatFilter];
tPKmi <- PKmi[tStatFilter];



rankp <- (1:length(ttStat))/length(ttStat);
rankp <- rankp[order(rankp, decreasing=TRUE)]

postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/logpvslogTPlusMicroDoF", as.character(DegreesOfFreedom), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);


#x11(width=11, height=8);

y <- log10(tp);
x <- log10(ttStat);


plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)), ylim=c(min(y), max(y)),
     cex=2, pch=1);

y <- log10(tSidak);
x <- log10(ttStat);

points(x, y, cex=2, pch=1, col="red");


y <- log10(tBHp);
x <- log10(ttStat);

points(x, y, cex=2, pch=1, col="green");


y <- log10(rankp);
x <- log10(ttStat);

points(x, y, cex=2, pch=1, col="blue");

legend(-3, -3.5, c("Phenom", "Sidak", "BH"), cex=2, pch=1, col=c("black", "red", "green"));

dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/pvsTMicroDoF", as.character(DegreesOfFreedom), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);


#x11(width=11, height=8);

y <- tp;
x <- ttStat;


plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)), ylim=c(min(y), max(y)),
     cex=2, pch=1);

y <- tSidak;
x <- ttStat;

points(x, y, cex=2, pch=1, col="red");


y <- tBHp;
x <- ttStat;

points(x, y, cex=2, pch=1, col="green");


legend(55, 1, c("Phenom", "Sidak", "BH"), cex=2, pch=1, col=c("black", "red", "green"));

dev.off()
##########
y <- log10(rankp);
x <- log10(ttStat);

points(x, y, cex=2, pch=1, col="blue");

y <- log10(exp(-1*ttStat));
x <- log10(ttStat);

points(x, y, cex=2, pch=1, col="orange");



                       ############################
                       #####  OUTPUT RESULTS  #####
                       ############################

postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/PgPhenomMicroDoF", as.character(DegreesOfFreedom), "DotThree.eps", sep=""),  horizontal=TRUE, onefile=TRUE);


library(grDevices);

rgb.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"),
                                space = "rgb");

Prob.palette <- rgb.palette(length(p));
ProbFilter <- order(p, decreasing=FALSE);
tProb <- p[ProbFilter];
 
tPK <- PKsub[ProbFilter];
tPg <- Pg[ProbFilter];

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

#Need to Reorganize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OutPlotData <- cbind(PlotData[,1], tStat, LogFoldChange, Intensity,
                     AdjStndErr, ArrayNumber, BHpValue, BonpValue,
                     PostBonpValue, BonPairpValue,
                     PlotData[, 2:length(PlotData[1,])]);




# Plot Data
FilterHeader1 <- paste("StatData File1 = ", StatDataGenome, sep="");
FilterHeader2 <- " ";
FilterHeader3 <- c(paste("Query Start = ", as.character(QueryStartID),
                                           sep=""),
                 paste("Query End = ", as.character(QueryEndID), sep=""));
FilterHeader4 <- paste("Window Limit = ", as.character(WindowLimit),
                                           sep="");   
FilterHeader5 <- c("Stat Param:",
                   paste("Lambda = ", as.character(StatLambda), sep=""));
FilterHeader6 <- " ";
FilterHeader7 <- "Data Generated by DeleuzeCrawl.R";
FilterHeader8 <- date();

##### Need to reorganize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

