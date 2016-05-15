
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

xtStatPos <- 3;
xpValuePos  <- 5;

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

Track <- seq(QueryStartPos, QueryEndPos, 100);

if (Crawl) {
    for (i in QueryStartPos:(QueryEndPos-1)) {
        if (as.logical(match(i, Track, nomatch=0))) {
            print(paste("Crawl at Position", i, "of", QueryEndPos, date(),
                        sep=" "))
        }
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

print(paste("Finished Crawling", date(), sep=" "));



GeName <- "blank";
GePosition <- 0;

GetStat <- 0;
GepValue <- 0;
GeLogFoldChng <- 0;
GeIntensity <- 0;
GeBioSampleNumber <- 0;
GeArrayNumber <- 0;

GeClStartID <- " ";
GeClEndID <- " ";
GeClStartPos <- 0;
GeClEndPos <- 0;
GeClPermProb <- 0;
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

DegreesOfFreedom <- as.numeric(StatLambda)-2;


Track <- seq(1, length(ClStartID), 100);

for (i in 1:length(ClStartID)) {
    if (as.logical(match(i, Track, nomatch=0))) {
        print(paste("At Cluster", i, "of", length(ClStartID), date(), sep=" "))
    }

  for (j in ClStartPos[i]:ClEndPos[i]) {
    if (tStat[j] > MintStat) {

        
      ClPos <- 1:ClLength[i];
    
      ReftStat <- tStat[j];

      tKsubSet <- tStat[ClStartPos[i]:ClEndPos[i]];

      tKsubSet <- tKsubSet[which(ClPos != (j-ClStartPos[i] + 1))];
    
      Ku <- c(Ku, sum(tStat[ClStartPos[i]:ClEndPos[i]]));

      Ksub <- c(Ksub, sum(tKsubSet));

  
      Pg <- c(Pg, 1-2*(1-pt(ReftStat, DegreesOfFreedom)));

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
      GePosition <- c(GePosition, j);
      GeClLength <- c(GeClLength, ClLength[i]);
      GetStat <- c(GetStat, tStat[j]);
      GepValue <- c(GepValue, pValue[j]);
      GeLogFoldChng <- c(GeLogFoldChng, LogFoldChange[j]);
      GeIntensity <- c(GeIntensity, Intensity[j]);
      GeBioSampleNumber <- c(GeBioSampleNumber, BioSampleNumber[j]);
      GeArrayNumber <- c(GeArrayNumber, ArrayNumber[j]);

      GeClStartID <- c(GeClStartID, ClStartID[i]);
      GeClEndID <- c(GeClEndID, ClEndID[i]);
      GeClStartPos <- c(GeClStartPos, ClStartPos[i]);
      GeClEndPos <- c(GeClEndPos, ClEndPos[i]);
      GeClPermProb <- c(GeClPermProb, ClPermProb[i]);


    }
  }
}

print(paste("Finished Cluster Calculations", date(), sep=" "));

GeName <- GeName[2:length(GeName)];
GePosition <- GePosition[2:length(GePosition)];
GetStat <- GetStat[2:length(GetStat)];
GepValue <- GepValue[2:length(GepValue)];
GeLogFoldChng <- GeLogFoldChng[2:length(GeLogFoldChng)];
GeIntensity <- GeIntensity[2:length(GeIntensity)];
GeBioSampleNumber <- GeBioSampleNumber[2:length(GeBioSampleNumber)];
GeArrayNumber <- GeArrayNumber[2:length(GeArrayNumber)];
GeClLength <- GeClLength[2:length(GeClLength)];

GeClStartID <- GeClStartID[2:length(GeClStartID)];
GeClEndID <- GeClEndID[2:length(GeClEndID)];
GeClStartPos <- GeClStartPos[2:length(GeClStartPos)];
GeClEndPos <- GeClEndPos[2:length(GeClEndPos)];
GeClPermProb <- GeClPermProb[2:length(GeClPermProb)];




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



pg <- 1-Pg;

pgKu <- 1-PgKu;
pgIKu <- 1-PgIKu;
pgEKu <- 1-PgEKu;

pgKsub <- 1-PgKsub;

pgKmi <- 1-PgKmi;
pgKet <- 1-PgKet;

p <- PKu -Pg;

cbind(GeClStartID, GeClEndID, GeClStartPos, GeClEndPos, GeClPermProb, GeClLength, GeName, GetStat, GepValue, GeLogFoldChng, GeIntensity, GeBioSampleNumber, GeArrayNumber, Ku, Ksub, Pg, PKsub, PKu, PKmi, PKet, PgKu, PgIKu, PgEKu, PgKsub, PgKmi, PgKet)

########### Edited up to here

OutPlotData <- cbind(GeClStartID, GeClEndID, GeClStartPos, GeClEndPos,
                     GeClPermProb, GeClLength, GeName, GePosition, GetStat,
                     GepValue, GeLogFoldChng, GeIntensity, GeBioSampleNumber,
                     GeArrayNumber, Ku, Ksub, Pg, PKsub, PKu, PKmi, PKet,
                     PgKu, PgIKu, PgEKu, PgKsub, PgKmi, PgKet); 




# Plot Data
FilterHeader1 <- paste("StatData File1 = ", StatDataGenome, sep="");
FilterHeader2 <- paste("Crawl = ", as.character(Crawl), sep="");
FilterHeader3 <- c(paste("Query Start = ", as.character(QueryStartID),
                                           sep=""),
                 paste("Query End = ", as.character(QueryEndID), sep=""));
FilterHeader4 <- c(paste("Limit of Min tStat = ", as.character(MintStat),
                         sep=""),
                   paste("Gap Limit = ", as.character(GapLimit),
                                           sep=""));   
FilterHeader5 <- paste("Total Permutations = ",
                         as.character(TotalPermutations), sep="");
FilterHeader6 <- c("Stat Param:",
                   paste("Lambda = ", as.character(StatLambda), sep=""));
FilterHeader7 <- "Data Generated by DeleuzeCrawl.R";
FilterHeader8 <- date();


FilterHeaderTable <- cbind("ClStartID", "ClEndID", "ClStartPos", "ClEndPos",
                           "ClPermProb", "ClLength", "Gene", "Gene Pos",
                           "tStat", "pValue", "LogFoldChng", "Intensity",
                           "BioSampleNumber", "ArrayNumber", "Ku", "Ksub",
                           "Pg", "PKsub", "PKu", "PKmi", "PKet", "PgKu",
                           "PgIKu", "PgEKu", "PgKsub", "PgKmi", "PgKet") 


cat(FilterHeader1, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n");
cat(FilterHeader2, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeader3, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeader4, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeader5, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeader6, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeader7, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeader8, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);
cat(FilterHeaderTable, file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", StatLambda,
                     ".GeClGC", sep=""), sep="\n", append=TRUE);

write.table(OutPlotData, file=paste(OutputPrefix, 
                               "L", StatLambda, ".GeClGC", sep=""),
                               quote=FALSE, sep="\t", col.names=FALSE,
                               row.names=FALSE, append=TRUE);



#rm(list=ls());

