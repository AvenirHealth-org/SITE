#///////////////////////////////////////////////////
require(openxlsx);#Required Package, to manipulate excel files
####R version of constants
SSnegative =  1 #//Susceptible sero-negative
IIncubates = 2 #//Infected, incubation stage
IActive = 3 #//Infected, active syphilis
ILatentSSy = 4 #//Infected, latent
RecovEarly = 5 #//Recovery after treatment
RecovLater = 6 #//Recovery without treatment
RecovSusc  = 7 #//Infected, primary syphilis
numbStages = 7 #//Number of stages
#Diagnosed = 2 #//diagnosed/treated individuals 
#UnDiagnos = 1#//Undiagnosed/untreated individuals
#//Sexual beharviour Constants
SyphLow = 1; #//Low risk/1 partner
SyphMed = 2 #//Medium risk
SyphHig =3  #//High risk
SyphNoSex =5  #//Men IDU
SyphNoSexW =4  #//Men IDU
#/*Men specific constants*/
SyphMSM =4 #//MSM
#/*WoMen specific constants*/
#SyphFSW =4 #//Female sex workers
SyphNumRisksMen =5 #//Number of risk groups for men
SyphNumRisksWom =4 #//Number of risk groups for women

#Prevalence output
RPRnegTPHApos=1 #// RPR- TPHA+
RPRposTPHApos=2 #// RPR+ TPHA+
TPHApos=3 #//TPHA+
RPRpos=4  #//RPR+
numPrevTypes=4#//Prevalence depending on the test?

#Auxiliary functions
#' @description Calculates the inverse of the logit transformation
#' @param x Vector of double
#' @return Vector of double, same size as x
expit <- function(x)
{
  result<- numeric();
  result[x<=0] = exp(x[x<=0])/(1+exp(x[x<=0]))
  result[x>0] = 1/(1+exp(-x[x>0]))
  return(result)
}

#' @description Calculates the logit for doubles in the range (0,1) 
#' @param x Vector of double
#' @return Vector of double, same size as x
logit <- function(x)
{
  result <- NaN;
  result[x>0 & x<1] <- log(x[x>0 & x<1]/(1-x[x>0 & x<1]))
  return(result)
}

#'@description  Estimates prevalence at equilibrium, given incidence
#'@param ValinitInci scalar, initial incidence (in the range 0 to 1)
#'@param nu treatment seeking probability   
#'@param tau scalar, Proportions of Primary+Secondary Syphilis that are symptomatic
#'@param phi scalar, Proportion of RPR Sero-negative immediately after treament 
#'@param psi Effectiveness of treatment scalar in the range (0,1)  
#'@param sigma1 scalar, 1/(Duration Untreated Primary+Secondary Syphilis) 
#'@param sigma2 scalar, 1/(Duration Untreated Primary+Secondary Syphilis) 
#'@param sigma3 scalar, 1/(Duration Latent Syphilis (Untreated)) 
#'@param sigma4 scalar, 1/(Duration Recovered after treatment (still RPR&TPHA-seropositive))
#'@param sigma5 scalar, 1/(Duration Incidentially cured recovered (still RPR&TPHA-seropositive))
#'@param screen_rate scalar, screening rate for TPHA+, RPR+
#'@param screen_rateTPHA scalar, screening rate for TPHA+
#'@return scalar, prevalence at the equilibrium       
calcRPRPrevR <- function(ValinitInci, nu, tau, phi, psi, sigma1, sigma2, sigma3, sigma4,
            sigma5, screen_rate, screen_rateTPHA)
{
  result <- .Call("calcRPRPrevCpp", as.double(ValinitInci), as.double(nu), as.double(tau), as.double(phi),
                 as.double(psi), as.double(sigma1), as.double(sigma2), as.double(sigma3), as.double(sigma4),
                 as.double(sigma5), as.double(screen_rate), as.double(screen_rateTPHA));
  return(result)
}

#'@description  Estimates incidence at equilibrium, given prevalence
#'@param previnit scalar, initial prevalence (in the range 0 to 1)
#'@param nu treatment seeking probability   
#'@param tau scalar, Proportions of Primary+Secondary Syphilis that are symptomatic
#'@param phi scalar, Proportion of RPR Sero-negative immediately after treament 
#'@param psi Effectiveness of treatment scalar in the range (0,1)  
#'@param sigma1 scalar, 1/(Duration Untreated Primary+Secondary Syphilis) 
#'@param sigma2 scalar, 1/(Duration Untreated Primary+Secondary Syphilis) 
#'@param sigma3 scalar, 1/(Duration Latent Syphilis (Untreated)) 
#'@param sigma4 scalar, 1/(Duration Recovered after treatment (still RPR&TPHA-seropositive))
#'@param sigma5 scalar, 1/(Duration Incidentially cured recovered (still RPR&TPHA-seropositive))
#'@param screen_rate scalar, screening rate for TPHA+, RPR+
#'@param screen_rateTPHA scalar, screening rate for TPHA+
#'@return scalar, prevalence at the equilibrium   
getinitinc <- function(previnit,nu, tau, phi, psi, sigma1, sigma2, sigma3, sigma4,
                       sigma5, screen_rate, screen_rateTPHA)
{
  fff <- function(x) sapply(x,function(yy) calcRPRPrevR(ValinitInci=yy, nu, tau, phi, psi, sigma1, sigma2, sigma3, sigma4,
                                                                    sigma5, screen_rate, screen_rateTPHA)-previnit)
  res <- vector("list",1);
  names(res) <- "root";
 if((previnit>0) & (previnit<1))
 {
   x0 = 1e-10;
   fx0 = fff(x0);
   x1 = x0;
   fx1 = fx0;
   step=0.1;
   finit=FALSE
   while(1)
   {
     while(1)
     {
       x0 = x1;
       fx0 = fx1;
       x1 = x0 + step;
       if(x1>=1) 
       { 
         break
       }
       fx1 = fff(x1);
       if(fx1*fx0<0)
       {
         finit=TRUE;
         break
       }
     }
     
     if(finit | step<1e-4) break
     step = step*step;
   }
   
   if(finit)
   {
     res <- uniroot(fff, interval=c(x0,x1),tol=1e-8)
   } else
   {
     length(res)=0;
   }
     
   
   if(length(res)==0)
   {
     res <- NA
     #res <- previnit/(sigma2*nu*phi*psi+sigma4*nu*(1-phi)+
     #                   sigma3*(1-nu)+sigma5*(1-nu))
   } else
   {
     res <- res$root[which.min(res$root)[1]]
   }
 }
 return(res)
}

#Default values for the natural history parameters
NatHistParameters <- list(Duration1=4/52,Duration2=26/52,Duration3=692/52, Duration4=23/52, 
                          Duration5=55/52, nu=0.21,psi=0.99,phi=0.4, FtoMtp=0.15, MtoFtp=0.25, MtoMtp=0.25, CondomEfficacy = 0.8)

NatHistParameters <- vector("list",12L)
names(NatHistParameters) <- c( "Duration1","Duration2","Duration3","Duration4","Duration5","nu","psi","phi","FtoMtp","MtoFtp","MtoMtp","CondomEfficacy")

#Default parameters will be disabled so that the code doesn't run if wrong input file is passed to R
varNatHistParameters <- list(Duration1=0, Duration2=4/52,Duration3=150/52, Duration4=8/52, 
                          Duration5=16/52, nu=0.1,psi=0.005,phi=0.1, FtoMtp=0.05, MtoFtp=0.1, MtoMtp=0.1, CondomEfficacy = 0)
varNatHistParameters <- NatHistParameters

nameDist <- c("NA","Gamma", "Gamma", "Gamma", "Gamma","Beta","Beta","Beta","Beta","Beta","Beta","Beta")
PDistNatHistParameters <- list(Mean=NatHistParameters, Var= varNatHistParameters, 
                               Dist=nameDist)

#Default, for testing are disabled in the final version
initInci <- c(0,0,1e-4,2e-4,0,0,0,1e-4,0)/(NatHistParameters$Duration2*NatHistParameters$nu*NatHistParameters$phi*NatHistParameters$psi+
                                             NatHistParameters$Duration4*NatHistParameters$nu*(1-NatHistParameters$phi)+
                                             NatHistParameters$Duration3*(1-NatHistParameters$nu)+NatHistParameters$Duration5*(1-NatHistParameters$nu))

initInci <- rep(NA,9L);

SimTimes=seq(1970,2050, by=1);
nyear=length(SimTimes);
MenTotal = rep(NA,nyear)
WoMenTotal = rep(NA,nyear)

MenDist = matrix(NA,nrow=SyphNumRisksMen,ncol=nyear)
MenCondomUse = matrix(NA,nrow=SyphNumRisksMen,ncol=nyear)
MenScreeningRate = matrix(NA,nrow=SyphNumRisksMen,ncol=nyear)
MenMeanNumPart = matrix(NA,nrow=SyphNumRisksMen,ncol=nyear)
MenMeanNumSex = matrix(NA,nrow=SyphNumRisksMen,ncol=nyear)

WomenDist = matrix(NA,nrow=SyphNumRisksWom,ncol=nyear)
WomenCondomUse = matrix(NA,nrow=SyphNumRisksWom,ncol=nyear)
WomenScreeningRate = matrix(NA,nrow=SyphNumRisksWom,ncol=nyear)
WomenMeanNumPart = matrix(NA,nrow=SyphNumRisksWom,ncol=nyear)
WomenMeanNumSex = matrix(NA,nrow=SyphNumRisksWom,ncol=nyear)

MenTotal[1] = 1000;
WoMenTotal[1] = MenTotal[1]/0.49*0.51 # Assume that 51% are women

#Replace with NAs
WoMenTotal[] <- MenTotal[] <- NA

growthrate <-  0.01;
growthrate <- NA;

for(t in 2:nyear)
{
  MenTotal[t] = MenTotal[t-1]*(1+growthrate);
  WoMenTotal[t] = WoMenTotal[t-1]*(1+growthrate);
}

MenDist[SyphLow,]=0.50;
MenDist[SyphMed,]=0.35;
MenDist[SyphHig,]=0.10;
MenDist[SyphMSM,]=0.05;
MenDist[SyphNoSex,]=0.0;
MenDist[] <- NA

MenCondomUse[SyphLow,] = 0.1;
MenCondomUse[SyphMed,] = 0.5;
MenCondomUse[SyphHig,] = 0.5;
MenCondomUse[SyphMSM,] = 0.1;
MenCondomUse[SyphNoSex,] = 0.0;
MenCondomUse[] <- NA

MenScreeningRate[SyphLow,] = 0.0;
MenScreeningRate[SyphMed,] = 0.0;
MenScreeningRate[SyphHig,] = 0.0;
MenScreeningRate[SyphMSM,] = 0.0;
MenScreeningRate[SyphNoSex,] = 0.0;
MenScreeningRate[] <- NA

MenMeanNumPart[SyphLow,] = 1;
MenMeanNumPart[SyphMed,] = 3;
MenMeanNumPart[SyphHig,] = 20;
MenMeanNumPart[SyphMSM,] = 5;
MenMeanNumPart[SyphNoSex,] = 0.0;
MenMeanNumPart[] <- NA

MenMeanNumSex[SyphLow,] = 156;
MenMeanNumSex[SyphMed,] = 156/3;
MenMeanNumSex[SyphHig,] = 156/20;
MenMeanNumSex[SyphMSM,] = 156/5;
MenMeanNumSex[SyphNoSex,] = 0.0;
MenMeanNumSex[] <- NA

WomenDist[SyphLow,]=0.50;
WomenDist[SyphMed,]=0.45;
WomenDist[SyphHig,]=0.05;
#WomenDist[SyphFSW,]=0.05;
WomenDist[SyphNoSexW,]=0.0;
WomenDist[] <- NA

WomenCondomUse[SyphLow,] = 0.1;
WomenCondomUse[SyphMed,] = 0.5;
WomenCondomUse[SyphHig,] = 0.5;
#WomenCondomUse[SyphFSW,] = 0.3;
WomenCondomUse[SyphNoSexW,] = 0.0;
WomenCondomUse[] <- NA

WomenScreeningRate[SyphLow,] = 0.0;
WomenScreeningRate[SyphMed,] = 0.0;
WomenScreeningRate[SyphHig,] = 0.0;
#WomenScreeningRate[SyphFSW,] = 0.0;
WomenScreeningRate[SyphNoSexW,] = 0.0;
WomenScreeningRate[] <- NA

WomenMeanNumPart[SyphLow,] = 1;
WomenMeanNumPart[SyphMed,] = 3;
WomenMeanNumPart[SyphHig,] = 200;
#WomenMeanNumPart[SyphFSW,] = 500;
WomenMeanNumPart[SyphNoSexW,] = 0.0;
WomenMeanNumPart[] <- NA

WomenMeanNumSex[SyphLow,] = 156;
WomenMeanNumSex[SyphMed,] = 156/3;
WomenMeanNumSex[SyphHig,] = 2 #156/20;
#WomenMeanNumSex[SyphFSW,] = 2;
WomenMeanNumSex[SyphNoSexW,] = 0.0;
WomenMeanNumSex[] <- NA

#TurnOver for FSW. They return to the medium risk population
turnoverFSW = rep(8, length(MenTotal));
turnoverFSW[] <- NA

MatPartnerChoice = matrix(0,ncol=SyphNumRisksMen+SyphNumRisksWom,nrow=SyphNumRisksMen+SyphNumRisksWom)
#Default mixing for men, based on Kenya
MatPartnerChoice[SyphLow,SyphNumRisksMen+SyphLow]=1/3;
MatPartnerChoice[SyphLow,SyphNumRisksMen+SyphMed]=1/3;
MatPartnerChoice[SyphLow,SyphNumRisksMen+SyphHig]=1/3;# Assumption here: partnership formation depends on the availability i.e marital statuses among men (low risks are all married)
MatPartnerChoice[SyphMed,SyphNumRisksMen+SyphLow]=0.527;#Percentatge married among Medium risk, Kenya
MatPartnerChoice[SyphMed,SyphNumRisksMen+SyphMed]=1-MatPartnerChoice[SyphMed,SyphNumRisksMen+SyphLow];
MatPartnerChoice[SyphHig,SyphNumRisksMen+SyphLow]=0.59;# Percentage of high risk men married in Kenya
MatPartnerChoice[SyphHig,SyphNumRisksMen+SyphHig] = 1-MatPartnerChoice[SyphHig,SyphNumRisksMen+SyphLow] 
MatPartnerChoice[SyphMSM,SyphNumRisksMen+SyphLow]=0.527;#Percentage of married men among MSM in Kenya
MatPartnerChoice[SyphMSM,SyphNumRisksMen+SyphMSM] = 1-MatPartnerChoice[SyphMSM,SyphNumRisksMen+SyphLow] 
#Default mixing for women, based on Kenya data
MatPartnerChoice[SyphNumRisksMen+SyphLow,SyphLow]=1/4;
MatPartnerChoice[SyphNumRisksMen+SyphLow,SyphMed]=1/4;
MatPartnerChoice[SyphNumRisksMen+SyphLow,SyphHig]=1/4;# Assumption here: partnership formation depends on the availability i.e marital statuses among men (low risks are all married)
MatPartnerChoice[SyphNumRisksMen+SyphLow,SyphMSM]=1/4;# Assumption here: partnership formation depends on the availability i.e marital statuses among men (low risks are all married)
MatPartnerChoice[SyphNumRisksMen+SyphMed,SyphLow]=0.43;#Percentage married women among Medium risk, Kenya
MatPartnerChoice[SyphNumRisksMen+SyphMed,SyphMed]=1-MatPartnerChoice[SyphNumRisksMen+SyphMed,SyphLow];
MatPartnerChoice[SyphNumRisksMen+SyphHig,SyphLow]=0.1;# Percentage of high risk women married in Kenya
MatPartnerChoice[SyphNumRisksMen+SyphHig,SyphHig] = 1-MatPartnerChoice[SyphNumRisksMen+SyphHig,SyphLow]
MatPartnerChoice[] <- NA

MatPartnerReferal = matrix(0,ncol=SyphNumRisksMen+SyphNumRisksWom,nrow=SyphNumRisksMen+SyphNumRisksWom)
#Default referal probabilities for men,
MatPartnerReferal[SyphLow,SyphNumRisksMen+SyphLow]=1; #Only occurs in stable partnerships
MatPartnerReferal[SyphLow,SyphNumRisksMen+SyphMed]=1;
MatPartnerReferal[SyphLow,SyphNumRisksMen+SyphHig]=1;
MatPartnerReferal[SyphMed,SyphNumRisksMen+SyphLow]=0;
MatPartnerReferal[SyphMed,SyphNumRisksMen+SyphMed]=0;
MatPartnerReferal[SyphHig,SyphNumRisksMen+SyphLow]=0;
MatPartnerReferal[SyphHig,SyphNumRisksMen+SyphHig]=0;
MatPartnerReferal[SyphMSM,SyphNumRisksMen+SyphLow]=0;
MatPartnerReferal[SyphMSM,SyphNumRisksMen+SyphMSM]=1;
#Default referal probabilities for women
MatPartnerReferal[SyphNumRisksMen+SyphLow,SyphLow]=1;
MatPartnerReferal[SyphNumRisksMen+SyphLow,SyphMed]=1;
MatPartnerReferal[SyphNumRisksMen+SyphLow,SyphHig]=1;
MatPartnerReferal[SyphNumRisksMen+SyphLow,SyphMSM]=1;
MatPartnerReferal[SyphNumRisksMen+SyphMed,SyphLow]=0;
MatPartnerReferal[SyphNumRisksMen+SyphMed,SyphMed]=0;
MatPartnerReferal[SyphNumRisksMen+SyphHig,SyphLow]=0;
MatPartnerReferal[SyphNumRisksMen+SyphHig,SyphHig]=0;
MatPartnerReferal[] <- NA

#' @description Run Syphilis Projection
#' @param in_MenDist matrix, distribution of Males populations
#' @param in_MenTotal Vector, total male population as a function of time
#' @param in_MenCondomUse matrix, condom use in male populations
#' @param in_MenScreeningRate matrix, screening rates in male populations
#' @param in_MenMeanNumPart matrix, number of partners per year in male populations
#' @param in_MenMeanNumSex matrix, number of sex per partner per year in male populations
#' @param in_WomenDist matrix, distribution of Females populations
#' @param in_WomenTotal Vector, total female population as a function of time
#' @param in_WomenCondomUse matrix, condom use in female populations
#' @param in_WomenScreeningRate matrix, screening rates in female populations
#' @param in_WomenMeanNumPart matrix, number of partners per year in female populations
#' @param in_WomenMeanNumSex matrix, number of sex per partner per year in female populations
#' @param in_SimTimes vector, years when simulations are run
#' @param in_MixingMatrix mixing matrix
#' @param in_MatPartnerReferal matrix specifying how Syphilis Positive individuals refer their partners to clinic 
#' @param in_turnoverFSW vector, turn-over rates for FSW
#' @param in_NatHistParameters list, of 12 numbers with names "Duration1","Duration2","Duration3","Duration4","Duration5","nu","psi","phi","FtoMtp","MtoFtp","MtoMtp","CondomEfficacy"
#' @param in_initInci vector, of 9 numbers corresponding to the incidence in each sub-population at the beginning of the simulation
#' @return object of class SyphilisPopulation, i.e. contains the following fields: \cr 
#' Men: results for men; list of length 5 containing results by men-sub-populations (Low risk, medium risk, high risk, MSM and those who are not sexually active)\cr
#' Women: results for women; list of length 4 containing results by men-sub-populations (Low risk, medium risk, high risk and those who are not sexually active)\cr
#' MixingMat: mixing matrix during the last year of simualtion\cr
#' nathistory: input paramters\cr
#' firstyear: first year of projection\cr
#' InitMixingMat: balanced partner exchange rates at the beginning of the epidemic\cr
#' FinalMixingMat: balanced parner exchange rates at the end (during the last year of simulations) of the epidemic\cr
#' InitBalancedPartners: Balanced number of partners during the first year of simulations\cr
#' FinalBalancedPartners: Balanced number of partners during the last year of simulations
RunSyphProj <- function(in_MenDist=MenDist,in_MenTotal=MenTotal, in_MenCondomUse=MenCondomUse, in_MenScreeningRate=MenScreeningRate, in_MenMeanNumPart=MenMeanNumPart, in_MenMeanNumSex=MenMeanNumSex,
                        in_WomenDist=WomenDist, in_WomenTotal=WoMenTotal, in_WomenCondomUse=WomenCondomUse, in_WomenScreeningRate=WomenScreeningRate, in_WomenMeanNumPart=WomenMeanNumPart, in_WomenMeanNumSex=WomenMeanNumSex,
                        in_SimTimes=SimTimes, in_MixingMatrix=MatPartnerChoice, in_MatPartnerReferal=MatPartnerReferal, in_turnoverFSW=turnoverFSW,in_NatHistParameters=NatHistParameters, in_initInci=initInci)
{
  result = .Call("RunSyphProj",as.matrix(in_MenDist), as.double(in_MenTotal), as.matrix(in_MenCondomUse), as.matrix(in_MenScreeningRate), as.matrix(in_MenMeanNumPart), as.matrix(in_MenMeanNumSex),
                 as.matrix(in_WomenDist),as.double(in_WomenTotal), as.matrix(in_WomenCondomUse), as.matrix(in_WomenScreeningRate), as.matrix(in_WomenMeanNumPart), as.matrix(in_WomenMeanNumSex),
                 as.double(in_SimTimes),as.matrix(in_MixingMatrix),as.matrix(in_MatPartnerReferal),as.double(in_turnoverFSW),as.list(in_NatHistParameters),as.double(in_initInci))
  class(result) <- "SyphilisPopulation";
  return(result)
}

#' @description Reads population sizes from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of two elements containing population sizes/distrutions for men and women
#######################Data base reading###################################################
getPopulationSizes <- function(infilename)
{
  result <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(88:104),rowNames=TRUE)
  mcolnames <- paste("year_",result[1,],sep="");
  names(result) <-mcolnames 
  result <- result[-c(3,9),]
  checkna <- sum(apply(result,1,sum));
  if(is.na(checkna))
  {
    stop("NA is not allowed as sample size")
    return(NULL)
  } else
  {
    checkneg <- min(apply(result,1,min));
    if(checkneg<0)
    {
      stop("Negative values are allowed as sample size")
      return(NULL)
    } else
    {
      checkover100 <- FALSE
      interm0 <- apply(result[3:6,],2,sum);
      if(((max(interm0)<=100.1) & (min(interm0)>=99.99))|((100*max(interm0)<=100.1) & (100*min(interm0)>=99.99)))
      {
        mcolsum <- apply(result[3:6,],2,sum);
        result[3:6,] <- result[3:6,]/matrix(rep(mcolsum,nrow(result[3:6,])), nrow=nrow(result[3:6,]), byrow = TRUE);
        interm0[] <- apply(result[3:6,],2,sum)
      } else
      {
        stop("Total sub-populations of women were entered. Propotions were calculated")
        mcolsum <- apply(result[3:6,],2,sum);
        result[3:6,] <- result[3:6,]/matrix(rep(mcolsum,nrow(result[3:6,])), nrow=nrow(result[3:6,]), byrow = TRUE);
        interm0[] <- apply(result[3:6,],2,sum)
      }
      
      interm1 <- apply(result[8:12,],2,sum);
      if(((max(interm1)<=100.1) & (min(interm1)>=99.99))|((100*max(interm1)<=100.1) & (100*min(interm1)>=99.99)))
      {
        mcolsum <- apply(result[8:12,],2,sum);
        result[8:12,] <- result[8:12,]/matrix(rep(mcolsum,nrow(result[8:12,])), nrow=nrow(result[8:12,]), byrow = TRUE);
        interm1[] <- apply(result[8:12,],2,sum)
      } else
      {
        stop("Total men sub-populations were entered. Propotions were calculated")
        mcolsum <- apply(result[8:12,],2,sum);
        result[8:12,] <- result[8:12,]/matrix(rep(mcolsum,nrow(result[8:12,])), nrow=nrow(result[8:12,]), byrow = TRUE);
        interm1[] <- apply(result[8:12,],2,sum)
      }
      
      if(all(abs(interm1-1)<1e-8) & all(abs(interm0-1)<1e-8))
      {
        Men <- matrix(0,nrow=SyphNumRisksMen+2,ncol=ncol(result));
        Men[1:6,] <- data.matrix(result[c(1,8:13),][1:6,])
        Men[7,] <- unlist(result[c(1,8:13),][7,])
        Women <- matrix(0,nrow=SyphNumRisksWom+2,ncol=ncol(result));
        Women[1:5,] <- data.matrix(result[c(1,3:8),][1:5,])
        Women[6,] <- unlist(result[c(1,3:7),][6,])
        return(list(Men=Men, Women=Women))
      } else
      {
        stop("Proportions of men and women risk groups should sum to 1. Please check the data")
        return(NULL)
      }
    }
  }
}

#' @description Read sexual behaviour from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of 5 elements containing the following \cr
#' Men: list of 4 elements;  number of partners, number of sex, referal probabilities and probability to refer an partners with active Syphilis\cr
#' Women: list of 4 elements;  number of partners, number of sex, referal probabilities and probability to refer an partners with active Syphilis\cr
#' MixingMatrix: matrix, mixing between sub-populations\cr
#' durSexWork: double duration of sexwork\cr
#' ReferralProbabilities: matrix, referral probabilities by sub-population
getSexualBehaviour <- function(infilename)
{
  resultall <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(1:27),rowNames=FALSE)
  mcolnames <- paste("year_",resultall[1,],sep="");
  names(resultall) <- mcolnames 
  
  resultWomenNumPart <- rbind(resultall[c(1,3:5),-1],0)
  resultMenNumPart <- rbind(resultall[c(1,7:10),-1],0)
  
  resultWomenNumSex <- rbind(resultall[c(1,13:15),-1],0)
  resultMenNumSex <- rbind(resultall[c(1,17:20),-1],0)
  
  rownames(resultWomenNumPart) <- c("Year","LRW", "MRW", "FSW","NoSex")
  rownames(resultWomenNumSex) <- c("Year","LRW", "MRW", "FSW","NoSex")
  
  rownames(resultMenNumPart) <- c("Year","LRM", "MRM", "HRM","MSM","NoSex")
  rownames(resultMenNumSex) <- c("Year","LRM", "MRM", "HRM","MSM","NoSex")
  
  part_mixingMatrix <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(46:53), rowNames=TRUE, cols=c(1:8))
  mixingMatrix <- MatPartnerChoice;
  colnames(mixingMatrix) <- c("LowRiskMen", "MediumRiskMen","HighRiskMen","MSM","NoRiskMen", "LowRiskWomen","MediumRiskWomen",    
                              "HighRiskWomen_FSW", "NoRiskWomen")
  rownames(mixingMatrix) <- colnames(mixingMatrix)
  mixingMatrix[] = 0;
  mixingMatrix[c(1:SyphMSM, SyphNumRisksMen+1:SyphHig),c(1:SyphMSM, SyphNumRisksMen+1:SyphHig)] = data.matrix(part_mixingMatrix)
  mixingMatrix[SyphNoSex,SyphNoSex] = 1; 
  mixingMatrix[SyphNumRisksMen+SyphNoSexW,SyphNumRisksMen+SyphNoSexW] = 1; 
  
  part_referral <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(61:68), rowNames=TRUE, cols=c(1:8))
  ReferralProbabilities <- 0*MatPartnerChoice;
  colnames(ReferralProbabilities) <- c("LowRiskMen", "MediumRiskMen","HighRiskMen","MSM","NoRiskMen", "LowRiskWomen","MediumRiskWomen",    
                                       "HighRiskWomen_FSW", "NoRiskWomen")
  rownames(ReferralProbabilities) <- colnames(ReferralProbabilities)
  ReferralProbabilities[c(1:SyphMSM, SyphNumRisksMen+1:SyphHig),c(1:SyphMSM, SyphNumRisksMen+1:SyphHig)] <- data.matrix(part_referral);
  
  
  all_sheet_names <- openxlsx::sheets(openxlsx::loadWorkbook(infilename))
  if(any(all_sheet_names=="ScreeningRates"))
  {
    part_referral_trend <- openxlsx::read.xlsx(infilename,sheet="ScreeningRates", rows=c(59:68), rowNames=TRUE)
  } else
  {
    nstart <- 121
    part_referral_trend <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(nstart+59:68), rowNames=TRUE)
  }
  
  
  
  
  resultWomenReferTrend <- resultWomenNumPart;
  resultWomenReferTrend[-1,] <- 0;
  resultWomenReferTrend[2:4,] <- part_referral_trend[1:3,]
  
  resultMenReferTrend <- resultMenNumPart;
  resultMenReferTrend[-1,] <- 0;
  resultMenReferTrend[2:5,] <- part_referral_trend[5:8,]
  
  part_infectious_detect <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(107:119), rowNames=TRUE)
  
  resultWomenInfPDetect <- resultWomenNumPart;
  resultWomenInfPDetect[-1,] <- 0;
  resultWomenInfPDetect[2:4,] <- part_infectious_detect[3:5,]
  
  resultMenInfPDetect <- resultMenNumPart;
  resultMenInfPDetect[-1,] <- 0;
  resultMenInfPDetect[2:5,] <- part_infectious_detect[7:10,]
  
  if(any(resultWomenInfPDetect[-1,]<0 | resultWomenInfPDetect[-1,]>1) | any(is.na(resultWomenInfPDetect[-1,])) )
  {
    stop("Please check NAs and values out of range (0,1) in propensity of identifying post-infectiousness of partners tables")
  } else if(any(resultMenInfPDetect[-1,]<0 | resultMenInfPDetect[-1,]>1) | any(is.na(resultMenInfPDetect[-1,])) )
  {
    stop("Please check NAs and values out of range (0,1) in propensity of identifying post-infectiousness of partners tables")
  }
  
  durSexWork <-  openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(56:57), rowNames=TRUE, cols=c(1:2))
  if(is.na(durSexWork[1,1]))
  {
    stop("NA is neither allowed as duration of female sex work. Please check.")
    return(NULL)
  } else if(durSexWork[1,1]<=0)
  {
    stop("The duration of sex work should be non negative. Please check.")
    return(NULL)
  }
  
  if((any(is.na(resultMenNumPart))) | (any(is.na(resultWomenNumPart))) | (any(is.na(resultMenNumSex)))| (any(is.na(resultWomenNumSex))))
  {
    stop("NA is neither allowed as number of sex acts nor partners. Please check.")
    return(NULL)
  } else
  {
    if((any(resultMenNumPart<0)) | (any(resultWomenNumPart<0)) | (any(resultMenNumSex<0))| (any(resultWomenNumSex<0)))
    {
      stop("Negative values are neither allowed as number of sex acts nor partners. Please check.")
      return(NULL)
    } else
    {
      
      if(any(mixingMatrix<0))
      {
        stop("Negative values are not allowed in the mixing matrix. Please check.")
        return(NULL)
      } else
      {
        sumRows = apply(mixingMatrix,1,sum)
        if(sum(sumRows==100)==length(sumRows))
        {
          for(ii in 1:length(sumRows))
          {
            mixingMatrix[ii,]=mixingMatrix[ii,]/100;
          }
          sumRows = apply(mixingMatrix,1,sum)
        }
        
        ReferralProbabilities[rownames(ReferralProbabilities)=="NoRiskMen",] <- 0;
        ReferralProbabilities[rownames(ReferralProbabilities)=="NoRiskWomen",] <- 0;
        ReferralProbabilities[,colnames(ReferralProbabilities)=="NoRiskMen"] <- 0;
        ReferralProbabilities[,colnames(ReferralProbabilities)=="NoRiskWomen"] <- 0;
        if(any(ReferralProbabilities>1))
        {
          ReferralProbabilities/100;
        }
        
        if(any(ReferralProbabilities<0) | any(ReferralProbabilities>1))
        {
          stop("Referral probabilities should be in the rane 0 to 1")
          return(NULL)
        } else if(any(is.na(ReferralProbabilities)))
        {
          stop("NA is not allow as referral probabilities. Please check.")
          return(NULL)
        } 
        
        if(any(is.na(resultWomenReferTrend))|any(is.na(resultMenReferTrend)))
        {
          stop("NA is not allow as referral trends. Please check.")
          return(NULL)
        } else if(any(resultWomenReferTrend[-1,]<0) | any(resultWomenReferTrend[-1,]>1) | any(resultMenReferTrend[-1,]<0) | any(resultMenReferTrend[-1,]>1) )
        {
          stop("Referral trends values should be in the range 0 to 1. Please check.")
          return(NULL)
        }
        
        return(list(Men=list(NumPart=resultMenNumPart, NumSex=resultMenNumSex, ReferralTrends=resultMenReferTrend,
                             PInfectDetect=resultMenInfPDetect), 
                    Women=list(NumPart=resultWomenNumPart,NumSex=resultWomenNumSex,ReferralTrends=resultWomenReferTrend,
                               PInfectDetect=resultMenInfPDetect), 
                    MixingMatrix=mixingMatrix, durSexWork=durSexWork[1,1],
                    ReferralProbabilities=ReferralProbabilities))
        
      }
    }
  }
}

#' @description Read screening rates from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of 6 elements containing the following \cr
#' Men: matrix, screening rates by male sub-populations as a function of time\cr
#' Women: matrix, screening rates by female sub-populations as a function of time\cr
#' Men_TPHA: matrix, screening rates of TPHA only by male sub-populations as a function of time\cr
#' Women_TPHA: matrix, screening rates of TPHA only by female sub-populations as a function of time\cr
#' Men_PropSympTreat: matrix, proportions of asymtomatic male that are treated by sub-populations as a function of time\cr
#' Women_PropSympTreat: matrix, proportions of asymtomatic female that are treated by sub-populations as a function of time\cr
getScreeningRates_init <- function(infilename)
{
  resTPHA <- openxlsx::read.xlsx(infilename,sheet="ScreeningRates", rows=c(17:28),rowNames=TRUE)
  Men_TPHA <- resTPHA[6:10,]
  rownames(Men_TPHA) <- c("Low Risk Men","Medium Risk Men", "High Risk Men", "MSM","No Risk")
  
  Women_TPHA <- resTPHA[1:4,]
  rownames(Women_TPHA) <- c("Low Risk Women","Medium Risk Women", "High Risk Women", "No Risk")
  
  checknaTPHA <- is.na(sum(apply(Men_TPHA,1,sum))+sum(apply(Women_TPHA,1,sum)));
  checknegTPHA <-  any(Men_TPHA<0) | any(Women_TPHA<0);
  
  #Proportion of symptomatic cases that are treated
  resPropSympTreat <- openxlsx::read.xlsx(infilename,sheet="ScreeningRates", rows=c(34:42),rowNames=TRUE,colNames=FALSE)
  Men_PropSympTreat <- rbind(resPropSympTreat[5:8,],0)
  rownames(Men_PropSympTreat) <- c("Low Risk Men","Medium Risk Men", "High Risk Men", "MSM","No Risk")
  
  Women_PropSympTreat <- rbind(resPropSympTreat[1:3,],0)
  rownames(Women_PropSympTreat) <- c("Low Risk Women","Medium Risk Women", "High Risk Women", "No Risk")
  
  checknaPropSympTreat <- is.na(sum(apply(Men_PropSympTreat,1,sum))+sum(apply(Women_PropSympTreat,1,sum)));
  checknegPropSympTreat <- any(Men_PropSympTreat<0 | Men_PropSympTreat>1 ) | any(Women_PropSympTreat<0 | Women_PropSympTreat>1) ;
  
  #########################################################################################################
  result <- openxlsx::read.xlsx(infilename,sheet="ScreeningRates", rows=c(1:15),rowNames=TRUE)
  mcolnames <- paste("year_",result[1,],sep="");
  
  names(Men_TPHA) <- mcolnames;
  names(Women_TPHA) <- mcolnames;
  names(Men_PropSympTreat) <- mcolnames;
  names(Women_PropSympTreat) <- mcolnames;
  
  names(result) <-mcolnames 
  result <- result[-c(2:3,8),]
  checkna <- sum(apply(result,1,sum));
  if(is.na(checkna) | checknaTPHA | checknaPropSympTreat)
  {
    stop("NA is not allowed as screening rates or probabilities.")
    return(NULL)
  } else
  {
    checkneg <- min(apply(result,1,min));
    if(checkneg<0 | checknegTPHA)
    {
      stop("Negative values are allowed as screening rates")
      return(NULL)
    } else
    {
      if(checknegPropSympTreat)
      {
        stop("Probabilities should be in the rage 0 to 1")
      } else
      {
        Men <- result[c(1,6:10),]
        rownames(Men) <- c(rownames(result[c(1,6:9),]),"No Risk")
        
        Women <- result[c(1,2:5),]
        rownames(Women) <- c(rownames(result[c(1,2:4),]),"No Risk")
        
        return(list(Men=Men,Women= Women, Men_TPHA=Men_TPHA, Women_TPHA=Women_TPHA, Men_PropSympTreat=Men_PropSympTreat,Women_PropSympTreat=Women_PropSympTreat))
      }
    }
  }
}


#' @description Read screening rates from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of 6 elements containing the following \cr
#' Men: matrix, screening rates by male sub-populations as a function of time\cr
#' Women: matrix, screening rates by female sub-populations as a function of time\cr
#' Men_TPHA: matrix, screening rates of TPHA only by male sub-populations as a function of time\cr
#' Women_TPHA: matrix, screening rates of TPHA only by female sub-populations as a function of time\cr
#' Men_PropSympTreat: matrix, proportions of asymtomatic male that are treated by sub-populations as a function of time\cr
#' Women_PropSympTreat: matrix, proportions of asymtomatic female that are treated by sub-populations as a function of time\cr
getScreeningRates <- function(infilename)
{
  all_sheet_names <- openxlsx::sheets(openxlsx::loadWorkbook(infilename))
  if(any(all_sheet_names=="ScreeningRates"))
  {
    return(getScreeningRates_init(infilename))
  } else
  {
    nsart <- 121;
    resTPHA <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(nsart+17:28),rowNames=TRUE)
    Men_TPHA <- resTPHA[6:10,]
    rownames(Men_TPHA) <- c("Low Risk Men","Medium Risk Men", "High Risk Men", "MSM","No Risk")
    
    Women_TPHA <- resTPHA[1:4,]
    rownames(Women_TPHA) <- c("Low Risk Women","Medium Risk Women", "High Risk Women", "No Risk")
    
    checknaTPHA <- is.na(sum(apply(Men_TPHA,1,sum))+sum(apply(Women_TPHA,1,sum)));
    checknegTPHA <-  any(Men_TPHA<0) | any(Women_TPHA<0);
    
    #Proportion of symptomatic cases that are treated
    resPropSympTreat <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(nsart+34:42),rowNames=TRUE,colNames=FALSE)
    Men_PropSympTreat <- rbind(resPropSympTreat[5:8,],0)
    rownames(Men_PropSympTreat) <- c("Low Risk Men","Medium Risk Men", "High Risk Men", "MSM","No Risk")
    
    Women_PropSympTreat <- rbind(resPropSympTreat[1:3,],0)
    rownames(Women_PropSympTreat) <- c("Low Risk Women","Medium Risk Women", "High Risk Women", "No Risk")
    
    checknaPropSympTreat <- is.na(sum(apply(Men_PropSympTreat,1,sum))+sum(apply(Women_PropSympTreat,1,sum)));
    checknegPropSympTreat <- any(Men_PropSympTreat<0 | Men_PropSympTreat>1 ) | any(Women_PropSympTreat<0 | Women_PropSympTreat>1) ;
    
    #########################################################################################################
    result <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(nsart+1:15),rowNames=TRUE)
    mcolnames <- paste("year_",result[1,],sep="");
    
    names(Men_TPHA) <- mcolnames;
    names(Women_TPHA) <- mcolnames;
    names(Men_PropSympTreat) <- mcolnames;
    names(Women_PropSympTreat) <- mcolnames;
    
    names(result) <-mcolnames 
    result <- result[-c(2:3,8),]
    checkna <- sum(apply(result,1,sum));
    if(is.na(checkna) | checknaTPHA | checknaPropSympTreat)
    {
      stop("NA is not allowed as screening rates or probabilities.")
      return(NULL)
    } else
    {
      checkneg <- min(apply(result,1,min));
      if(checkneg<0 | checknegTPHA)
      {
        stop("Negative values are allowed as screening rates")
        return(NULL)
      } else
      {
        if(checknegPropSympTreat)
        {
          stop("Probabilities should be in the rage 0 to 1")
        } else
        {
          Men <- result[c(1,6:10),]
          rownames(Men) <- c(rownames(result[c(1,6:9),]),"No Risk")
          
          Women <- result[c(1,2:5),]
          rownames(Women) <- c(rownames(result[c(1,2:4),]),"No Risk")
          
          return(list(Men=Men,Women= Women, Men_TPHA=Men_TPHA, Women_TPHA=Women_TPHA, Men_PropSympTreat=Men_PropSympTreat,Women_PropSympTreat=Women_PropSympTreat))
        }
      }
    }
  }
}


#' @description Read Condom use from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of 2 elements containing the following \cr
#' Men: matrix, condom use by male sub-populations as a function of time\cr
#' Women: matrix, condom use by female sub-populations as a function of time
getCondomUse <- function(infilename)
{
  result <- openxlsx::read.xlsx(infilename,sheet="RReadData", rows=c(72:85),rowNames=TRUE)
  mcolnames <- paste("year_",result[1,],sep="");
  names(result) <-mcolnames 
  result <- result[-c(2,6),]
  checkna <- sum(apply(result,1,sum));
  if(is.na(checkna))
  {
    stop("NA is not allowed as condom usage")
    return(NULL)
  } else
  {
    checkneg <- min(apply(result,1,min));
    if(checkneg<0)
    {
      stop("Negative values are allowed as condom usage")
      return(NULL)
    } else
    {
      checkmax <- max(apply(result[-1,],1,max));
      if(checkmax>1)
      {
        result[-1,] <- result[-1,]/100
        checkmax <- max(apply(result[-1,],1,max));
      }
       
      if(checkmax>1)
      {
        stop("Please check condom usage")
        return(NULL)
      } else
      {
        Men <- rbind(result[c(1,5:8),],0)
        rownames(Men) <- c(rownames(result[c(1,5:8),]),"No Risk")
        
        Women <- rbind(result[c(1,2:4),],0)
        rownames(Women) <- c(rownames(result[c(1,2:4),]),"No Risk")
        return(list(Men=Men, Women=Women))
      }
    }
  }
}

#' @description Read Natural history parameters from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of 4 elements containing the following \cr
#' NatHistParam: list of 8 elements representing duration in compartments and probability to change compartments\cr
#' InitInci: vector of double. Incidence in each sub-population at during the first year of simulation\cr
#' PerSexTPs: vector of double. Male to female and female to male per sex transmission probabilities\cr
#' SexMixParams: vector of double; weith and assortativity parameters\cr
getNaHistParam <- function(infilename)
{
  PrevDataFrame <- openxlsx::read.xlsx(infilename,sheet="NatHistoryParameters", rows=c(10:17),rowNames=TRUE)
  InitInci <- NULL
  if(sum(is.na(PrevDataFrame))>=1)
  {
    stop("NA values are not allowed as parameters")
    return(NULL)
  } else
  {
    if(sum(PrevDataFrame[,1]<0)>=1)
    {
      stop("Negative values are not allowed as parameters")
      return(NULL)
    } else
    {
      if(sum(PrevDataFrame[,1]>1)>=1)
      {
        stop("nu, psi and phi should be in the range 0 to 1")
        return(NULL)
      } else
      {
        InitInci <- unlist(PrevDataFrame[,1])
        InitInci <- c(InitInci[4:7],0, InitInci[1:3],0)
      }
    }
  }
  
  PSeaTPDataFrame <- openxlsx::read.xlsx(infilename,sheet="NatHistoryParameters", rows=c(18:22),rowNames=TRUE)
  PerSexTPs <- NULL
  if(sum(is.na(PSeaTPDataFrame))>=1)
  {
    stop("NA values are not allowed as parameters")
    return(NULL)
  } else
  {
    if(any(PSeaTPDataFrame[,1]<0))
    {
      stop("Negative values are not allowed as parameters")
      return(NULL)
    } else
    {
      if(any(PSeaTPDataFrame[,1]>1))
      {
        stop("Per-sex-act transmission probabilities should be in the range 0 to 1")
        return(NULL)
      } else
      {
        PerSexTPs <- unlist(PSeaTPDataFrame[,1])
      }
    }
  }
  
  SexMixDataFrame <- openxlsx::read.xlsx(infilename,sheet="NatHistoryParameters", rows=c(23:26),rowNames=TRUE)
  SexMixParams <- NULL
  if(sum(is.na(SexMixDataFrame))>=1)
  {
    stop("NA values are not allowed as parameters")
    return(NULL)
  } else
  {
    if(any(SexMixDataFrame[,1]<0))
    {
      stop("Negative values are not allowed as parameters")
      return(NULL)
    } else
    {
      if(any(SexMixDataFrame[,1]>1))
      {
        stop("Both weight and assortativity parameters should be in the range 0 to 1")
        return(NULL)
      } else
      {
        SexMixParams <- unlist(SexMixDataFrame[,1])
        names(SexMixParams) <- c("w","assortativity","MRtoHR_assortativ")
      }
    }
  }
  
  resultDataFrame <- openxlsx::read.xlsx(infilename,sheet="NatHistoryParameters", rows=c(1:9),rowNames=TRUE)
  if(sum(is.na(resultDataFrame))>=1)
  {
    stop("NA values are not allowed as parameters")
    return(NULL)
  } else
  {
    if(sum(resultDataFrame[,1]<0)>=1)
    {
      stop("Negative values are not allowed as parameters")
      return(NULL)
    } else
    {
      if(sum(resultDataFrame[6:8,1]>1)>=1)
      {
        stop("nu, psi and phi should be in the range 0 to 1")
        return(NULL)
      } else
      {
        NatHistParam <- list(Duration1 = resultDataFrame[1,1]/52, Duration2 = resultDataFrame[2,1]/52, Duration3 = resultDataFrame[3,1]/52, 
                             Duration4 = resultDataFrame[4,1]/52, Duration5 = resultDataFrame[5,1]/52, nu = resultDataFrame[6,1],
                       psi=resultDataFrame[7,1], phi = resultDataFrame[8,1])
        result <- list(NatHistParam=NatHistParam,InitInci=InitInci,PerSexTPs=PerSexTPs, SexMixParams=SexMixParams)
        return(result)
      }
    }
  }
} #end getNaHistParam

#Get Prevalence data
#' @description Read observed prevalence in surveys from the input file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list 
getSyphPrevFromDb <- function(infilename)
{
  names_all <- c("Women + Men ","Women","LowRiskWomen","MediumRiskWomen","HighRiskWomen","Men","LowRiskMen","MediumRiskMen","HighRiskMen","MSM")
  
  res_all <- vector("list",10)
  countdeb = 0; 
  
  for(ii in 1:10)
  {
    resultBS <- openxlsx::read.xlsx(infilename,sheet="PrevalenceData", rows=c(countdeb+1:5),rowNames=TRUE, colNames=(countdeb<=0))
    
    if(ncol(resultBS)>0)
    {
      mcolnames <- paste("year_",resultBS[1,],sep="");
      names(resultBS) <-mcolnames 
      checkna <- sum(apply(resultBS,1,sum));
      if(is.na(checkna))
      {
        stop("NA is not allowed as Prevalence data. The corresponding survey(s) will be ignored.")
        resultBS <- resultBS[,which(!is.na(checkna))]
      } else
      {
        checkneg <- apply(resultBS,2,min);
        if(any(checkneg<0))
        {
          stop("Negative values as Prevalence data. The corresponding survey(s) will be ignored.")
          resultBS <- resultBS[,which(checkneg>0)]
        } 
        
        if(any(resultBS[2,]>1))
        {
          stop("Prevalence should be in the range zero to 1. Survey(s) will be ignored.")
          resultBS <- resultBS[,which(checkneg>0)]
        }
      }
      
      res_all[[ii]] <- resultBS;
    }
    
    countdeb = countdeb+5
  }
  names(res_all) <- names_all
  res_all
}

#Get Prevalence data
#' @description Read all input data from the file
#' @param infilename string, name of the xlsx file where the data have been saved
#' @return list of class "Syphdata i.e. containing the following fields:\cr
#' MenDist: matrix, men distribution as a function of time\cr 
#' MenTotal: vector total number of men as a function of time\cr
#' MenCondomUse: matrix, condom use among men as a function of time, by sub-population\cr
#' MenScreeningRate: list, screening rates among men as a function of time, by sub-population\cr
#' MenMeanNumPart: matrix, number of partners among men, by sub-populations\cr
#' MenMeanNumSex: number of sex per partner for men, by sub-population\cr 
#' WomenDist: matrix, women distribution as a function of time\cr 
#' WoMenTotal: vector total number of women as a function of time\cr
#' WomenCondomUse: matrix, condom use among women as a function of time, by sub-population\cr
#' WomenScreeningRate: list, screening rates among women as a function of time, by sub-population\cr
#' WomenMeanNumPart: matrix, number of partners among women, by sub-populations\cr
#' WomenMeanNumSex: number of sex per partner for women, by sub-population\cr 
#' SimTimes: vector, years when simulations are run\cr
#' MatPartnerChoice: mixing matrix\cr
#' ReferralProbabilities: matrix, referral probabilities by sub-population\cr
#' turnoverFSW: vector, turn-over rates for FSW\cr
#' NatHistParam: list, of 12 numbers with names "Duration1","Duration2","Duration3","Duration4","Duration5","nu","psi","phi","FtoMtp","MtoFtp","MtoMtp","CondomEfficacy"\cr 
#' initInci: vector, incidence in sub-populations during the first year of simulation
#' PrevData: list, prevalence data collected in surveys
getSyphDataFromDb <- function(infilename)
{
 SSes <- getPopulationSizes(infilename)
 SexBe <- getSexualBehaviour(infilename)
 ScRat <- getScreeningRates(infilename)
 Cuse <- getCondomUse(infilename)
 PrevData <- getSyphPrevFromDb(infilename)
 MenDist <- data.matrix(SSes$Men[-c(1,7),]) 
 MenTotal <- unlist(SSes$Men[7,] )
 MenCondomUse <- data.matrix(Cuse$Men[-1,])
 MenScreeningRate <- data.matrix(rbind(ScRat$Men[-1,],ScRat$Men_TPHA, 
                                       ScRat$Men_PropSympTreat,
                                       SexBe$Men$ReferralTrends[-1,],
                                       SexBe$Men$PInfectDetect[-1,]))
 MenMeanNumPart <- data.matrix(SexBe$Men$NumPart[-1,])
 MenMeanNumSex <- data.matrix(SexBe$Men$NumSex[-1,])
 
 WomenDist <- data.matrix(SSes$Women[-c(1,6),]) 
 WoMenTotal <- unlist(SSes$Women[6,] )
 WomenCondomUse <- data.matrix(Cuse$Women[-1,])
 WomenScreeningRate <- data.matrix(rbind(ScRat$Women[-1,],ScRat$Women_TPHA,
                                         ScRat$Women_PropSympTreat,SexBe$Women$ReferralTrends[-1,], 
                                         SexBe$Women$PInfectDetect[-1,]))
 WomenMeanNumPart <- data.matrix(SexBe$Women$NumPart[-1,])
 WomenMeanNumSex <- data.matrix(SexBe$Women$NumSex[-1,])
 SimTimes <- 0.5*(unlist(Cuse$Men[1,])+unlist(Cuse$Women[1,])) # Just to check the length here

 MatPartnerChoice <- SexBe$MixingMatrix;
 ReferralProbabilities <- SexBe$ReferralProbabilities;
 turnoverFSW <- rep(SexBe$durSexWork,length(SimTimes))
 
 NatHistParam <- getNaHistParam(infilename)
 NatHistParam$NatHistParam$FtoMtp <- NatHistParam$PerSexTPs[1];
 NatHistParam$NatHistParam$MtoFtp <- NatHistParam$PerSexTPs[2];
 NatHistParam$NatHistParam$MtoMtp <- NatHistParam$PerSexTPs[3];
 NatHistParam$NatHistParam$CondomEfficacy <- NatHistParam$PerSexTPs[4];
 NatHistParam$NatHistParam$w_assortativeness <- NatHistParam$SexMixParams[names(NatHistParam$SexMixParams)=="w"]
 NatHistParam$NatHistParam$assortativity <- NatHistParam$SexMixParams[names(NatHistParam$SexMixParams)=="assortativity"]
 NatHistParam$NatHistParam$MRtoHR_assortativ <- NatHistParam$SexMixParams[names(NatHistParam$SexMixParams)=="MRtoHR_assortativ"]
 
 try(
 NatHistParam$InitInci[1:4] <- sapply(1:4, function(ll)
 {
   getinitinc(previnit=NatHistParam$InitInci[ll],NatHistParam$NatHistParam$nu, MenScreeningRate[2*SyphNumRisksMen+ll,1],
              NatHistParam$NatHistParam$phi, NatHistParam$NatHistParam$psi, 1/NatHistParam$NatHistParam$Duration1,
              1/NatHistParam$NatHistParam$Duration2, 1/NatHistParam$NatHistParam$Duration3,
              1/NatHistParam$NatHistParam$Duration4, 1/NatHistParam$NatHistParam$Duration5, MenScreeningRate[ll,1],
              MenScreeningRate[SyphNumRisksMen+ll,1])
 }),silent=FALSE)
 
 try(
 NatHistParam$InitInci[5+1:3] <- sapply(1:3, function(ll)
 {
   getinitinc(previnit=NatHistParam$InitInci[5+ll],NatHistParam$NatHistParam$nu, WomenScreeningRate[2*SyphNumRisksWom+ll,1],
              NatHistParam$NatHistParam$phi, NatHistParam$NatHistParam$psi, 1/NatHistParam$NatHistParam$Duration1,
              1/NatHistParam$NatHistParam$Duration2, 1/NatHistParam$NatHistParam$Duration3,
              1/NatHistParam$NatHistParam$Duration4, 1/NatHistParam$NatHistParam$Duration5, 
              WomenScreeningRate[ll,1], WomenScreeningRate[SyphNumRisksWom+ll,1])
 }),silent=FALSE)
 
 if(any((is.na(NatHistParam$InitInci)|(NatHistParam$InitInci<0))))
 {
    stop("Initial prevalence RPR+ might not be possible. Please check")  
 }
 
 result <- list( MenDist = MenDist, MenTotal = MenTotal, MenCondomUse=MenCondomUse,
                 MenScreeningRate = MenScreeningRate, MenMeanNumPart=MenMeanNumPart,
                 MenMeanNumSex = MenMeanNumSex, WomenDist=WomenDist, WoMenTotal = WoMenTotal,
                 WomenCondomUse = WomenCondomUse, WomenScreeningRate = WomenScreeningRate,
                 WomenMeanNumPart = WomenMeanNumPart, WomenMeanNumSex = WomenMeanNumSex,
                 SimTimes = SimTimes, MatPartnerChoice = MatPartnerChoice, ReferralProbabilities=ReferralProbabilities,
                 turnoverFSW=turnoverFSW, NatHistParam=NatHistParam$NatHistParam, initInci=NatHistParam$InitInci, PrevData=PrevData)
  any_null <-  sum(unlist(lapply(result,is.null)))
  any_na <- sum(unlist(lapply(result,function(x) any(is.na(x)))))
  result$filename <- infilename
  
  if(any_na==0 & any_null==0 & length(result)!=0)
  {
    class(result) <- "SyphData"
  }
  return(result)
}

#' @description Run Syphilis simulation
#' @param list of class "Syphdata" i.e. containing the following fields:\cr
#' MenDist: matrix, men distribution as a function of time\cr 
#' MenTotal: vector total number of men as a function of time\cr
#' MenCondomUse: matrix, condom use among men as a function of time, by sub-population\cr
#' MenScreeningRate: list, screening rates among men as a function of time, by sub-population\cr
#' MenMeanNumPart: matrix, number of partners among men, by sub-populations\cr
#' MenMeanNumSex: number of sex per partner for men, by sub-population\cr 
#' WomenDist: matrix, women distribution as a function of time\cr 
#' WoMenTotal: vector total number of women as a function of time\cr
#' WomenCondomUse: matrix, condom use among women as a function of time, by sub-population\cr
#' WomenScreeningRate: list, screening rates among women as a function of time, by sub-population\cr
#' WomenMeanNumPart: matrix, number of partners among women, by sub-populations\cr
#' WomenMeanNumSex: number of sex per partner for women, by sub-population\cr 
#' SimTimes: vector, years when simulations are run\cr
#' MatPartnerChoice: mixing matrix\cr
#' ReferralProbabilities: matrix, referral probabilities by sub-population\cr
#' turnoverFSW: vector, turn-over rates for FSW\cr
#' NatHistParam: list, of 12 numbers with names "Duration1","Duration2","Duration3","Duration4","Duration5","nu","psi","phi","FtoMtp","MtoFtp","MtoMtp","CondomEfficacy"\cr 
#' initInci: vector, incidence in sub-populations during the first year of simulation
#' PrevData: list, prevalence data collected in surveys
#' @return  object of class SyphilisPopulation, i.e. contains the following fields: \cr 
#' Men: results for men; list of length 5 containing results by men-sub-populations (Low risk, medium risk, high risk, MSM and those who are not sexually active)\cr
#' Women: results for women; list of length 4 containing results by men-sub-populations (Low risk, medium risk, high risk and those who are not sexually active)\cr
#' MixingMat: mixing matrix during the last year of simualtion\cr
#' nathistory: input paramters\cr
#' firstyear: first year of projection\cr
#' InitMixingMat: balanced partner exchange rates at the beginning of the epidemic\cr
#' FinalMixingMat: balanced parner exchange rates at the end (during the last year of simulations) of the epidemic\cr
#' InitBalancedPartners: Balanced number of partners during the first year of simulations\cr
#' FinalBalancedPartners: Balanced number of partners during the last year of simulations
RunSyphProjFromSyphData <- function(x)
{
  n_odev <- dev.list();
  if(length(n_odev)>=6)
  {
    n_odev <- closealldevices() #Close all the devices if there are any devices opened
  }
  result=NULL;
  if(class(x)=="SyphData")
  {
    result <- RunSyphProj(in_MenDist=x$MenDist,in_MenTotal=x$MenTotal, in_MenCondomUse=x$MenCondomUse, 
                        in_MenScreeningRate=x$MenScreeningRate, in_MenMeanNumPart=x$MenMeanNumPart, 
                        in_MenMeanNumSex=x$MenMeanNumSex, in_WomenDist=x$WomenDist, in_WomenTotal=x$WoMenTotal, 
                        in_WomenCondomUse=x$WomenCondomUse, in_WomenScreeningRate=x$WomenScreeningRate, 
                        in_WomenMeanNumPart=x$WomenMeanNumPart, in_WomenMeanNumSex=x$WomenMeanNumSex,
                        in_SimTimes=x$SimTimes, in_MixingMatrix=x$MatPartnerChoice, in_MatPartnerReferal=x$ReferralProbabilities,in_turnoverFSW=x$turnoverFSW,
                        in_NatHistParameters=x$NatHistParam,in_initInci=x$initInci)
    namesmixgmats <- c("LRMen","MRMen","HRMen","MSM","NoRiskMen","LRWomen","MRWomen","HRWomen","NoRiskWomen")
    rownames(result$InitMixingMat$CombinedMat) <- colnames(result$InitMixingMat$CombinedMat) <- namesmixgmats
    rownames(result$FinalMixingMat$CombinedMat) <- colnames(result$FinalMixingMat$CombinedMat) <- namesmixgmats
    
    rownames(result$InitBalancedPartners$BalancedPart) <- colnames(result$InitBalancedPartners$BalancedPart) <- namesmixgmats
    rownames(result$FinalBalancedPartners$BalancedPart) <- colnames(result$FinalBalancedPartners$BalancedPart) <- namesmixgmats
    result$filename <- gsub(".xlsx","_out.xlsx",x$filename)
  } else
  {
    stop("x must be of class SyphData")
  }
  return(result);
}

########################SUMMARIZING THE RESULTS############################################
#***********************PREVALENCE*********************************************************#
#' @description Summarize the results
#' @param list of class "SyphilisPopulation"
#' @return list of length 3 containing the following fields:
#' summary: dataframe containing prevalences among both sexes\cr
#' men: dataframe containing prevalences among both men\cr
#' women: dataframe containing prevalences among both women
CalcPrev <- function(result)
{
  if(class(result)!="SyphilisPopulation")
  {
    stop("result should of type SyphilisPopulation")
  } else
  {
    ListPrevMen <- list()
    ListPrevWoMen <- list()
    ListPrevAll <- vector("list",3);
    for(ii in 1:2)
    {
      n_subpop = length(result[[ii]]);
      if(ii==1) ListPrevMen <- vector("list",n_subpop) else ListPrevWoMen <- vector("list",n_subpop)
      num_active_all <- rep(0,length(result[[ii]][[1]]$SSneg));
      num_latent_all <- num_active_all
      num_RPRnegTPHAPos_all <- num_active_all
      num_RPRTPHAPos_all <- num_active_all
      num_screened_all <- num_active_all
      num_screenedPositive_all <- num_active_all
      num_screenedTPHA_all <- num_active_all
      num_screenedTPHAPositive_all <- num_active_all;
      num_symptomatic_all <- num_active_all;
      num_symptomaticTreat_all <- num_active_all;
      
      num_referred_all <- num_active_all;
      num_referredTPHAPos_all <- num_active_all;
      num_referredRPRPos_all <- num_active_all;
      
      num_EarlyRecov_all <- num_active_all;
      num_LateRecov_all <- num_active_all;
      
      num_total_pop <- num_active_all
      for(jj in 1:n_subpop)
      {
        part_result= data.frame(year=result$firstyear-1+1:length(num_active_all));
        part_result$PrevActive <- result[[ii]][[jj]]$IActive/result[[ii]][[jj]]$Total
        part_result$PrevLatent <- result[[ii]][[jj]]$ILatent/result[[ii]][[jj]]$Total
        
        part_result$PrevActive[result[[ii]][[jj]]$Total==0] <- 0
        part_result$PrevLatent[result[[ii]][[jj]]$Total==0] <- 0
        
        part_result$PrevRPRnegTPHAPos <- result[[ii]][[jj]]$allPrevalenceTypes[RPRnegTPHApos,] 
        part_result$PrevRPRTPHAPos <- result[[ii]][[jj]]$allPrevalenceTypes[RPRposTPHApos,]
        
        part_result$NumScreened <- result[[ii]][[jj]]$NumberScreened[2,]
        part_result$NumScreenedPositive <- result[[ii]][[jj]]$NumTreatedaftScreened[2,]
        part_result$NumScreenedTPHA <- result[[ii]][[jj]]$NumberScreened[1,]
        part_result$NumScreenedTPHAPositive <- result[[ii]][[jj]]$NumTreatedaftScreened[1,]
        
        part_result$Symptomatic <- result[[ii]][[jj]]$NumSymptomaticClinTreat
        part_result$SymptomaticTreat <- result[[ii]][[jj]]$NumSymptomaticClinTreat #NOT CORRECT
        
        part_result$Referred <- result[[ii]][[jj]]$NumContactTraced[1,]
        part_result$ReferredTPHAPos <- result[[ii]][[jj]]$NumContactTraced[2,]
        part_result$ReferredRPRPos <- result[[ii]][[jj]]$NumContactTraced[3,]
        
        part_result$PrevLateRecov <- result[[ii]][[jj]]$LateRecov/result[[ii]][[jj]]$Total
        part_result$PrevEarlyRecov <- result[[ii]][[jj]]$EarlyRecov/result[[ii]][[jj]]$Total
        
        part_result$PrevLateRecov[result[[ii]][[jj]]$Total==0] <- 0
        part_result$PrevEarlyRecov[result[[ii]][[jj]]$Total==0] <- 0
        
        part_result$Total <- result[[ii]][[jj]]$Total;
        
        num_active_all <- num_active_all+result[[ii]][[jj]]$IActive;
        num_latent_all <- num_latent_all+result[[ii]][[jj]]$ILatent
        num_RPRnegTPHAPos_all <- num_RPRnegTPHAPos_all+result[[ii]][[jj]]$RecovSSneg
        num_RPRTPHAPos_all <- num_RPRTPHAPos_all+result[[ii]][[jj]]$EarlyRecov+result[[ii]][[jj]]$LateRecov+result[[ii]][[jj]]$IActive+result[[ii]][[jj]]$ILatent
        num_screened_all <- num_screened_all+result[[ii]][[jj]]$NumberScreened[2,]
        num_screenedPositive_all <- num_screenedPositive_all+result[[ii]][[jj]]$NumTreatedaftScreened[2,]
        num_screenedTPHA_all <- num_screenedTPHA_all+result[[ii]][[jj]]$NumberScreened[1,]
        num_screenedTPHAPositive_all <- num_screenedTPHAPositive_all+result[[ii]][[jj]]$NumTreatedaftScreened[1,];
        
        num_symptomatic_all <- num_symptomatic_all+result[[ii]][[jj]]$NumSymptomaticClinTreat#NA*result[[ii]][[jj]]$NumIncidentCases;
        num_symptomaticTreat_all <- num_symptomaticTreat_all+result[[ii]][[jj]]$NumSymptomaticClinTreat;
        
        num_referred_all <- num_referred_all + result[[ii]][[jj]]$NumContactTraced[1,];#result[[ii]][[jj]]$NumReferred;
        num_referredTPHAPos_all <- num_referredTPHAPos_all + result[[ii]][[jj]]$NumContactTraced[2,];#result[[ii]][[jj]]$NumReferredTestedPositive[1,];
        num_referredRPRPos_all <- num_referredRPRPos_all + result[[ii]][[jj]]$NumContactTraced[3,];#result[[ii]][[jj]]$NumReferredTestedPositive[2,];
        
        num_EarlyRecov_all <- num_EarlyRecov_all + result[[ii]][[jj]]$LateRecov;
        num_LateRecov_all <- num_LateRecov_all + result[[ii]][[jj]]$EarlyRecov;
        
        num_total_pop <- num_total_pop+result[[ii]][[jj]]$Total
        
        if(ii==1)
        {
          ListPrevMen[[jj]] <- part_result;
          names(ListPrevMen)[jj] = names(result[[ii]])[[jj]]
        } else
        {
          ListPrevWoMen[[jj]] <- part_result;
          names(ListPrevWoMen)[jj] = names(result[[ii]])[[jj]]
        }
        
        if(jj==n_subpop)
        {
          part_result= data.frame(year=result$firstyear-1+1:length(num_active_all));
          part_result$PrevActive <- num_active_all/num_total_pop
          part_result$PrevLatent <- num_latent_all/num_total_pop
          part_result$PrevRPRnegTPHAPos <- num_RPRnegTPHAPos_all/num_total_pop
          part_result$PrevRPRTPHAPos <- num_RPRTPHAPos_all /num_total_pop
          
          part_result$NumScreened <-  num_screened_all
          part_result$NumScreenedPositive <-  num_screenedPositive_all
          part_result$NumScreenedTPHA <-  num_screenedTPHA_all
          part_result$NumScreenedTPHAPositive <-  num_screenedTPHAPositive_all
          
          part_result$Symptomatic <- num_symptomatic_all;
          part_result$SymptomaticTreat <- num_symptomaticTreat_all;
          
          part_result$Referred <- num_referred_all;
          part_result$ReferredTPHAPos <- num_referredTPHAPos_all;
          part_result$ReferredRPRPos <- num_referredRPRPos_all;
          
          part_result$PrevLateRecov <- num_EarlyRecov_all/num_total_pop;
          part_result$PrevEarlyRecov <- num_LateRecov_all/num_total_pop;
          
          part_result$PrevLateRecov[num_total_pop==0] <- 0;
          part_result$PrevEarlyRecov[num_total_pop==0] <- 0;
          
          part_result$Total <- num_total_pop;
          ListPrevAll[[ii]] <- part_result
        }
      }#End n_subpop
      
    }#End for ii
    names(ListPrevAll) =c("Men", "Women")
    part_result= data.frame(year=result$firstyear-1+1:length(result[[1]][[1]][[1]]));
    total <- (ListPrevAll$Men$Total+ListPrevAll$Women$Total);
    part_result$PrevActive <- (ListPrevAll$Men$PrevActive*ListPrevAll$Men$Total+ListPrevAll$Women$PrevActive*ListPrevAll$Women$Total)/total
    part_result$PrevLatent <- (ListPrevAll$Men$PrevLatent*ListPrevAll$Men$Total+ListPrevAll$Women$PrevLatent*ListPrevAll$Women$Total)/total
    part_result$PrevRPRnegTPHAPos <- (ListPrevAll$Men$PrevRPRnegTPHAPos*ListPrevAll$Men$Total+ListPrevAll$Women$PrevRPRnegTPHAPos*ListPrevAll$Women$Total)/total
    part_result$PrevRPRTPHAPos <- (ListPrevAll$Men$PrevRPRTPHAPos*ListPrevAll$Men$Total+ListPrevAll$Women$PrevRPRTPHAPos*ListPrevAll$Women$Total)/total
    part_result$NumScreened <- ListPrevAll$Men$NumScreened+ListPrevAll$Women$NumScreened
    part_result$NumScreenedPositive <- ListPrevAll$Men$NumScreenedPositive+ListPrevAll$Women$NumScreenedPositive
    part_result$NumScreenedTPHA <- ListPrevAll$Men$NumScreenedTPHA+ListPrevAll$Women$NumScreenedTPHA
    part_result$NumScreenedTPHAPositive <- ListPrevAll$Men$NumScreenedTPHA+ListPrevAll$Women$NumScreenedTPHAPositive
    
    part_result$PrevActive[total==0] <- 0
    part_result$PrevLatent [total==0] <- 0
    part_result$PrevRPRnegTPHAPos [total==0] <- 0
    part_result$PrevRPRTPHAPos [total==0] <- 0
    
    part_result$Symptomatic <- ListPrevAll$Men$Symptomatic+ListPrevAll$Women$Symptomatic;
    part_result$SymptomaticTreat <- ListPrevAll$Men$SymptomaticTreat+ListPrevAll$Women$SymptomaticTreat;
    
    part_result$Referred <- ListPrevAll$Men$Referred+ListPrevAll$Women$Referred;
    part_result$ReferredTPHAPos <- ListPrevAll$Men$ReferredTPHAPos+ListPrevAll$Women$ReferredTPHAPos;
    part_result$ReferredRPRPos <- ListPrevAll$Men$ReferredRPRPos+ListPrevAll$Women$ReferredRPRPos;
    
    part_result$PrevLateRecov <- (ListPrevAll$Men$PrevLateRecov*ListPrevAll$Men$Total+ListPrevAll$Women$PrevLateRecov*ListPrevAll$Women$Total)/total;
    part_result$PrevEarlyRecov <- (ListPrevAll$Men$PrevEarlyRecov*ListPrevAll$Men$Total+ListPrevAll$Women$PrevEarlyRecov*ListPrevAll$Women$Total)/total;
    
    part_result$Total <- total;
    ListPrevAll[[3]] <- part_result;
    names(ListPrevAll)[3]="Women + Men ";
    
    class(ListPrevAll)="SummarySyphilisPopulation";
    class(ListPrevMen) = "SummaryMen"
    class(ListPrevWoMen) = "SummaryWoMen"
    
    result <- list(summary=ListPrevAll,men=ListPrevMen,women=ListPrevWoMen)
    class(result) = "SyphilisSummary"
    return(result)
  }
}

#' @description Display results in chart
#' @param list of class "SyphilisSummary"
myplotSyph <- function(X, fmain="Syphilis Prevalence")
{
  maxprev = max(X[,2:5],na.rm=T);
  pow=1
  if(maxprev>0)
  {
    while(maxprev<0.1)
    {
      pow=pow*10;
      maxprev=maxprev*pow;
    }
  }
  
  fylim = c(min(pow*X[,2:5],na.rm=T),max(max(pow*X[,2:5],na.rm=T),max(pow*X[,5]+pow*X[,4],na.rm=T)))
  plot(x=X$year,y=pow*X$PrevActive,xlab="Time (years)",ylab=paste("Prevalence x", pow),type="l",cex=1.2, col="blue",ylim=fylim, main=fmain)
  
  lines(X$year,pow*X$PrevLatent, col="green")
  lines(X$year,pow*X$PrevRPRnegTPHAPos+pow*X$PrevRPRTPHAPos, col="red")
  lines(X$year,pow*X$PrevRPRTPHAPos)
  legend("topright",c("Primary + Secondary Syphilis", "Latent","TPHA+","RPR+TPHA+"),col=c("blue","green","red","black"),lty=c(1,1,1,1),box.col=NA)
}

#' @description Display prevalence results in chart
#' @param result: list of class "SyphilisSummary"
#' @param pop: string, population type. Should be either "all" (default) or "Men" or "Women"
plot.SummarySyphilisPopulation <- function(result,pop="all")
{
  if(class(result)=="SummarySyphilisPopulation")
  {
    if(pop=="all")
    {
      GX = result$'Women + Men '
      myplotSyph(GX,fmain="Syphilis, Women + Men ")
    } else if(pop=="Men")
    {
      GX = result$Men
      myplotSyph(GX,fmain="Syphilis among Men")
    }else if(pop=="Women")
    {
      GX = result$Women
      myplotSyph(GX,fmain="Syphilis among Women")
    } else
    {
      windows()
      par(mfrow=c(3,1), mar=c(4.5,4.0,2.0,2.0),cex=0.8)
      GX = result$'Women + Men '
      myplotSyph(GX,fmain="Syphilis, Women + Men ")
      GX = result$Men
      myplotSyph(GX,fmain="Syphilis among Men")
      GX = result$Women
      myplotSyph(GX,fmain="Syphilis among Women")
    }
  } else
  {
    stop("class of result should be SummarySyphilisPopulation")
  }
}

#' Display prevalence results in chart
#' @param result: list of class "SyphilisSummaryMen"
#' @param pop: string, population type. Should be either "all" or men sub-populations
plot.SummarySyphilisMen<- function(result,pop="all")
{
  if(class(result)=="SummaryMen")
  {
    newpop= which(names(result)==pop)
    if(length(newpop)==1)
    {
      fname_out=paste("Syphilis among men,", pop,sep=" ")
      GX = result[[newpop]]
      myplotSyph(GX,fmain=fname_out)
    } else
    {
      windows()
      par(mfrow=c(2,2), mar=c(4.5,4.0,2,2.0),cex=0.8)
      for(ii in 1:(length(result)-1))
      {
        rg <- names(result)[[ii]];
        if(rg=="LowRisk") 
        {
          rg <- "Low Risk"
        } else if (rg=="MediumRisk") 
        {
          rg <- "Medium Risk"
        } else if (rg=="HighRisk")
        {
          rg <- "High Risk"
        }
        fname_out=paste("Syphilis among men,", rg,sep=" ")
        GX = result[[ii]]
        myplotSyph(GX,fmain=fname_out)
      }
    }
  } else
  {
    stop("class of result should be SummaryMen")
  }
}

#' @description Display results in chart
#' @param result: list of class "SyphilisSummaryWomen"
#' @param pop: string, population type. Should be either "all" or women sub-populations
plot.SummarySyphilisWomen<- function(result,pop="all")
{
  if(class(result)=="SummaryWoMen")
  {
    newpop= which(names(result)==pop)
    if(length(newpop)==1)
    {
      fname_out=paste("Syphilis among women,", pop,sep=" ")
      GX = result[[newpop]]
      myplotSyph(GX,fmain=fname_out)
    } else
    {
      windows()
      par(mfrow=c(3,1), mar=c(4.5,4.0,2.0,2.0),cex=0.8)
      for(ii in 1:(length(result)-1))
      {
        rg <- names(result)[[ii]];
        if(rg=="LowRisk") 
        {
          rg <- "Low Risk"
        } else if (rg=="MediumRisk") 
        {
          rg <- "Medium Risk"
        } else if (rg=="HighRisk")
        {
          rg <- "High Risk"
        }
        
        fname_out=paste("Syphilis among women,", rg,sep=" ")
        GX = result[[ii]]
        myplotSyph(GX,fmain=fname_out)
      }
    }
  } else
  {
    stop("class of result should be SummaryMen")
  }
}

#***********************INCIDENCE*********************************************************#
#' @description Summarize the results
#' @param result: list of class "SyphilisPopulation"
#' @return list of length 3 containing the following fields:
#' summary: dataframe containing Incidence rates among both sexes\cr
#' men: dataframe containing Incidence rates among both men\cr
#' women: dataframe containing Incidence rates among both women
CalcIncidence <- function(result)
{
  if(class(result)!="SyphilisPopulation")
  {
    stop("result should of type SyphilisPopulation")
  } else
  {
    ListInciMen <- list()
    ListInciWoMen <- list()
    ListInciAll <- vector("list",3);
    for(ii in 1:2)
    {
      n_subpop = length(result[[ii]]);
      if(ii==1) ListInciMen <- vector("list",n_subpop) else ListIncivWoMen <- vector("list",n_subpop)
      num_susceptibles_all <- rep(0,length(result[[ii]][[1]]$SSneg));
      num_newinfections_all <- num_susceptibles_all
      num_total_pop <- num_susceptibles_all
      for(jj in 1:n_subpop)
      {
        part_result= data.frame(year=result$firstyear-1+1:length(num_susceptibles_all));
        part_result$Incidence <- result[[ii]][[jj]]$Incidence
        part_result$NewCases <- result[[ii]][[jj]]$NumIncidentCases; #result[[ii]][[jj]]$Incidence*(result[[ii]][[jj]]$SSneg+result[[ii]][[jj]]$RecovSSneg)
        part_result$Susceptibles <- result[[ii]][[jj]]$SSneg+result[[ii]][[jj]]$RecovSSneg
        part_result$Total <- result[[ii]][[jj]]$Total;
        
        num_susceptibles_all <- num_susceptibles_all+(result[[ii]][[jj]]$SSneg+result[[ii]][[jj]]$RecovSSneg);
        num_newinfections_all <- num_newinfections_all+part_result$NewCases
        num_total_pop <- num_total_pop+part_result$Total
        
        if(ii==1)
        {
          ListInciMen[[jj]] <- part_result;
          names(ListInciMen)[jj] = names(result[[ii]])[[jj]]
        } else
        {
          ListInciWoMen[[jj]] <- part_result;
          names(ListInciWoMen)[jj] = names(result[[ii]])[[jj]]
        }
        
        if(jj==n_subpop)
        {
          part_result= data.frame(year=result$firstyear-1+1:length(num_susceptibles_all));
          part_result$Incidence <- num_newinfections_all/num_susceptibles_all
          part_result$NewCases <- num_newinfections_all
          part_result$Susceptibles <- num_susceptibles_all
          part_result$Total <- num_total_pop;
          
          ListInciAll[[ii]] <- part_result
        }
      }#End n_subpop
      
    }#End for ii
    names(ListInciAll) =c("Men", "Women")
    part_result= data.frame(year=result$firstyear-1+1:length(result[[1]][[1]][[1]]));
    total <- (ListInciAll$Men$Total+ListInciAll$Women$Total);
    part_result$Incidence <- (ListInciAll$Men$NewCases+ListInciAll$Women$NewCases)/(ListInciAll$Men$Susceptibles+ListInciAll$Women$Susceptibles)
    part_result$NewCases <- (ListInciAll$Men$NewCases+ListInciAll$Women$NewCases)
    part_result$Susceptibles <- (ListInciAll$Men$Susceptibles+ListInciAll$Women$Susceptibles)
    part_result$Total <- total;
    ListInciAll[[3]] <- part_result;
    names(ListInciAll)[3]="Women + Men ";
    
    class(ListInciAll)="SummarySyphilisPopulationIncidence";
    class(ListInciMen) = "SummaryMenIncidence"
    class(ListInciWoMen) = "SummaryWoMenIncidence"
    
    result <- list(summary=ListInciAll,men=ListInciMen,women=ListInciWoMen)
    class(result) = "SyphilisIncidenceSummary"
    return(result)
  }
}

#' @description Display results in chart
#' @param list of class "SyphilisIncidenceSummary"
myplotSyphInci <- function(X, fmain="Syphilis Incidence")
{
  maxinci = max(X$Incidence, na.rm=T);
  pow=1
  if(maxinci>0)
  {
    while(maxinci<0.1)
    {
      pow=pow*10;
      maxinci=maxinci*pow;
    }
  }
  
  fylim = c(min(pow*X$Incidence,na.rm=T),max(pow*X$Incidence,na.rm=T))
  plot(x=X$year,y=pow*X$Incidence,xlab="Time (years)",ylab= paste("Incidence x",pow, ", py", sep=""),type="l",cex=1.2,col="blue",ylim=fylim, main=fmain)
}

#' @description Display Incidence results in chart
#' @param result: list of class "SyphilisIncidenceSummary"
#' @param pop: string, population type. Should be either "all" (default) or "Men" or "Women"
plot.SummarySyphilisPopulationIncidence <- function(result,pop="all")
{
  if(class(result)=="SummarySyphilisPopulationIncidence")
  {
    if(pop=="all")
    {
      GX = result$'Women + Men '
      myplotSyphInci(GX,fmain="Syphilis, Women + Men ")
    } else if(pop=="Men")
    {
      GX = result$Men
      myplotSyphInci(GX,fmain="Syphilis among men")
    }else if(pop=="Women")
    {
      GX = result$Women
      myplotSyphInci(GX,fmain="Syphilis among women")
    } else
    {
      windows()
      par(mfrow=c(3,1))
      GX = result$'Women + Men '
      myplotSyphInci(GX,fmain="Syphilis, Women + Men ")
      GX = result$Men
      myplotSyphInci(GX,fmain="Syphilis among men")
      GX = result$Women
      myplotSyphInci(GX,fmain="Syphilis among women")
    }
  } else
  {
    stop("class of result should be SummarySyphilisPopulationIncidence")
  }
}

#' @description Display prevalence results in chart
#' @param result: list of class "SummaryMenIncidence"
#' @param pop: string, population type. Should be either "all" or men sub-populations
plot.SummarySyphilisMenIncidence<- function(result,pop="all")
{
  if(class(result)=="SummaryMenIncidence")
  {
    newpop= which(names(result)==pop)
    if(length(newpop)==1)
    {
      fname_out=paste("Syphilis among men,", pop,sep=" ")
      GX = result[[newpop]]
      myplotSyphInci(GX,fmain=fname_out)
    } else
    {
      windows()
      par(mfrow=c(2,2))
      for(ii in 1:(length(result)-1))
      {
        
        rg <- names(result)[[ii]];
        if(rg=="LowRisk") 
        {
          rg <- "Low Risk"
        } else if (rg=="MediumRisk") 
        {
          rg <- "Medium Risk"
        } else if (rg=="HighRisk")
        {
          rg <- "High Risk"
        }
        
        fname_out=paste("Syphilis among men,", rg,sep=" ")
        GX = result[[ii]]
        myplotSyphInci(GX,fmain=fname_out)
      }
    }
  } else
  {
    stop("class of result should be SummaryMenIncidence")
  }
}

#' @description Display prevalence results in chart
#' @param result: list of class "SummaryWoMenIncidence"
#' @param pop: string, population type. Should be either "all" or women sub-populations
plot.SummarySyphilisWomenIncidence<- function(result,pop="all")
{
  if(class(result)=="SummaryWoMenIncidence")
  {
    newpop= which(names(result)==pop)
    if(length(newpop)==1)
    {
      fname_out=paste("Syphilis among women,", pop,sep=" ")
      GX = result[[newpop]]
      myplotSyphInci(GX,fmain=fname_out)
    } else
    {
      windows()
      par(mfrow=c(1,3))
      for(ii in 1:(length(result)-1))
      {
        
        rg <- names(result)[[ii]];
        if(rg=="LowRisk") 
        {
          rg <- "Low Risk"
        } else if (rg=="MediumRisk") 
        {
          rg <- "Medium Risk"
        } else if (rg=="HighRisk")
        {
          rg <- "High Risk"
        }
        
        fname_out=paste("Syphilis among women,", rg,sep=" ")
        GX = result[[ii]]
        myplotSyphInci(GX,fmain=fname_out)
      }
    }
  } else
  {
    stop("class of result should be SummaryWomenIncidence")
  }
}


#' @description Save Syphilis simulation results to xlsx file
#' @param x: list of class "SyphilisPopulation"
#' @param infilename: output file name
#' @param allRes: boolean, specifying if all results should be printed in the file
#' @return NULL
WriteSyphResToDb <- function(x,infilename=NULL,allRes=FALSE)
{
  #print("aa") 
  infilename <- infilename
  if(is.null(infilename))
  {
    infilename <- x$filename
  }
  
  if(class(x)=="SyphilisPopulation")
  {
   # {# wb2 <- openxlsx::loadWorkbook("data-raw/Template_ResultFile_v3.xlsx"); default_outFile <- wb2; save(file="data/ElineWb.RData","default_outFile")
      # wb2 <- openxlsx::loadWorkbook("data/OUTPUT_v1.xlsx"); default_outFile <- wb2; save(file="data/ElineWb.RData","default_outFile")
      # wb <- openxlsx::copyWorkbook(default_outFile);   openxlsx::saveWorkbook(wb, file="qq.xlsx")
    
    wb <- openxlsx::copyWorkbook(default_outFile)
    n_allsheets <- openxlsx::sheets(wb);
    if(is.element("ProjectionResults",n_allsheets))
    {
      data_na <- data.frame(matrix(NA, nrow=206, ncol=1000))
      openxlsx::writeData(wb,"ProjectionResults",data_na, startCol = 3, startRow = 2, rowNames = FALSE, colNames = FALSE);

    } else
    {
      openxlsx::addWorksheet(wb,"ProjectionResults");
    }


    if(is.element("ProjSyph",n_allsheets))
    {
      openxlsx::removeWorksheet(wb,"ProjSyph")
    }
    if(allRes)  openxlsx::addWorksheet(wb,"ProjSyph") ;

    currind=1;

    AllPrevalence <- CalcPrev(x);
    AllIncidence <- CalcIncidence(x);

    subt <- names(AllPrevalence);
    subt[1] <- "Overall";
    subt[2] <- "Men";
    subt[3] <- "Women"

    hs1 <- openxlsx::createStyle(fgFill = "#DCE6F1", halign = "CENTER", textDecoration = "italic", fontColour = "indianred4", bgFill="#F0E68C",  border = "Bottom")

    hs2 <- openxlsx::createStyle(fontColour = "#556B2F", fgFill = "#4F80BD",
                       halign = "center", valign = "center", textDecoration = "bold",  border = "TopBottomLeftRight")

    hs3 <- openxlsx::createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center",
                               fgFill = "#4F81BD", border="TopBottom", borderColour = "#4F81BD")

    GrossiereteDeline <- NULL;
    nomsGrossiereteDeline <- NULL;
    VectYears <- NULL;

    for(ll in 1:3)
    {
      if(allRes)
      {
        openxlsx::writeData(wb,"ProjSyph",data.frame(x=subt[ll]),startRow =currind,colNames = FALSE, rowNames = FALSE, borders="rows", headerStyle = hs1 )
        openxlsx::addStyle(wb, sheet = "ProjSyph", hs3, rows = currind, cols = 1:(nrow(AllPrevalence[[ll]][[1]])+2), gridExpand = TRUE)
      }#End if(allRes)

      currind = currind+1
      for(ii in 1:length(AllPrevalence[[ll]]))
      {
        xPrev <- t(data.matrix(AllPrevalence[[ll]][[ii]]));
        xInci <- t(data.matrix(AllIncidence[[ll]][[ii]]));
        subsubtl <- names(AllPrevalence[[ll]])[ii]
        if((ll>1))
        {
          if(subsubtl=="LowRisk")
          {
            subsubtl= paste(subt[ll], "Low Risk", sep=", ")
          } else if(subsubtl=="MediumRisk")
          {
            subsubtl= paste(subt[ll], "Medium Risk",  sep=", ")
          } else if(subsubtl=="MSM")
          {
            subsubtl= "Men who have sex with men (MSM)"
          } else if((ll==2) & (subsubtl=="HighRisk") )
          {
            subsubtl= paste(subt[ll],"High Risk",  sep=", ")
          } else if((ll==3) & (subsubtl=="HighRisk") )
          {
            subsubtl= "Women, High Risk (FSW)"
          } else if(subsubtl=="NoSex")
          {
            subsubtl= paste(subt[ll], "Not Sexually Active" , sep=", ")
          }
        }

        if(allRes)
        {
          openxlsx::writeData(wb,"ProjSyph",data.frame(x=subsubtl),startRow =currind,colNames = FALSE, rowNames = FALSE, borders="rows", headerStyle = hs2)
          openxlsx::addStyle(wb, sheet = "ProjSyph", hs1, rows = currind, cols = 1:(nrow(AllPrevalence[[ll]][[ii]])+2), gridExpand = TRUE)
        }

        currind=currind+1
        finalResult <- rbind(xInci[c(1,3,2),],xPrev[-c(1,15,16),])
        newnames <- rownames(finalResult)
        newnames[c(2:17)] <- c("Incident Cases",
                               "Incidence Rate, per 1000",
                               "Prevalence %, Primary+Secondary",
                               "Prevalence %, Latent",
                               "Prevalence %, TPHA+",
                               "Prevalence %, RPR+, TPHA+",
                               "Number Screened with RPR" ,
                               "Number Screened Positive with RPR",
                               "Number Screened with TPHA only" , "Number Screened Positive with TPHA only",
                               "Number clinically treated", "Number of treatment success among clinically treated",
                               "Number of contacts traced","Number of contacts traced and found Prim_Sec","Number of contact traced in Latent",
                               "Total Population")

        #NEW UNDER ELINE's REQUEST
        finalResult[6,] <- finalResult[6,]+finalResult[7,]
        rownames(finalResult) <- newnames;

        if(allRes)
        {
          extranames <- data.frame(x=rep(subsubtl,nrow(finalResult)))
          openxlsx::writeData(wb,"ProjSyph",extranames, startCol = 1, startRow = currind, rowNames = FALSE, colNames = FALSE);
          openxlsx::writeData(wb,"ProjSyph",data.frame(finalResult), startCol = 2, startRow = currind, rowNames = TRUE, colNames = FALSE);
        }

        prevPrimSecAmongTraced <- finalResult[15,]/finalResult[14,]
        prevPrimSecAmongTraced[finalResult[14,]==0] <- 0

        rownames(xPrev)[c(15,16)] = c("Prevalence %, Incidentally cured (Compartment 5)", "Prevalence %, Recovered after treatment (Compartment 4)")

        part_GrossiereteDeline <- rbind(finalResult,prevPrimSecAmongTraced,xPrev[c(15,16),]);
        part_GrossiereteDeline <- part_GrossiereteDeline[c(2:5,19:20,6:12, 14:16,18,17),]
        part_GrossiereteDeline <- data.matrix(part_GrossiereteDeline);

        pivotname1 <- rep(subsubtl,nrow(part_GrossiereteDeline))
        pivotname2 <- rownames(part_GrossiereteDeline)
        pivotname2[17] <- "Prevalence %, Prim_Sec among contacts traced (by index patient's group)"
        part_nomGrossiereteDeline <- data.frame(cbind(pivotname1,pivotname2),row.names = NULL)

        GrossiereteDeline <- rbind(GrossiereteDeline,part_GrossiereteDeline)

        nomsGrossiereteDeline <- rbind(nomsGrossiereteDeline,part_nomGrossiereteDeline)
        if(is.null(VectYears)) VectYears <- finalResult[1,]

        currind = currind+length(newnames)+1;
      }
      currind = currind+1;
    }#End for(ll in 1:3)

    goodind <- which(nomsGrossiereteDeline[,2]!="Number of contact traced in Latent")

    indMen <- which(nomsGrossiereteDeline[goodind,]$pivotname1=="Men");
    indWomen <- which(nomsGrossiereteDeline[goodind,]$pivotname1=="Women");
    indBoth <- which(nomsGrossiereteDeline[goodind,]$pivotname1==unique(nomsGrossiereteDeline[goodind,]$pivotname1)[3]);
    indMRG <- which(is.element(nomsGrossiereteDeline[goodind,]$pivotname1,
                               c(unique(nomsGrossiereteDeline[goodind,]$pivotname1)[4:8])));
    indWRG <- which(is.element(nomsGrossiereteDeline[goodind,]$pivotname1,
                               c(unique(nomsGrossiereteDeline[goodind,]$pivotname1)[9:12])));
    reOrderEl <- c(indWomen,indMen,indBoth,indWRG,indMRG)

    years <- data.frame(t(VectYears))
    rownames(years) <- "years"
    openxlsx::writeData(wb,"ProjectionResults",years, startCol = 2, startRow = 1, rowNames = TRUE, colNames = FALSE);
    
    ind_inc <- grep("Incidence Rate, per 1000",nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2);
    GrossiereteDeline[goodind,][reOrderEl,][ind_inc,] <- 1000*GrossiereteDeline[goodind,][reOrderEl,][ind_inc,];
    
    ind_prev <- which(grepl("Prevalence %,",nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2,fixed=TRUE));
    GrossiereteDeline[goodind,][reOrderEl,][ind_prev,] <- 100*GrossiereteDeline[goodind,][reOrderEl,][ind_prev,];
    
    openxlsx::writeData(wb,"ProjectionResults",nomsGrossiereteDeline[goodind,][reOrderEl,], startCol = 1, startRow = 2, rowNames = FALSE, colNames = FALSE);
    openxlsx::writeData(wb,"ProjectionResults",data.frame(GrossiereteDeline[goodind,][reOrderEl,]), startCol = 3, startRow = 2, rowNames = FALSE, colNames = FALSE);

    if(is.element("MixingMatrices",n_allsheets))
    {
      openxlsx::removeWorksheet(wb,"MixingMatrices");
    }
    openxlsx::addWorksheet(wb,"MixingMatrices");

    copy_of_first_year = x$firstyear;
    copy_of_last_year = VectYears[length(VectYears)];

    head1 <- data.frame(x=paste("Mixing Matrix, end of ", copy_of_first_year));
    openxlsx::writeData(wb,"MixingMatrices",head1, startCol = 1, startRow = 1, rowNames = FALSE, colNames = FALSE);

    temp_res <- data.frame(x$InitMixingMat$CombinedMat);
    temp_res$Sum = apply(temp_res[,],1,sum)
    openxlsx::writeData(wb,"MixingMatrices", temp_res[-c(5,9),-c(5,9)], startCol = 2, startRow = 2, rowNames = TRUE, colNames = TRUE);

    head2 <- data.frame(x=paste("Mixing Matrix, end of ", copy_of_last_year));
    openxlsx::writeData(wb,"MixingMatrices",head2, startCol = 1, startRow = 2+1+ncol(x$InitMixingMat$CombinedMat), rowNames = FALSE, colNames = FALSE);

    temp_res <- data.frame(x$FinalMixingMat$CombinedMat);
    temp_res$Sum = apply(temp_res[,],1,sum);
    openxlsx::writeData(wb,"MixingMatrices", temp_res[-c(5,9),-c(5,9)], startCol = 2, startRow = 2+2+ncol(x$InitMixingMat$CombinedMat), rowNames = TRUE, colNames = TRUE);

    head3 <- data.frame(x=paste("Balanced Yearly Number of Partners, end of ", copy_of_first_year));
    openxlsx::writeData(wb,"MixingMatrices",head3, startCol = 1, startRow = 2+2+1+2*ncol(x$InitMixingMat$CombinedMat), rowNames = FALSE, colNames = FALSE);

    temp_res <- data.frame(x$InitBalancedPartners$BalancedPart);
    temp_res$Sum = apply(temp_res[,],1,sum);
    openxlsx::writeData(wb,"MixingMatrices",temp_res[-c(5,9),-c(5,9)], startCol = 2, startRow = 2+2+2+2*ncol(x$InitMixingMat$CombinedMat), rowNames = TRUE, colNames = TRUE);

    head4 <- data.frame(x=paste("Balanced Yearly Number of Partners, end of ", copy_of_last_year));
    openxlsx::writeData(wb,"MixingMatrices",head4, startCol = 1, startRow = 1+2+2+2+3*ncol(x$InitMixingMat$CombinedMat), rowNames = FALSE, colNames = FALSE);

    temp_res <- data.frame(x$FinalBalancedPartners$BalancedPart);
    temp_res$Sum = apply(temp_res[,],1,sum);
    openxlsx::writeData(wb,"MixingMatrices",temp_res[-c(5,9),-c(5,9)], startCol = 2, startRow = 2+2+2+2+3*ncol(x$InitMixingMat$CombinedMat), rowNames = TRUE, colNames = TRUE);

    ##Extra Graphs
    if(!is.element("ExtraCharts",n_allsheets))
    {
      openxlsx::addWorksheet(wb,"ExtraCharts");
    }
    
    ind_wlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    ind_wmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    ind_whr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, High Risk (FSW)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    
    #Men
    ind_mlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    ind_mmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    ind_mhr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, High Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    ind_mmsm_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men who have sex with men (MSM)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incident Cases")
    
    m_dat_plot <- data.frame(Years=as.numeric(unlist(years[1,])))
    m_dat_plot$'Women, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wlr_nc,]
    m_dat_plot$'Women, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wmr_nc,]
    m_dat_plot$'Women, High Risk (FSW)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_whr_nc,]
    m_dat_plot$'Men, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mlr_nc,]
    m_dat_plot$'Men, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmr_nc,]
    m_dat_plot$'Men, High Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mhr_nc,]
    m_dat_plot$'Men who have sex with men (MSM)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmsm_nc,]
    
    m_names=c("Women, Low Risk","Women, Medium Risk","Women, High Risk (FSW)",
              "Men, Low Risk", "Men, Medium Risk","Men, High Risk",
              "Men who have sex with men (MSM)")
    m_colors=c("pink","brown","orange3","royalblue4","blue","navyblue",
               "saddlebrown")
    
    xmin = min(m_dat_plot$Years,na.rm=T);
    xmax = max(m_dat_plot$Years,na.rm=T)+5;
    
    ymin = 0;
    ymax = max(m_dat_plot[,which(names(m_dat_plot)!="Years")],na.rm=T);
    
    tiff("Incidence0.tiff", w=3.42*2, h=2.44*2, units="in", pointsize=8, compression="lzw", res=300)
    par(oma=c(0, 0, 0, 15.5))
    plot.new();
    plot.window(xlim=c(xmin,xmax), ylim=c(ymin,ymax), ann=FALSE, xaxs='i', yaxs='i')
    
    abline(v=axTicks(1)[axTicks(1)<=max(m_dat_plot$Years,na.rm=T)],col="grey",lty=2)
    abline(h=axTicks(2),col="grey",lty=2)
    
    with(m_dat_plot, lines(Years,get("Women, Low Risk"),col="pink", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Women, Medium Risk"),col="brown", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Women, High Risk (FSW)"),col="orange3", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men, Low Risk"),col="royalblue4", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men, Medium Risk"),col="blue", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men, High Risk"),col="navyblue", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men who have sex with men (MSM)"),col="saddlebrown", cex=1.25, lwd=2))
    
    m_xabs <- seq(xmin,xmax, by=10)
    m_yabs <- axTicks(2)
    
    m_yabs[length(m_yabs)+1] <- m_yabs[length(m_yabs)]+(m_yabs[length(m_yabs)]-m_yabs[length(m_yabs)-1])
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=+0.0, at=m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], pos=0)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.0, at=m_yabs, labels=sprintf("%s", m_yabs),pos=xmin)
    
    legend(par('usr')[2], par('usr')[4], legend = m_names,
           bty = "n", xpd=NA,
           col = m_colors,
           lty = 1, lwd=2 #c(NA,NA,2,2),
    )
    
    title(xlab="Year", line=1.75)
    title(ylab="Incident new Syphilis infections", line=1.75)
    dev.off()
    openxlsx::insertImage(wb, "ExtraCharts","Incidence0.tiff",width = 3.42*2,height = 2.44*2, startRow = 1,startCol = 1,units = "in",dpi = 300)
    
    ##########################################################################################################
    ########################################
    ind_wlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    ind_wmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    ind_whr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, High Risk (FSW)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    
    #Men
    ind_mlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    ind_mmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    ind_mhr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, High Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    ind_mmsm_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men who have sex with men (MSM)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Incidence Rate, per 1000")
    
    m_dat_plot_2 <- data.frame(Years=as.numeric(unlist(years[1,])))
    m_dat_plot_2$'Women, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wlr_nc,]
    m_dat_plot_2$'Women, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wmr_nc,]
    m_dat_plot_2$'Women, High Risk (FSW)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_whr_nc,]
    m_dat_plot_2$'Men, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mlr_nc,]
    m_dat_plot_2$'Men, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmr_nc,]
    m_dat_plot_2$'Men, High Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mhr_nc,]
    m_dat_plot_2$'Men who have sex with men (MSM)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmsm_nc,]
    
    ymin = 0;
    ymax = max(m_dat_plot_2[,which(names(m_dat_plot_2)!="Years")],na.rm=T);
    
    tiff("Incidence1.tiff", w=3.42*2, h=2.44*2, units="in", pointsize=8, compression="lzw", res=300)
    par(oma=c(0, 0, 0, 15.5))
    plot.new();
    plot.window(xlim=c(xmin,xmax), ylim=c(ymin,ymax), ann=FALSE, xaxs='i', yaxs='i')
    
    abline(v=axTicks(1)[axTicks(1)<=max(m_dat_plot_2$Years,na.rm=T)],col="grey",lty=2)
    abline(h=axTicks(2),col="grey",lty=2)
    
    with(m_dat_plot_2, lines(Years,get("Women, Low Risk"),col="pink", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Women, Medium Risk"),col="brown", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Women, High Risk (FSW)"),col="orange3", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men, Low Risk"),col="royalblue4", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men, Medium Risk"),col="blue", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men, High Risk"),col="navyblue", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men who have sex with men (MSM)"),col="saddlebrown", cex=1.25, lwd=2))
    
    m_xabs <- seq(xmin,xmax, by=10)
    m_yabs <- axTicks(2)
    
    m_yabs[length(m_yabs)+1] <- m_yabs[length(m_yabs)]+(m_yabs[length(m_yabs)]-m_yabs[length(m_yabs)-1])
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=+0.0, at=m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], pos=0)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.0, at=m_yabs, labels=sprintf("%s", m_yabs),pos=xmin)
    
    legend(par('usr')[2], par('usr')[4], legend = m_names,
           bty = "n", xpd=NA,
           col = m_colors,
           lty = 1, lwd=2 #c(NA,NA,2,2),
    )
    
    title(xlab="Year", line=1.75)
    title(ylab=paste("Incident rate of new syphilis,\n", "infections, per 1000 person-years",sep=""), line=1.75)
    dev.off()
    openxlsx::insertImage(wb, "ExtraCharts","Incidence1.tiff",width = 3.42*2,height = 2.44*2, startRow = 1,startCol = 12,units = "in",dpi = 300)
    
    
    ###############
    ind_wlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    ind_wmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    ind_whr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, High Risk (FSW)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    
    #Men
    ind_mlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    ind_mmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    ind_mhr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, High Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    ind_mmsm_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men who have sex with men (MSM)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, Primary+Secondary")
    
    m_dat_plot <- data.frame(Years=as.numeric(unlist(years[1,])))
    m_dat_plot$'Women, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wlr_nc,]
    m_dat_plot$'Women, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wmr_nc,]
    m_dat_plot$'Women, High Risk (FSW)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_whr_nc,]
    m_dat_plot$'Men, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mlr_nc,]
    m_dat_plot$'Men, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmr_nc,]
    m_dat_plot$'Men, High Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mhr_nc,]
    m_dat_plot$'Men who have sex with men (MSM)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmsm_nc,]
    
    m_names=c("Women, Low Risk","Women, Medium Risk","Women, High Risk (FSW)",
              "Men, Low Risk", "Men, Medium Risk","Men, High Risk",
              "Men who have sex with men (MSM)")
    m_colors=c("pink","brown","orange3","royalblue4","blue","navyblue",
               "saddlebrown")
    
    xmin = min(m_dat_plot$Years,na.rm=T);
    xmax = max(m_dat_plot$Years,na.rm=T)+5;
    
    ymin = 0;
    ymax = max(m_dat_plot[,which(names(m_dat_plot)!="Years")],na.rm=T);
    
    tiff("prevalence0.tiff", w=3.42*2, h=2.44*2, units="in", pointsize=8, compression="lzw", res=300)
    par(oma=c(0, 0, 0, 15.5))
    plot.new();
    plot.window(xlim=c(xmin,xmax), ylim=c(ymin,ymax), ann=FALSE, xaxs='i', yaxs='i')
    
    abline(v=axTicks(1)[axTicks(1)<=max(m_dat_plot$Years,na.rm=T)],col="grey",lty=2)
    abline(h=axTicks(2),col="grey",lty=2)
    
    with(m_dat_plot, lines(Years,get("Women, Low Risk"),col="pink", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Women, Medium Risk"),col="brown", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Women, High Risk (FSW)"),col="orange3", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men, Low Risk"),col="royalblue4", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men, Medium Risk"),col="blue", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men, High Risk"),col="navyblue", cex=1.25, lwd=2))
    with(m_dat_plot, lines(Years,get("Men who have sex with men (MSM)"),col="saddlebrown", cex=1.25, lwd=2))
    
    m_xabs <- seq(xmin,xmax, by=10)
    m_yabs <- axTicks(2)
    
    m_yabs[length(m_yabs)+1] <- m_yabs[length(m_yabs)]+(m_yabs[length(m_yabs)]-m_yabs[length(m_yabs)-1])
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=+0.0, at=m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], pos=0)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.0, at=m_yabs, labels=sprintf("%s", m_yabs),pos=xmin)
    
    legend(par('usr')[2], par('usr')[4], legend = m_names,
           bty = "n", xpd=NA,
           col = m_colors,
           lty = 1, lwd=2 #c(NA,NA,2,2),
    )
    
    title(xlab="Year", line=1.75)
    title(ylab="Prevalence %, Primary+Secondary", line=1.75)
    dev.off()
    openxlsx::insertImage(wb, "ExtraCharts","prevalence0.tiff",width = 3.42*2,height = 2.44*2, startRow = 25,startCol = 1,units = "in",dpi = 300)
    
    ##########################################################################################################
    ########################################
    ind_wlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    ind_wmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    ind_whr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Women, High Risk (FSW)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    
    #Men
    ind_mlr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Low Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    ind_mmr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, Medium Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    ind_mhr_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men, High Risk" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    ind_mmsm_nc <- which(nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname1=="Men who have sex with men (MSM)" & nomsGrossiereteDeline[goodind,][reOrderEl,]$pivotname2=="Prevalence %, RPR+, TPHA+")
    
    m_dat_plot_2 <- data.frame(Years=as.numeric(unlist(years[1,])))
    m_dat_plot_2$'Women, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wlr_nc,]
    m_dat_plot_2$'Women, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_wmr_nc,]
    m_dat_plot_2$'Women, High Risk (FSW)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_whr_nc,]
    m_dat_plot_2$'Men, Low Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mlr_nc,]
    m_dat_plot_2$'Men, Medium Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmr_nc,]
    m_dat_plot_2$'Men, High Risk' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mhr_nc,]
    m_dat_plot_2$'Men who have sex with men (MSM)' <- GrossiereteDeline[goodind,][reOrderEl,][ind_mmsm_nc,]
    
    ymin = 0;
    ymax = max(m_dat_plot_2[,which(names(m_dat_plot_2)!="Years")],na.rm=T);
    
    tiff("prevalence1.tiff", w=3.42*2, h=2.44*2, units="in", pointsize=8, compression="lzw", res=300)
    par(oma=c(0, 0, 0, 15.5))
    plot.new();
    plot.window(xlim=c(xmin,xmax), ylim=c(ymin,ymax), ann=FALSE, xaxs='i', yaxs='i')
    
    abline(v=axTicks(1)[axTicks(1)<=max(m_dat_plot_2$Years,na.rm=T)],col="grey",lty=2)
    abline(h=axTicks(2),col="grey",lty=2)
    
    with(m_dat_plot_2, lines(Years,get("Women, Low Risk"),col="pink", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Women, Medium Risk"),col="brown", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Women, High Risk (FSW)"),col="orange3", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men, Low Risk"),col="royalblue4", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men, Medium Risk"),col="blue", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men, High Risk"),col="navyblue", cex=1.25, lwd=2))
    with(m_dat_plot_2, lines(Years,get("Men who have sex with men (MSM)"),col="saddlebrown", cex=1.25, lwd=2))
    
    m_xabs <- seq(xmin,xmax, by=10)
    m_yabs <- axTicks(2)
    
    m_yabs[length(m_yabs)+1] <- m_yabs[length(m_yabs)]+(m_yabs[length(m_yabs)]-m_yabs[length(m_yabs)-1])
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=+0.0, at=m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], m_xabs[m_xabs<=max(m_dat_plot$Years,na.rm=T)], pos=0)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.0, at=m_yabs, labels=sprintf("%s", m_yabs),pos=xmin)
    
    legend(par('usr')[2], par('usr')[4], legend = m_names,
           bty = "n", xpd=NA,
           col = m_colors,
           lty = 1, lwd=2 #c(NA,NA,2,2),
    )
    
    title(xlab="Year", line=1.75)
    title(ylab=paste("Prevalence %, RPR+, TPHA+",sep=""), line=1.75)
    dev.off()
    openxlsx::insertImage(wb, "ExtraCharts","prevalence1.tiff",width = 3.42*2,height = 2.44*2, startRow = 25,startCol = 12,units = "in",dpi = 300)
    ##End Extra graphs
    
    openxlsx::removeWorksheet(wb, "Trend graphs");
    openxlsx::removeWorksheet(wb, "Summary+Cost-Effect");
    openxlsx::removeWorksheet(wb, "Readme");
    openxlsx::removeWorksheet(wb, "Dictionary");
    openxlsx::removeWorksheet(wb, "Lists dropdown menus");
    
    
    #openxlsx::dataValidation(wb, "Trend graphs", col = 12, rows = 2, type = "list", value = "'ProjectionResults'!$C$1:$CE$1")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 12, rows = 3, type = "list", value = "'ProjectionResults'!$C$1:$CE$1")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 12, rows = 4, type = "list", value = "'Prevalence'!$N$17:$N$33")

    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 7, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 8, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 9, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 10, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 11, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 12, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    #openxlsx::dataValidation(wb, "Trend graphs", col = 11, rows = 13, type = "list", value = "'Trend graphs'!$K$21:$K$31")
    
    #openxlsx::dataValidation(wb, "Summary+Cost-Effect", col = 3, rows = 1, type = "list", value = "'Lists dropdown menus'!$E$2:$E$82")
    #openxlsx::dataValidation(wb, "Summary+Cost-Effect", col = 11, rows = 2, type = "list", value = "'Lists dropdown menus'!$A$2:$A$3")
    #openxlsx::dataValidation(wb, "Summary+Cost-Effect", col = 12, rows = 2, type = "list", value = "'Lists dropdown menus'!$B$2:$B$3")
    
    openxlsx::saveWorkbook(wb,file=infilename,overwrite=TRUE)
    unlink("Incidence0.tiff")
    unlink("Incidence1.tiff")
    if(file.exists("Incidence0.tiff"))
    {
      file.remove("Incidence0.tiff")
    }
    if(file.exists("Incidence1.tiff"))
    {
      file.remove("Incidence1.tiff")
    }
    
    unlink("prevalence0.tiff")
    unlink("prevalence1.tiff")
    if(file.exists("Incidence0.tiff"))
    {
      file.remove("prevalence0.tiff")
    }
    if(file.exists("prevalence1.tiff"))
    {
      file.remove("prevalence1.tiff")
    }
    TRUE
  } else
  {
    stop("class(x) must be SyphilisPopulation")
    NULL
  }
}

closealldevices <- function()
{
  while (!is.null(dev.list()))  dev.off()
}


