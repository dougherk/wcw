####################################################################################################################
#				THE PROBABILITY OF SELECTING A WEAK CONDORCET WINNER (wCW)
# 				      Written by Keith Dougherty <dougherk@uga.edu>
#                                                   Used in
#                                          Diss, Dougherty, Heckelman
#                   "When Ties are Possible: Weak Condorcet Winners and Arrovian Rationality"
#
# DESCRIPTION
#   This program calculates the probability that various preference aggregation procedures (in this case voting 
#   rules) will select a weak Condorcet Winner (wCW) among four types of transitive social rankings.  
# ASSUME
# 1- A=3.  This limits the wCW types in the study.
# 2- Preferences drawn from IAC, IC, or a weighted distribution of the strict linear preference order.
#
# NOTES
# 1- Use -9 for missing in this program.
# 2- Rcpp functions are listed in iterate2.cpp.
# 3- The program works for PAR's that produce socially strict preferences or indifference over pairs of alternatives
#    (that that produce weak preference, i.e. greater than or equal to as a category are not allowed)
# 4- currently the program counts a PAR as succesffully chosing a wCW if it 1) picks any of the weak Condorcet 
#    winners and 2) it does not pick a non-wCW.  A PAR need not pick ALL weak Condorcet winners.
###################################################################################################################

#----------------  -----------------  --------------#
# dt==1 means impartial anonymous culture condition #
#----------------  -----------------  --------------#

rm(list=ls())
library(devtools)
library(Rcpp)
library(gtools)
setwd('C:/Users/dougherk/Dropbox/weak_condorcet/programs/')
sourceCpp("c_functions/iterate2.cpp")

# Declare some primitives (if you want to do one case at a time; otherwise use input file)
set.seed(021114)
N<- 4				# the number of individuals
A<- 3				# the number of alternatives.			
dt<-1				# preference distribution type: =0 if IC, =1 if IAC, =2 if weighted.
T<-10				# the number of trials in the simulation.
V<-c(0,0,0,0,0,0)		# preference rates for dt==2 (use 0s for dt==0 and dt==1 so we don't get confused)

# Read input file
  #input<-read.delim('test_in.txt', header = TRUE, sep = "\t")
  input<-read.delim('input.txt', header = TRUE, sep = "\t")

for(tt in 1:nrow(input)){	# for each set of inputs in the input files (tt is a row in the matrix of inputs)
    N<-input[tt,"N"]
    A<-input[tt,"A"]
    dt<-input[tt,"dt"]
    T<-input[tt,"T"]

  # Generate all possible permutations of an n-tuple using gtools (i.e., all possible strict linear preference orders)
    perm <-permutations(n = A, r = A, v = 0:(A-1))	

  # RUN THE SIMULATIONS FOR T TRIALS
    out<-iterate(perm, N, T, dt, V)			
  
  # CALCULATE THE PROPORTION OF VIOLATIONS IN THE T TRIALS
    type <- out$type
    out2 <- rbind(out$count_prr,out$count_aprr,out$count_borda,out$count_nanson,out$count_hare,out$count_pmr,out$count_copeland,out$count_coombs,out$count_dowdall,out$count_sk,out$count_baldwin,out$count_rp)
    type_tot <- sum(type[1:4])
    out2_tot<-rowSums(out2)
    out3<-t( t(out2) / type[1:4])

  # PRINT VOTING RULE OUTPUT
    file1<-sprintf("output/out_totals_N%s_A%s_dt%s_T%s.csv", as.character(N), as.character(A), as.character(dt), as.character(T))
    header1<-"The probability of selecting a Weak Condorcet Winner from wcw.R and iterate2.cpp."
    header2<-sprintf("     Assuming: N=%s; A=%s; dt=%s; Trials=%s.", as.character(N), as.character(A), as.character(dt), as.character(T))
    header3<-sprintf("     Pref Weights=(%s,%s,%s,%s,%s,%s).\n", as.character(V[1]), as.character(V[2]), as.character(V[3]), as.character(V[4]),  as.character(V[5]), as.character(V[6]))
    header4<-sprintf("The total number of intransitive draws are:")
    header5<-sprintf("%s.\n", as.character(out$num_intransitive))
     if(dt==2){
      write.table(rbind(header1,header2,header3,header4,header5), file1, quote = FALSE, row.names = FALSE, col.names=FALSE)
    } 
    if(dt<2){
      write.table(rbind(header1,header2,header4,header5), file1, quote = FALSE, row.names = FALSE, col.names=FALSE)
    } 
    # append the numeric output to the output file
      write.table("The frequency that each PAR selects the corresponding wCW types:", file=file1, append = TRUE, quote = FALSE, row.names = FALSE, col.names=FALSE)
      colnames(out3)<-c("type 1","type 2","type 3","type 4")
      rownames(out3)<-c("Plurality Rank","Anti-PRR","Borda Count","Nanson Rule","Hare Rule","Pairwise Maj","Copeland","Coombs", "Dowdall","Simpson-Kramer","Baldwin","Ranked Pairs")  
      out3 <- cbind(PAR = rownames(out3), out3)
      rownames(out3)<- NULL
      write.table(out3, file=file1, append = TRUE, sep=",", quote = FALSE, row.names = FALSE, col.names=TRUE)
      write.table("\nThe total number of wCW types 1-4 drawn out of 1 million are:", file=file1, append = TRUE, quote = FALSE, row.names = FALSE, col.names=FALSE)
      type<- c("",type)
      names(type)<-c("","type 1","type 2","type 3","type 4","type 5","type 6","type 7")
      write.table(t(type), file=file1, append = TRUE, sep=",", quote = FALSE, row.names = FALSE, col.names=TRUE)
  
  # PRINT TYPE OUTPUT IN SINGLE TABLE
    file2<-sprintf("output/out_types_various_N_A%s_dt%s_T%s.csv", as.character(A), as.character(dt), as.character(T))
    if(tt==1){  
      header6<-sprintf("     Assuming: A=%s; dt=%s; Trials=%s.", as.character(A), as.character(dt), as.character(T))
      header7<-"\nThe total number of wCW types 1-7 drawn out of T trials are:"
      write.table(rbind(header1,header6,header7), file2, quote = FALSE, row.names = FALSE, col.names=FALSE)
      header8 <- c("N","type 1","type 2","type 3","type 4","type 5","type 6","type 7")
      write.table(t(header8), file=file2, append = TRUE, sep=",", quote = FALSE, row.names = FALSE, col.names=FALSE)
    }    	
    # append the numeric output to the output file
      ntype<- c(N, type[2:length(type)])
      names(ntype)<-c("N","type 1","type 2","type 3","type 4","type 5","type 6","type 7")
      write.table(t(ntype), file=file2, append = TRUE, sep=",", quote = FALSE, row.names = FALSE, col.names=FALSE)

}





