/* ITERATE

	This function manages the iterations of the program as well as the order in which various functions are called.  Presumably, they will be loaded in R
	outside this function.

   NOTE:
   F	 is the frequency distribution on perm (i.e., the number of people with each of the preference orders listed left to right in perm).

   VERSION: This version of iterate (version 2) is slower than the version without a suffix.  The version w/o suffix skips the evaluation of a normative
   criteria (transitivity or IIA) if an earlier normative criteria was violated (like Pareto or transitivity for the case of IIA).  Version 2, the one used
   here, continues on if there is a contradition so we can have a more accurate picture of the probabiliyt of violating Pareto, transitivity, and IIA
   sepearately.

   The two "spread the exalted pair" algorithms, iia_move_two() and iia_move_twoB(), are included.
   Also included are the coombs, dowdall, and simpson-kramer methods.

   INPUT:
   perm	 is an Af x A matrix of individual preferences, ordered most preferred (left) to least preferred (right).
   N	 is the number of voters.
   T	 is the number of iterations in the simulation.
   dt	 is the distribution type: dt=0 for IC; dt=1 for IAC.

*/


#include <Rcpp.h>
using namespace Rcpp;

/* ---------------------------------------- PREFERENCE DISTRIBUTIONS ---------------------------------------------------- */

/* ICC: IMPARTIAL CULTURE CONDITION */

// [[Rcpp::export]]

	IntegerVector icc(IntegerMatrix P, int N) {
			int n = P.nrow();					/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
			IntegerVector F(n,0) ;				/* F is vector, length n, that records the number of individuals with the respective row preference in P. */
			NumericVector bins(n) ;				/* bins records the upper bound of each bin, which must be between 0 and 1 */
		    double ran ;
			double h = 1/ (double) n;

		    /* create n bins in the [0,1] interval */
			for(int j=0; j < n; j++){
			  bins(j) = h*(j+1);
			};

			/* randomly draw N individuals and put them in the appropriate bin for F */
			for(int i=0; i < N; i++){			/* for each individual in that vector */
		      ran = R::runif(0,1);				/* randomly draw a number between 0 and 1 inclusive */
		      if( ran <= bins(0) ){
				F(0)++;							/* increase F in the first position if you found a match in the first bin */
			  }
			  else{
		  	    for(int j=1; j < n; j++){		/* go through the bins and see which bin the number belongs to */
		    	  if( (ran <= bins(j)) & (ran > bins(j-1)) ){
					F(j)++;						/* increase F in the j position if you found a match with bin j */
				  }
			    };
		  	  }
			};

	return F;
	}

/* IAC: IMPARTIAL ANONYMOUS CULTURE CONDITION */

// [[Rcpp::export]]


	Rcpp::IntegerVector iac(IntegerMatrix P, int N)
	{
	   int n = P.nrow();					/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
	   Rcpp::NumericVector p(n-1);
  	   IntegerVector F(n,0) ;				/* F is vector, length n, that records the number of individuals with the respective row preference in P. */
	   double ran ;

    /* Generate a probability distribution, q, of each of the n orderings using the broken stick method */
       /* randomly draw a vector of numbers, p, and sort */
       p = runif(n-1);						/* randomly draw a vector of numbers between 0 and 1 */
	   std::sort(p.begin(), p.end());		/* sort the numbers in ascending order */

	/* Randomly draw N individuals and put them in the appropriate bin for F (following the broken stick rules) */
	   for(int i=0; i < N; i++){			/* for each individual */
	     ran = R::runif(0,1);				/* randomly draw a number (scalar) between 0 and 1 inclusive */
	     if( ran <= p(0) ){					/* first element match */
			F(0)++;							/* increase F in the first position if you found a match in the first bin */
		 }
		 else{
			if( ran > p(n-2) ){			/* last element match */
			  F(n-1)++;
			}
			else{
			   for(int j=1; j < n-1; j++){
			     if( (ran <= p(j)) & (ran > p(j-1)) ){		/* middle element match */
				   F(j)++;
				 }
		  	   };
			}
  	     }
	   };

	return(F);
	}

/* WEIGHTED PREFERENCE DISTRIBUTION, dt=2, (ASSUMES A=3) */

// [[Rcpp::export]]

	IntegerVector weighted_pref(IntegerMatrix P, int N, NumericVector V) {
			int n = P.nrow();					/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
			IntegerVector F(n,0) ;				/* F is vector, length n, that records the number of individuals with the respective row preference in P. */
			NumericVector bins(n) ;				/* bins records the upper bound of each bin, which must be between 0 and 1 */
		    double ran ;
			double h = 1/ (double) n;

  			if ( V.length() != n ){
			  Rcpp::stop("The length of the weight vector is not equal to A!");
			}

		    /* create n bins in the [0,1] interval: which are cummulative additions of each weight */
		    bins(0) = V(0);
			for(int j=1; j < n; j++){
			  bins(j) = V(j) + bins(j-1);
			};

			/* randomly draw N individuals and put them in the appropriate bin for F */
			for(int i=0; i < N; i++){			/* for each individual in that vector */
		      ran = R::runif(0,1);				/* randomly draw a number between 0 and 1 inclusive */
		      if( ran <= bins(0) ){
				F(0)++;							/* increase F in the first position if you found a match in the first bin */
			  }
			  else{
		  	    for(int j=1; j < n; j++){		/* go through the bins and see which bin the number belongs to */
		    	  if( (ran <= bins(j)) & (ran > bins(j-1)) ){
					F(j)++;						/* increase F in the j position if you found a match with bin j */
				  }
			    };
		  	  }
			};

	return F;
	}

/* --------------------------------------------- PREFERENCE AGGREGATION RULES ----------------------------------------------- */

/* PLURALITY RANKING RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix prr(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerVector V(A,0);					/* vector of first place votes for each alternative (a row index) */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */

	    for (int i = 0; i < n; i++){ 			/* go through each row of preference profile P */
		  V( P(i,0) ) += F(i);					/* add the total number of first place votes in V */
		}

		/* assign pairwise social preferences based on plurality votes */
	    for (int i = 0; i < A; i++){
  	      for (int j = 0; j < A; j++){
			if( V(i) > V(j) ){
			  soc(i,j)= 1;
			}
			if( V(j) > V(i) ){
			  soc(j,i)= 1;
			}
			if( V(i) == V(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  }
		}


	return soc;
	}


/* ANTI-PLURALITY RANKING RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix aprr(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerVector V(A,0);					/* vector of first place votes for each alternative (a row index) */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */

	    for (int i = 0; i < n; i++){ 			/* go through each row of preference profile P */
		  V( P(i,A-1) ) += F(i);				/* add the total number of last place votes in V */
		}

		/* assign pairwise social preferences based on last-place votes (more votes, less preferred) */
	    for (int i = 0; i < A; i++){
  	      for (int j = 0; j < A; j++){
			if( V(i) < V(j) ){
			  soc(i,j)= 1;
			}
			if( V(j) < V(i) ){
			  soc(j,i)= 1;
			}
			if( V(i) == V(j) ){
			  soc(i,j)= 0;
			}
		  }
		}

	return soc;
	}


/* BORDA COUNT */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix borda(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
     	IntegerVector BC(A, 0);					/* borda count for candidate i */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
     	int points;								/* the borda counts assigned to a particular level of preference (1st, 2nd, etc.) */

		/* assign borda count to each candidate (candidates names are row indeces in BC) */
     	for (int i = 0; i < n; i++){ 			/* for each row */
          points = A;							/* start with the maximum number of Borda points when starting each row */
          for (int j = 0; j < A; j++){			/* go through the columns and determine the borda points at each level of preference */
	    	BC( P(i,j) ) += points * F(i);
   	    	points--;
          };
   	    };


		/* assign pairwise social preferences based on Borda Count */
	    for (int i = 0; i < A; i++){
  	      for (int j = 0; j < A; j++){
			if( BC(i) > BC(j) ){
			  soc(i,j)= 1;
			}
			if( BC(j) > BC(i) ){
			  soc(j,i)= 1;
			}
			if( BC(i) == BC(j) ){
			  soc(i,j)= 0;
			}
		  }
		}

	return soc;
	}


/* MEAN WITH GATE (USED IN NANSON).  TAKES AVERAGE OF SECOND COLUMN IN AN Ax2 MATRIX */

// [[Rcpp::export]]
	double avg_w_gate(IntegerMatrix X, IntegerVector G) {
		int A = X.nrow();						/* X is the Ax2 matrix. G is a vector, length A. */
		double sums=0;

        for (int j = 0; j < A; j++){			/* go through the rows of X */
		  if(G(j) == 1){
 		    sums += X(j,1);
		  }
		}

	return (sums / sum(G));
	}


/* NANSON BORDA. RETURNS MATRIX OF ALTERNATIVES AND THEIR BORDA COUNTS (USED BY NANSON AND BALDWIN) */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix nborda(IntegerVector F, IntegerMatrix P, IntegerVector G) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
     	IntegerMatrix BC1(A, 2);				/* borda count for alternative i: column 0 = alternative name, column 1 = borda score */
		std::fill(BC1.begin(), BC1.end(), 0);
     	int points;								/* the borda counts assigned to a particular level of preference (1st, 2nd, etc.) */

     	for (int i = 0; i < A; i++){
		  BC1(i,0) = i;
		}

		/* assign borda count to each candidate (candidates names are row indeces in BC) */
     	for (int i = 0; i < n; i++){ 			/* for each row in P */
          points = A;							/* start with the maximum number of Borda points when starting each row */
          for (int j = 0; j < A; j++){			/* go through the columns and determine the borda points at each level of preference */
			for (int k = 0; k < A; k++){		/* go through the indeces of BC1 to find a matching name */
			  /* if the alternative name matches P(i,j) and the gate for P(i,j) is open (i.e., this alternative is still in the election), then ... */
			  if( (BC1(k,0) == P(i,j)) & (G( P(i,j) ) == 1)){
	 		    BC1(k,1) += points * F(i);
   	    	    points--;
			  }
			};
          };
   	    };


	return BC1;
	}


/* NANSON */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix nanson(Rcpp::IntegerVector F, Rcpp::IntegerMatrix perm) {
		IntegerMatrix P(clone(perm));
		int A = P.ncol();						/* A is the number of alternatives.  P is the permutation matrix, perm.   */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
		double m;								/* the mean of the borda count scores */
		IntegerMatrix BC(A,2);				 	/* borda count for alternative i: column 0 = alternative name, column 1 = borda score */
     	IntegerVector relim(A, A);				/* relim records the round an alt was eliminated with the vector index used as the alt name (initially A) */
    	IntegerVector G(A, 1);					/* G is a gate length A. G = 1 if alternative at that index is included in BC; 0= if eliminated  */
		int LessA = 0;

		/* run borda count and calcualte mean of scores */
		int r = 1;								/* r is the round, used to keep track of round that each alternative is eliminated */
		do {									/* repeat until only one alternative remains */
		  BC = nborda(F,P,G);
		  m = avg_w_gate(BC, G);				/* calculate the average borda count values (avg of column 1) for alternatives with open gate (col2==1) */
	      for (int i = 0; i < A; i++){			/* walk through each alternative */
		    if(BC(i,1) <= m){					/* for alternatives AT OR BELOW the mean borda score, we */

			  /* if not previously eliminated: record the round eliminated (value) in relim for each alternative (vector index), */
			  if( G(BC(i,0)) == 1 ){
			    relim( BC(i,0) ) = r;
			  }
		 	  G(BC(i,0)) = 0;					/* then mark the alternative as eliminated in the gate for the next round, */
			  LessA++;							/* and keep track of the number of alternatives eliminated so we know when to stop. */

		    }
		  };
		r++;
		} while(A - LessA > 1);					/* if more than one alternative remains, we do it again */


		/* assign pairwise social preferences based on larger r (i.e., later round eliminated) */
	    for (int i = 0; i < A; i++){
  	      for (int j = i; j < A; j++){
			if( relim(i) > relim(j) ){
			  soc(i,j)= 1;
			}
			if( relim(j) > relim(i) ){
			  soc(j,i)= 1;
			}
			if( relim(i) == relim(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  };
		};

	return soc;
	}



/* HARE PLURALITY (USED BY HARE RULE). RETURNS THE PLURALITY VOTE FOR EACH ALTERNATIVE (note: P cannot shrink) */

// [[Rcpp::export]]

	Rcpp::IntegerVector hplurality(IntegerVector F, IntegerMatrix P) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerVector V(A,0);					/* vector of first place votes for each alternative (a row index) */

	    for (int i = 0; i < n; i++){ 			/* go through each row of preference profile P */
		  V( P(i,0) ) += F(i);					/* add the total number of first place votes in V */
		}

	return V;
	}

/* MOVE_BOTTOM (USED BY HARE & COOMBS). MOVES ALTERNATIVE "a" TO THE BOTTOM OF EACH PREFERENCE ORDER */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix move_bottom(IntegerMatrix P, int a) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix P2(n, A);					/* new permutation matrix */
		int s=0;								/* s=0 indicates straight copy, s=1 indicates copy one column to the right */

		/* walk through P and copy alternatives into P2 that are not a (moving everything less preferred to "a" up one spot) */
     	for (int i = 0; i < n; i++){ 			/* for each row */
     	s=0;
          for (int j = 0; j < A-1; j++){		/* go through the columns and determine if the alternative in that (i,j) is a */
			if(P(i,j) == a){					/* as soon as we hit alternative a, we will copy next preference down the list (column to the right) */
			  s=1;
			}
			P2(i,j) = P(i,j+s);
		  };
		P2(i,A-1) = a;							/* finally put "a" in least preferred position, before going to a new row */
		};


	return P2;
	}

/* SECOND_SMALLEST (USED BY HARE RULE). IDENTIFIES THE SECOND SMALLEST INTEGER IN A VECTOR OF INTEGERS */

// [[Rcpp::export]]
	int second_smallest(IntegerVector x) {
		int n = x.length() ;
		IntegerVector y = clone(x);
		std::sort(y.begin(), y.end());
		int second = y(0);

  	    for (int i = 0; i < n; i++){
		  if(y(i) >  second){
			  second = y(i);
			  break;
	  	  }
		}

	return second;
    }


/* HARE RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix hare(Rcpp::IntegerVector F, Rcpp::IntegerMatrix perm) {
		IntegerMatrix P(clone(perm));
		int A = P.ncol();						/* A is the number of alternatives.  P is the permutation matrix, perm.   */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
	    IntegerVector V2(A);						/* records the plurality vote in a given round */
		int m;									/* the min plurality vote in a given round */
     	IntegerVector elim(A, A);				/* records the round an alt was eliminated with the vector index used as the alt name */

		/* run plurality rule and determine alternative with least number of votes */
		int r = 1;								/* r is the round, used to keep track of the round that each alternative is eliminated */
		int j = 0;								/* j keeps track of how many alternatives have been eliminated (i.e., moved to bottom) */
		do {
		  V2 = hplurality(F,P);					/* determine plurality votes for each alternative */
	      if(r==1){
			m = min(V2);						/* in round 1, determine the alternative with the smallest plurality vote */
		  }
		  else{
			m = second_smallest(V2); 			/* in later rounds, determine the alternative with second smallest vote (because those eliminate in r1 will have 0 votes) */
		  }
  	      for (int i = 0; i < A; i++){			/* walk through each alternative */
		    if(V2(i) == m){						/* if the alternative has the least number of votes, we */
		  	  elim(i) = r;						/* record the round eliminated (value), */
			  P = move_bottom(P, i);			/* then move that alternative to the bottom of the preference profile */
			  j++;
		    }
		  };
		r++;
		} while(j < A-1);							/* stop when only one alternative remains */


		/* assign pairwise social preferences based on larger r (i.e., later round eliminated) */
	    for (int i = 0; i < A; i++){
  	      for (int j = i; j < A; j++){
			if( elim(i) > elim(j) ){
			  soc(i,j)= 1;
			}
			if( elim(j) > elim(i) ){
			  soc(j,i)= 1;
			}
			if( elim(i) == elim(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  };
		};

	return soc;
	}

/* HEAD-TO_HEAD (NAME OF WINNER: USED BY PAIRWISE MAJORITY RULE, COPELAND, & RANKED PAIRS) */

// [[Rcpp::export]]

	int head_to_head(IntegerVector F, IntegerMatrix P, int x, int y) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerVector V(2,0);					/* vector of votes for x (index 0) versus y (index 1) */

	    for (int i = 0; i < n; i++){
  	      for (int j = 0; j < A; j++){
			if(P(i,j) == x){					/* if we arrive at x first (i.e. further left in P) we will assign x the votes in that row */
			  V(0) += F(i);
			  break;
			}
			if(P(i,j) == y){					/* if we arrive at y first (i.e. further left in P) we will assign y the votes in that row */
			  V(1) += F(i);
			  break;
			}
		  };
		};

		if( V(0) > V(1) ){
		  return x;
		}
		if( V(0) < V(1) ){
		  return y;
		}
		if( V(0) == V(1) ){
		  return -8;
		}
	}

/* PAIRWISE MAJORITY RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix pair_maj_rule(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
		int winner;

		/* compare each pair of alternatives pairwise head-to-head in pairwise majority rule contests */
	    for (int i = 0; i < A-1; i++){
  	      for (int j = i+1; j < A; j++){
			winner = head_to_head(F,P,i,j);
			if (winner == i){
			  soc(i,j)= 1;
			}
			else{
			  if (winner == j){
			    soc(j,i)= 1;
			  }
			  else{
			    soc(i,j)= 0;
			    soc(j,i)= 0;
		  	  }
			}
		  };
		};

	return soc;
	}

/* COLSUMS (not used) */

// [[Rcpp::export]]

  IntegerVector colSumsC(IntegerMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    IntegerVector out(ncol);

    for (int j = 0; j < ncol; j++) {
      double total = 0;
      for (int i = 0; i < nrow; i++) {
        total += x(i, j);
      }
      out[j] = total;
    }
    return out;
  }

/* ROWSUMS (for integers, used with wCW) */

// [[Rcpp::export]]

  IntegerVector rowSumsC(IntegerMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    IntegerVector out(nrow);

    for (int i = 0; i < nrow; i++) {
      double total = 0;
      for (int j = 0; j < ncol; j++) {
        total += x(i, j);
      }
      out[i] = total;
    }
    return out;
  }

/* ROWSUMS N (for numeric, used with Copeland and wCW_type) */

// [[Rcpp::export]]

  NumericVector rowSumsN(NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector out(nrow);

    for (int i = 0; i < nrow; i++) {
      double total = 0;
      for (int j = 0; j < ncol; j++) {
        total += x(i, j);
      }
      out[i] = total;
    }
    return out;
  }


/* COPELAND */

// [[Rcpp::export]]

Rcpp::IntegerMatrix copeland(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, returned */
		IntegerMatrix socp(clone(soc));			/* tournament matrix which is just 1's, 0's, and .5's */
		std::fill(soc.begin(), soc.end(), -9);	/* make off diagonal elements in soc  -9 for missing */
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
		std::fill(socp.begin(), socp.end(), 0);	/* make all elements in socp 0 */
		IntegerVector V(A);						/* total number of pairwise wins for alternative i */
		int winner;

		/* create tournament matrix by comparing each pair of alternatives head-to-head */
	    for (int i = 0; i < A-1; i++){
  	      for (int j = i+1; j < A; j++){
			winner = head_to_head(F,P,i,j);

			if (winner == i){
			  socp(i,j)= 1;
			}
			if (winner == j){
			  socp(j,i)= 1;
			}
			if (winner == -8){					/* i and j tie */
			  socp(i,j)= .5;
			  socp(j,i)= .5;
			}
		  };
		};

		/* sum the rows of the tournament matrix */
		V = rowSumsC(socp);

		/* assign pairwise social preferences based on larger V (larger row sum) */
	    for (int i = 0; i < A; i++){
  	      for (int j = i+1; j < A; j++){
			if( V(i) > V(j) ){
			  soc(i,j)= 1;
			}
			if( V(j) > V(i) ){
			  soc(j,i)= 1;
			}
			if( V(i) == V(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  };
		};

return soc;
}

/* REMOVE COLUMN (used in Coombs) */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix rmcol(Rcpp::IntegerMatrix X, int a) {
		int A = X.ncol();
		IntegerMatrix X2(X.nrow(),A-1);
		int s=0;								/* s=0 indicates straight copy, s=1 indicates copy one column to the right */

          for (int j = 0; j < A-1; j++){		/* go through the columns and determine if the column is "a" */
			if(j == a){							/* as soon as we hit a, we will copy next preference down the list (column to the right) */
			  s=1;
			}
			X2(_,j) = X(_,j+s);
		  };

	return X2;
	}


/* COOMBS ANTI-PLURALITY. RETURNS THE ANTI-PLURALITY VOTES (used in Coombs) */

// [[Rcpp::export]]

	IntegerVector cantiplural(IntegerVector F, IntegerMatrix P, int A) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int AA = P.ncol();						/* P is the permutation matrix, perm.  AA is the number of alternatives in the "potentially shrunk" matrix */
		IntegerVector V(A,0);					/* vector of last place votes for each alternative (a row index) */
		int i;

	    for (i = 0; i < n; i++){ 				/* go through each row of preference profile P */
		  V( P(i,AA-1) ) += F(i);				/* add the total number of last place votes in V.  Because 0 is the starting index, AA-1 is last place */
		};

	return V;
	}


/* MOVE_TOP (used in Coombs). MOVES ALTERNATIVE "a" TO THE TOP OF EACH PREFERENCE ORDER */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix move_top(IntegerMatrix P, int a) {
		int A = P.ncol();							/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix P2 = move_bottom(P, a);		/* taking advantage of some existing code */
		IntegerMatrix P3 = rmcol(P2, P2.ncol()); 	/* remove the last column of P2 so that P3 has one less column than P */
		IntegerVector V= rep(a,P3.nrow());			/* create a column of the alternative name that should be on top */

	return cbind(V,P3);								/* combine the alternative on top with the A-1 column P3 that has the other preferences in order */
	}

/* COOMBS RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix coombs(Rcpp::IntegerVector F, Rcpp::IntegerMatrix perm) {
		IntegerMatrix P(clone(perm));
		int A = P.ncol();						/* A is the number of alternatives.  P is the permutation matrix, perm.   */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
	    IntegerVector V2(A);					/* records the anti-plurality vote in a given round */
		int m;									/* the max anti-plurality vote in a given round */
     	IntegerVector elim(A, A);				/* records the round an alt was eliminated in with the vector index used as the alt name */
     	IntegerVector move(A, 0);				/* records the alts that have to be moved to the top each round */

		/* determine which alternative has the most least place votes */
		int r = 1;								/* r is the round, used to keep track of the round that each alternative is eliminated */
		int j = 0;								/* j keeps track of how many alternatives have been eliminated (i.e., moved to top) */
		do {
		  V2 = cantiplural(F,P,A);				/* determine the number of last place votes for each alternative */
  		  m = max(V2);							/* determine the greatest number of last place votes */

  	      for (int i = 0; i < A; i++){			/* walk through each alternative */
		    if(V2(i) == m){						/* if the alternative has the greatest number of last place votes, we */
		  	  elim(i) = r;						/* record the round eliminated (value) */
		  	  move(i) = 1;						/* record which alternatives have to be moved up this round */
			  j++;
		    }
		  };

   	      for (int i = 0; i < A; i++){			/* move all alternatives to the top of the preference profile that need to be moved this round */
			if(move(i)==1){
			  P = move_top(P, i);
			}
		  };
   	      for (int i = 0; i < A; i++){			/* reset indicator */
			move(i)=0;
		  };

		r++;
		} while(j < A);							/* stop when no alternatives remain */


		/* assign pairwise social preferences based on larger r (i.e., later round eliminated is more preferred) */
	    for (int i = 0; i < A; i++){
  	      for (int j = i; j < A; j++){
			if( elim(i) > elim(j) ){
			  soc(i,j)= 1;
			}
			if( elim(j) > elim(i) ){
			  soc(j,i)= 1;
			}
			if( elim(i) == elim(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  };
		};

	return soc;
	}

/* DOWDALL RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix dowdall(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
     	Rcpp::NumericVector DC(A, 0.0);			/* dowdall count for candidate i */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */

		/* assign dowdall count to each candidate (candidates names are row indeces in DC) */
     	for (int i = 0; i < n; i++){ 			/* for each row in P */
          for (double j = 0; j < A; j++){		/* go through the columns and determine the dowdall points at each level of preference */
	    	DC( P(i,j) ) += A/(j+1) * F(i);		/* we multiply by A to avoid dividing by 3 (infinitely repeating) for A=3; we divide by j+1 because our index starts at 0 */
          };
   	    };

		/* assign pairwise social preferences based on Dowdall Count */
	    for (int i = 0; i < A; i++){
  	      for (int j = 0; j < A; j++){
			if( DC(i) > DC(j) ){
			  soc(i,j)= 1;
			}
			if( DC(j) > DC(i) ){
			  soc(j,i)= 1;
			}
			if( DC(i) == DC(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  }
		}

	return soc;
	}

/* ROW MINIMUM (used in Simpson-Kramer) */

// [[Rcpp::export]]

	NumericVector rowMins(NumericMatrix X) {
	  int n = X.nrow();
	  NumericVector V(n);
	  for (int i=0; i<n; i++) {
	     NumericVector W = X.row(i);
	     V[i] = *std::min_element(W.begin(), W.end());  // from the STL
	  }
	  return(V);
	}

/* VOTE_MARGIN (used in Simpson-Kramer & ranked pairs) */

// [[Rcpp::export]]

	int vote_margin(IntegerVector F, IntegerMatrix P, int x, int y) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerVector V(2,0);					/* vector of votes for x (index 0) versus y (index 1) */

	    for (int i = 0; i < n; i++){
  	      for (int j = 0; j < A; j++){
			if(P(i,j) == x){					/* if we arrive at x first (i.e. further left in P) we will assign x the votes in that row */
			  V(0) += F(i);
			  break;
			}
			if(P(i,j) == y){					/* if we arrive at y first (i.e. further left in P) we will assign y the votes in that row */
			  V(1) += F(i);
			  break;
			}
		  };
		};

	return V(0) - V(1);							/* returns vote margin of x  w.r.t.  y */
	}


/* SIMPSON-KRAMER RULE */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix sk(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, returned */
		IntegerMatrix socp(clone(soc));			/* tournament matrix which includes vote margins of i over j */
		std::fill(soc.begin(), soc.end(), -9);	/* make off diagonal elements in soc  -9 for missing */
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
		std::fill(socp.begin(), socp.end(), 0);	/* make all elements in socp 0 */
		NumericMatrix VM(A ,A);					/* matrix of pairwise vote margins for i over j (can be negative) */
		NumericVector V(A);						/* total number of pairwise wins for alternative i */

		/* create tournament matrix by comparing each pair of alternatives head-to-head */
	    for (int i = 0; i < A; i++){
  	      for (int j = 0; j < A; j++){
			if (i != j){
			  VM(i,j) = vote_margin(F,P,i,j);	/* vote margin w.r.t. i (can be positive or negative) */
			}
		  };									/* close for j */
		};										/* close for i */

		/* record the minimum vote margin for each alternative i */
		V = rowMins(VM);

		/* assign pairwise social preferences based on minimum vote margin V */
	    for (int i = 0; i < A; i++){
  	      for (int j = 0; j < A; j++){
			if( V(i) > V(j) ){
			  soc(i,j)= 1;
			}
			if( V(j) > V(i) ){
			  soc(j,i)= 1;
			}
			if( V(i) == V(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  };									/* close for j */
		};										/* close for i */

return soc;
}

/* BUBBLE SORT (USED IN RANKED PAIRS) */

// [[Rcpp::export]]

	IntegerMatrix bubbleSort(IntegerMatrix X, int col, int a=0)
	{
		/* a==1 for acending sort; otherwise descending sort (the default) */
		int n = X.nrow();
		IntegerMatrix M(clone(X));
		IntegerVector temp(M.ncol());				/* temp will be a row of M */

		if(a==1){									/* ascending sort */
  	      for(int i = 0; i < n-1; i++){
	        // Last i elements are already in place
	        for (int j = 0; j < n-i-1; j++){
	          if ( M(j,col) > M(j+1,col) ){
	            temp = M(j,_);
	            M(j,_) = M(j+1,_);
	            M(j+1,_) = temp;
	  	      }
		    };										/* close for i */
		  };										/* close for j */
		}
		else{										/* descending sort */
	      for(int i = 0; i < n; i++){
	        for (int j = i+1; j < n; j++){
	          if ( M(i,col) < M(j,col) ){
	            temp = M(i,_);
	            M(i,_) = M(j,_);
	            M(j,_) = temp;
	  	      }
		    };										/* close for i */
		  };										/* close for j */
		}											/* close if */

	return M;
    }


/* SMALLEST WITH OPEN GATE (USED BY BALDWIN). IDENTIFIES THE SMALLEST INTEGER IN COL 1 OF A MATRIX IF GATE IS OPEN */

// [[Rcpp::export]]
	int smallest_g(IntegerMatrix B, IntegerVector G){
		IntegerVector V(B.nrow(), NA_INTEGER);
		if( B.nrow() != G.length() ){
			Rcpp::stop("Error. In smallest_g(), the number of rows of B (Borda Scores) must equal the length of G (gate).");
		}

  	    for (int i = 0; i < B.nrow(); i++){
		  if( G(i)==1 ){						/* If the gate is open we will record the Borda score from B in V. */
			V(i) = B(i,1);
	  	  }
		};
	return min(na_omit(V));
    }


/* BALDWIN */

// [[Rcpp::export]]

	Rcpp::IntegerMatrix baldwin(Rcpp::IntegerVector F, Rcpp::IntegerMatrix perm) {
		int A = perm.ncol();					/* A is the number of alternatives.  P is the permutation matrix, perm.   */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, starting contents -9 (i.e., missing) */
		std::fill(soc.begin(), soc.end(), -9);
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
		int m;									/* the min borda count in a given round */
		IntegerMatrix BC(A,2);				 	/* borda count for alternative i: column 0 = alternative name, column 1 = borda score */
     	IntegerVector relim(A, A);				/* relim records the round an alt was eliminated with the vector index used as the alt name (initially A) */
    	IntegerVector G(A, 1);					/* G is a gate length A. G = 1 if alternative at that index is included in BC; 0= if eliminated  */
		int LessA = 0;

		/* run borda count and determine smallest borda score */
		int r = 1;								/* r is the round, used to keep track of round that each alternative is eliminated */
		do {									/* repeat until only one alternative remains */
		  BC = nborda(F,perm,G);				/* BC (Ax2) with col 0 = alternative name, col 1 = borda score */
		  m = smallest_g(BC,G);					/* determine the smallest borda score among all non-eliminated alternatives */
	      for (int i = 0; i < A; i++){			/* walk through each alternative */
		    if(BC(i,1) == m){					/* for all alternatives with the minimum score, we */

			  /* if not previously eliminated: record the round eliminated (value) in relim for each alternative (vector index), */
			  if( G(BC(i,0)) == 1 ){
			    relim( BC(i,0) ) = r;
			  }
		 	  G(BC(i,0)) = 0;					/* then mark the alternative as eliminated in the gate for the next round, */
			  LessA++;							/* and keep track of the number of alternatives eliminated so we know when to stop. */

		    }
		  };
		r++;
		} while(A - LessA > 1);					/* if more than one alternative remains, we do it again */


		/* assign pairwise social preferences based on larger r (i.e., later round eliminated) */
	    for (int i = 0; i < A; i++){
  	      for (int j = i; j < A; j++){
			if( relim(i) > relim(j) ){
			  soc(i,j)= 1;
			}
			if( relim(j) > relim(i) ){
			  soc(j,i)= 1;
			}
			if( relim(i) == relim(j) ){
			  soc(i,j)= 0;
			  soc(j,i)= 0;
			}
		  };
		};

	return soc;
	}

/* TRANSITIVITY EVALUATION (used in ranked pairs and as a criteria) */

// [[Rcpp::export]]
  int transitivity_eval(Rcpp::IntegerMatrix soc
	                   ) {
  int A = soc.nrow() ;							/* A is the number of alternatives */
  int contradiction = 0;

      for(int i=0; i < A; i++) {         		/* row index for first alternative in transitive condition */
         for(int j=0; j < A; j++) {         	/* column index for second alternative in transitive condition */
            if(soc(i,j) >= 0) {					/* if a>= b, i.e. P, R, or I, then ... */
              for(int k=0; k < A; k++) {  		/* for all possible c (the column index for third alternative) */
                 if(soc(j,k) >= 0) {			/* check whether if b>=c. If it is, we need to test for transitivity ... */

					switch(soc(i,k)) {
					  case -9 :
						contradiction=1;
						goto end_function;
					  case 0 :
					    if( (soc(i,j) != 0) | (soc(j,k) != 0) ){
						  contradiction=1;
					      goto end_function;
					    }
						break;
					  case 1 :
					    if( (soc(i,j) != 1) & (soc(j,k) != 1) ){
						  contradiction=1;
					      goto end_function;
						}
					  	break;
					return 0;
					}							/* close switch */

                  };    						/* close if (j,k) */
              };       							/* close for k */
            };         							/* close if (i,j) */
         };            							/* close for j */
      };               							/* close for i */

    end_function:
	return contradiction;						/* identified a three-tuple in an intransitive relationship (using goto) */
}


/* TRANSITIVITY UPDATE (USED IN RANKED PAIRS) */

// [[Rcpp::export]]

  Rcpp::IntegerMatrix transitivity_update(Rcpp::IntegerMatrix soc
	                   ) {
  int A = soc.nrow() ;							/* A is the number of alternatives */
  int Ttrue = 0;

  do{
    Ttrue = 0;                                  /* We are going through this matrix mulitiple times because we are updating cells
                                                   based on transitivity before being confident that transitivity is not violated */
      for(int i=0; i < A; i++) {         		/* row index for first alternative in transitive condition */
         for(int j=0; j < A; j++) {         	/* column index for second alternative in transitive condition */
            if(soc(i,j) >= 0) {					/* if a>= b, i.e. P, R, or I, then ... */
              for(int k=0; k < A; k++) {  		/* for all possible c (the column index for third alternative) */
                 if(soc(j,k) >= 0) {			/* check whether if b>=c. If it is, we need to test for transitivity ... */

					switch(soc(i,k)) {
					  case -9 :
					    if( soc(k,i) == -9 ){	/* in this case we can update soc */
						   Ttrue = 1;			/* indicates we are updating soc and need to loop through again */
						   if( (soc(i,j) == 0) & (soc(j,k) == 0) ){
						     soc(i,k) = 0;
						     soc(k,i) = 0;
						   }
						   else{
						     soc(i,k) = 1;
						     soc(k,i) = -9;
						   }
						}
						else{					/* a contradiction because soc(k,i) is not missing, even though soc(i,k) is missing */
						  goto end_function;
						}
						break;
					  case 0 :					/* in this case both antecedents should be indifferent or we have a contradiction */
					    if( (soc(i,j) != 0) | (soc(j,k) != 0) ){
					      goto end_function;
					    }
						break;
					  case 1 :					/* in this case at least one of the antecedents need to be strict or we have a contradiction */
					    if( (soc(i,j) != 1) & (soc(j,k) != 1) ){
					      goto end_function;
						}
					  	break;
					return 0;
					}							/* close switch */

                  };    						/* close if (j,k) */
              };       							/* close for k */
            };         							/* close if (i,j) */
         };            							/* close for j */
      };               							/* close for i */
   } while(Ttrue != 0);


    end_function:
	return (soc);								/* returns soc after all possible transitivity updates */
}

/* DUDE MISSING (USED IN RANKED PAIRS) */

// [[Rcpp::export]]

	int dude_missing(IntegerVector X)
	{
		if ( X.length() != 3 ){
		  Rcpp::stop("Yo, stop. dude_missing() is designed to input a vector length 3 for A=3.");
		}

		if( (X(0) == 0 & X(1) == 1) |  (X(0) == 1 & X(1) == 0) ){
		  return 2;
		}
		if( (X(0) == 0 & X(1) == 2) |  (X(0) == 2 & X(1) == 0) ){
		  return 1;
		}
		if( (X(0) == 1 & X(1) == 2) |  (X(0) == 2 & X(1) == 1) ){
		  return 0;
		}
   	}



/* RANKED PAIRS */

// [[Rcpp::export]]

Rcpp::IntegerMatrix ranked_pairs(Rcpp::IntegerVector F, Rcpp::IntegerMatrix P) {
		int n = P.nrow();						/* P is the permutation matrix, perm.  n is the number of strict linear orders, A! */
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix soc(A, A);				/* social preference between each binary pair, possibly returned */
		std::fill(soc.begin(), soc.end(), -9);	/* set off-diagonal elements in soc to -9 for missing */
		soc.fill_diag(0);						/* diagonal elements indifferent (because they are same alternatives) */
		IntegerMatrix socp(clone(soc));			/* tournament matrix which includes vote margins of i over j */
		std::fill(socp.begin(), socp.end(), 0);	/* make all elements in socp 0 */
		IntegerVector V(A);						/* total number of pairwise wins for alternative i */
		IntegerMatrix VM( (((A*A)-A)/2) , 3);	/* matrix of vote margins of i's pairwise wins over j's pairwise wins: col 0 = i, col 1= j, col 3= vote margins. */
		int winner;
		int dude;

		if ( A > 3 ){							/* to make ranked pairs apply to A>3 we must add a test for transitivity at KD below */
		  Rcpp::stop("The ranked pairs function is currently designed for no more than three alternatives!");
		}

		/* create tournament matrix by comparing each pair of alternatives head-to-head */
		int k=0;
	    for (int i = 0; i < A-1; i++){
  	      for (int j = i+1; j < A; j++){
			winner = head_to_head(F,P,i,j);

			if (winner == i){
			  VM(k,0) = i;
			  VM(k,1) = j;
			  VM(k,2) = abs(vote_margin(F,P,i,j));	/* vote margin is a positive number */
			}
			if (winner == j){
			  VM(k,0) = j;
			  VM(k,1) = i;
			  VM(k,2) = abs(vote_margin(F,P,i,j));	/* vote margin is a positive number */
			}
			if (winner == -8){						/* i and j tie */
			  VM(k,0) = i;
			  VM(k,1) = j;
			  VM(k,2) = 0;
			}
		  k++;
		  };										/* close for j */
		};											/* close for i */


		/* sort VM by vote margin (largest to smallest) */
		VM = bubbleSort(VM,2,0);

		/* Assign social preferences for all three pairs based on larger VM. */
	    for(int i = 0; i < A; i++){
			if( VM(i,2) > 0 ){
			  soc( VM(i,0), VM(i,1) )= 1;
			}
			if( VM(i,2) == 0 ){
			  soc( VM(i,0), VM(i,1) )= 0;
			  soc( VM(i,1), VM(i,0) )= 0;
			}
		};

		/* FOUR cases (note case 5, sent to Jac, was combined with case 2) */

		  /* 1- the pairwise relationships created by VM(0), VM(1), and VM(2) are transitive. */
		  if(transitivity_eval(soc) == 0){
		     return soc;
		  }
		  else{
		  /* 2- the pairwise relationships created by VM(0), VM(1), and VM(2) are intransitive and VM(0) >= VM(1) > VM(2). */
			 if( (VM(0,2) >= VM(1,2)) & (VM(1,2) > VM(2,2)) ){
		       /* action: lock in VM(0) and VM(1), already in soc, then apply transitivity to third pair. */
			   soc( VM(2,0), VM(2,1) )= -9;				/* erase preference among alternatives in VM(2) to assure transitivity updates correctly. */
			   soc( VM(2,1), VM(2,0) )= -9;
			   return transitivity_update(soc);
		 	 }
		  /* 3- the pairwise relationships created by VM(0), VM(1), and VM(2) are intransitive and VM(0) > VM(1) = VM(2). */
			 if( (VM(0,2) > VM(1,2)) & (VM(1,2) == VM(2,2)) ){
		       /* action: lock in VM(0), make third alternative indifferent to the winner of VM(0), according to Tideman quotes. */
			   std::fill(soc.begin(), soc.end(), -9);	/* erase all preferences previously recorded in soc */
			   soc.fill_diag(0);
			   soc( VM(0,0), VM(0,1) )= 1;				/* lock in VM(0) */
			   dude= dude_missing(VM(0,_));				/* identify the alternative not in the VM(0) pair.  Call him dude. */
			   soc( VM(0,0), dude )= 0;					/* make dude indifferent to VM(0,0) */
			   soc( dude, VM(0,0) )= 0;
			   soc( dude, VM(0,1) )= 1;					/* make dude preferred to VM(0,1) */
			   return soc;
		 	 }
		  /* 4- the pairwise relationships created by VM(0), VM(1), and VM(2) are intransitive and VM(0) = VM(1) = VM(2). */
			 if( (VM(0,2) == VM(1,2)) & (VM(1,2) == VM(2,2)) ){
		       /* action: report indifference among all three alternatives */
			   std::fill(soc.begin(), soc.end(), 0);	/* that makes all values in soc 0s */
			   return soc;
			 }
		  }												/* close else */
}



/* ------------------------------------------- FUNCTIONS THAT EVALUATE WEAK CONDORCET WINNERS ---------------------------------------------------- */

/* WHICH_MAX (used with wCW() ) */

// [[Rcpp::export]]

   IntegerVector which_max(NumericVector v) {
    double current_max = v[0];
    int n = v.size();
    std::vector< int > res;
    res.push_back( 0 );
    int i;

    for( i = 1; i < n; ++i) {
      double x = v[i];
      if( x > current_max ) {
        res.clear();
        current_max = x;
        res.push_back( i );
      } else if ( x == current_max ) {
        res.push_back( i );
      }
    }
    Rcpp::IntegerVector iv( res.begin(), res.end() );
   return iv;
  }

/* wCW() -- identifies Weak Condorcet Winners (only works for 3 transitive alternatives) */

// [[Rcpp::export]]

   IntegerVector wCW(IntegerMatrix soc) {
	int A = soc.ncol();						/* soc is a matrix of pairwise majority rule relationships.  A is the number of alternatives */
	IntegerVector S(A);						/* column sums of soc */
	IntegerVector V(A,0);					/* vector indicating whether alternative is a wCW */

	   S = rowSumsC(soc);

	/* determine which alternatives are wCW and store in V (determine some types as well) */
	   for(int i = 0; i < A; i++){
	     if(S(i)==2){						/* i is the only wCW (type I or II) */
	       V(i) = 1;
	       break;
	     }
	     if(S(i)==0){						/* all three are wCW.  Type IV. */
	       V = {1,1,1};
	     }
	     if(S(i)==1){						/* could be type II or III */
	       if(max(V)==2){					/* type II */
	         V(which_max(S)) = 1;			/* assign 1 to index location with 2 */
	         break;
	       }
	       else{							/* type III */
	         for(int j = 0; j < A; j++){	/* assign 1 to two index locations with 1 */
	           if(S(j) == 1){V(j)=1;}
		     };
	         break;
	       }								/* close else */
	     }									/* close if(S(i)==1) */
	   };									/* close for loop */

   return(V);
   }

/* wCW_compare - compares decision set to wCW */

// [[Rcpp::export]]

   int wCW_compare(IntegerMatrix soc, IntegerMatrix socg) {
	  int A = soc.ncol();						/* soc is created by a PAR; socg is created by pmr */
	  IntegerVector V = wCW(soc);				/* applied to soc, wCW() determines the social decision set from the PAR */
	  IntegerVector W = wCW(socg);				/* applied to gsoc, wCW() finds the weak Condorcet winner(s) */
	  int tally = 0;							/* the number of times a PAR selects an alternative that is also a wCW */

      for(int i = 0; i < A; i++){
	    if((V(i)==1) & (W(i)==1)){				/* count the number of times the PAR picked "a" weak Condorcet winner */
	      tally++;
	    }
      };
      if(sum(V) == tally){						/* if those were the only alternatives chosen by the PAR, then the PAR picks wCWs */
		return(1);
      }

   return(0);									/* otherwise the PAR does not pick weak Condorcet winners in this case */
   }


/* wCW_type() -- determines Weak Condorcet Winner type for any type (only works for three alternatives, trans or intrans) */

// [[Rcpp::export]]

   int wCW_type(IntegerMatrix soc) {
	int A = soc.ncol();						/* soc is a matrix of pairwise majority rule relationships.  A is the number of alternatives */
	IntegerVector S(A);						/* column sums of soc */
	int SS;									/* the sum of vector S */
	int type=0;								/* one of four wCW types (see paper) */

	/* sum the columns of soc.  Sums will indicate wCW type and winners as follows */
	   /* 	(2) won both;				(-8) won one, lost one; 	(1) won one, tied one;
	      	(-9) tied one, lost one;	(0) tied both; 				(-18) lost both.
	   */
	   S = rowSumsC(soc);
	   SS = sum(S);

	/* determine the wCW type using values in S and SS */
	     if(SS== -24 & min(S)== -18){				/* type 1 */
	       type = 1;
	     }
	     if((SS== -16) & (max(S)==2)){					/* type 2 */
	       type = 2;
	     }
	     if((SS== -16) & (min(S)== -18)){				/* type 3 */
	       type = 3;
	     }
	     if(SS== 0){									/* type 4 */
	       type = 4;
	     }
	     if((SS== -16) & (max(S)==1) & (min(S)== -9)){	/* type 5 */
	       type = 5;
	     }
	     if(SS== -8){									/* type 6 */
	       type = 6;
	     }
	     if(SS== -24 & min(S)== -8){					/* type 7 */
	       type = 7;
	     }

	/* warning */
	if(type==0){
	   Rcpp::stop("Error: Type was not properly assigned by wCW_type().");
 	}

   return(type);
   }


/* wCW Type B (only works for 3 transitive alternatives).  Cannot differentiate types V & VI using this method.  Used to test wCW_type() */

// [[Rcpp::export]]

   int wCW_typeB(IntegerVector F, IntegerMatrix P) {
		int A = P.ncol();						/* P is the permutation matrix, perm.  A is the number of alternatives */
		IntegerMatrix soc(A, A);				/* tournament matrix which is just 1's, 1.5's, and 0's */
		NumericMatrix socp(A, A);				/* tournament matrix which is just 1's, 0's, and .5's */

		std::fill(socp.begin(), socp.end(), 0);	/* make all elements in socp 0 */
		NumericVector V(A);						/* total number of pairwise wins for alternative i */
		int winner;
		int type = 0;

		/* create tournament matrix by comparing each pair of alternatives head-to-head */
	    for (int i = 0; i < A-1; i++){
  	      for (int j = i+1; j < A; j++){
			winner = head_to_head(F,P,i,j);

			if (winner == i){
			  socp(i,j)= 1;
			}
			if (winner == j){
			  socp(j,i)= 1;
			}
			if (winner == -8){					/* i and j tie */
			  socp(i,j)= .5;
			  socp(j,i)= .5;
			}
		  };
		};

		/* sum the rows of the tournament matrix */
		V = rowSumsN(socp);

  		/* TYPES AND COPELAND VALUES
  		   I.   A:2,   B:1,   C:0
  		   II.  A:2,   B:.5,  C:.5
   		   III. A:1.5, B:1.5, C:0
  		   IV.  A:1,   B:1,   C:1
		*/

		/* evaluate the type of wCW using V */
        	   if((max(V) == 2) & (min(V) == 0)){
    		     type = 1;
    		   }
    		   if((max(V) == 2) & (min(V) == .5)){
    		     type = 2;
    		   }
    		   if(max(V) == 1.5){
        	     type = 3;
        	   }
    		   if(max(V) == 1){
    		     type = 4;
    		   }

   return(type);
   }


/* ITERATE: THE PARENT FUNCTION */

// [[Rcpp::export]]

	Rcpp::List iterate(IntegerMatrix perm, int N, int T, int dt, NumericVector V) {
	  int n = perm.nrow();						/* perm is the permutation matrix; n is the number of strict linear orders, A! */
	  int A = perm.ncol();						/* A is the number of alternatives */
  	  IntegerVector F(n,0) ;					/* F is vector, length n, that records the number of individuals with the respective row preference in P. */

	  /* "count" tallies the total number of times the PAR selects a wCW among each of the four transitive wCW types */
	  IntegerVector count_prr(4,0);
	  IntegerVector count_aprr(4,0);
	  IntegerVector count_borda(4,0);
	  IntegerVector count_nanson(4,0);
	  IntegerVector count_hare(4,0);
	  IntegerVector count_pmr(4,0);
	  IntegerVector count_copeland(4,0);
	  IntegerVector count_coombs(4,0);
	  IntegerVector count_dowdall(4,0);
	  IntegerVector count_sk(4,0);
	  IntegerVector count_baldwin(4,0);
	  IntegerVector count_rp(4,0);
	  int count_coombs_rp = 0;					/* counts the number of times coombs and ranked pairs produce different orders */
	  int count_intransitive = 0;				/* total number of intransitive rankings across all iterations (we skip these) */
	  IntegerVector typeA(7,0);					/* tally of wCW types as determined by wCW_type() */
	  IntegerMatrix socg(A,A);				    /* matrix of pairwise majority rule relationships -- g for goal */
	  List ret;
	  List TBL;
	  int TA, TB;								/* indicates type, measured in two different ways */
	  int times_through=0;			/* TEST */

 	  if(dt < 0 || dt > 2 || A != 3){
	    Rcpp::stop("Error: It must be the case that A=3 and dt = 0, 1, or 2 (note: distribution type (dt) =0 if IC; =1 if IAC; = 2 if weighted).");
 	  }

	/* --- BEGIN TRIAL (i.e. ITERATION) LOOP: TRIAL t IN THE PROBABILITY EXPERIMENT --- */

	 for(int t=1; t <= T; t++){

		/* DRAW PREFERENCES */
     	   if(dt == 0) {									/* if IC preferences requested, draw IC on perm */
		     F=icc(perm, N);
		   }
    	   if(dt == 1) {									/* if IAC preferences requested, draw IAC on perm */
		     F=iac(perm, N);
		   }
    	   if(dt == 2) {									/* if weighted preferences requested, draw preferences using weights V */
		     F=weighted_pref(perm, N, V);
		   }

		/* DETERMINE PAIRWISE-MAJORITY RELATIONSHIPS */
		   socg=pair_maj_rule(F, perm);

		/* DETERMINE THE wCW TYPE (any type) AND INCREMENT COUNTER */
	   	   TA = wCW_type(socg);
  		   switch( TA ){
  			case 1:
  			  typeA(0)++;							/* type 1 */
  			  break;
  			case 2:
  			  typeA(1)++;							/* type 2 */
  			  break;
  			case 3:
  			  typeA(2)++;							/* type 3 */
  			  break;
  			case 4:
  			  typeA(3)++;							/* type 4 */
  			  break;
  			case 5:
  			  typeA(4)++;							/* type 5 */
  			  break;
  			case 6:
  			  typeA(5)++;							/* type 6 */
  			  break;
  			case 7:
  			  typeA(6)++;							/* type 7 */
  			  break;
  		   }

		/* EVALUATE TRANSITIVITY OF PAIRWISE RELATIONSHIPS (1 = intransitive) */
	       if(transitivity_eval(socg)==1){
			   count_intransitive = count_intransitive + 1;
			   goto end_iteration;					/* if intransitive, count case and move to the next iteration */
	   	   }

		/* DETERMINE THE wCW TYPE AGAIN (done a second way for comparison, works for transitive relations only) */
	   	   TB = wCW_typeB(F,perm);
 	       if( TA != TB ){
	  	     Rcpp::stop("Error: wCW_type() and wCW_typeB() accessed different wCW types from the same transitive profile.");
 	       }

		/* DETERMINE WHETHER THE PAR SELECTED A WEAK CONDORCET WINNER AND ADD 1 OR 0 TO THE COUNT FOR THAT PAR BY TYPE (wCW_compare = 1 if decision set is wCW; = 0 otherwise) */
		   count_prr(TA-1) = count_prr(TA-1) + wCW_compare(prr(F,perm),socg);							  	/* plurality ranking rule */
		   count_aprr(TA-1) = count_aprr(TA-1) + wCW_compare(aprr(F,perm),socg);						  	/* anti-plurality ranking rule */
		   count_borda(TA-1) = count_borda(TA-1) + wCW_compare(borda(F,perm),socg);						  	/* borda count */
		   count_nanson(TA-1) = count_nanson(TA-1) + wCW_compare(nanson(F,perm),socg);					  	/* nanson's rule */
		   count_hare(TA-1) = count_hare(TA-1) + wCW_compare(hare(F,perm),socg);						  	/* hare's rule */
		   count_pmr(TA-1) = count_pmr(TA-1) + wCW_compare(pair_maj_rule(F,perm),socg);						/* pairwise majority rule */
		   count_copeland(TA-1) = count_copeland(TA-1) + wCW_compare(copeland(F,perm),socg);				/* copeland's rule */
		   count_coombs(TA-1) = count_coombs(TA-1) + wCW_compare(coombs(F,perm),socg);				  	 	/* coombs method */
		   count_dowdall(TA-1) = count_dowdall(TA-1) + wCW_compare(dowdall(F,perm),socg);					/* dowdall method */
		   count_sk(TA-1) = count_sk(TA-1) + wCW_compare(sk(F,perm),socg);							 		/* simpson-kramer method */
		   count_baldwin(TA-1) = count_baldwin(TA-1) + wCW_compare(baldwin(F,perm),socg);					/* baldwin method */
		   count_rp(TA-1) = count_rp(TA-1) + wCW_compare(ranked_pairs(F,perm),socg);  						/* ranked pairs method */

		/* COMPARE OUTPUT FROM COOMBS AND RANKED PAIRS */
   		   if( sum(coombs(F,perm)==ranked_pairs(F,perm)) != (A*A) ){			/* If the outputs are not identical */
             count_coombs_rp++;
		   }

     end_iteration:
     times_through++;
	 };
	/* --- END TRIAL LOOP --- */


	ret["dt"] = dt;
	ret["F"] = F;
	ret["perm"] = perm;
	ret["num_intransitive"] = count_intransitive;
	ret["times_through"] = times_through;
	ret["type"] = typeA;
	ret["count_prr"] = count_prr;
	ret["count_aprr"] = count_aprr;
	ret["count_borda"] = count_borda;
	ret["count_nanson"] = count_nanson;
	ret["count_hare"] = count_hare;
	ret["count_pmr"] = count_pmr;
	ret["count_copeland"] = count_copeland;
	ret["count_coombs"] = count_coombs;
	ret["count_dowdall"] = count_dowdall;
	ret["count_sk"] = count_sk;
	ret["count_baldwin"] = count_baldwin;
	ret["count_rp"] = count_rp;
	ret["contra_coombs_rp"] = count_coombs_rp;

	return(ret);
	}


