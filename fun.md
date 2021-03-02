(Bayesian?) Statistics with the Symmetric Group
===============================================

 

We are interested in a finite set of things:

General Case: N = {1,2,…,n-1,n}

Binary Example: M = {Match, Non-Match}

Real-World Example: C = {AUS, BRA, ESP, FJI, FRA, IDN, ITA, JPN, NZL, PRT, USA,
ZAF}

 

Number of things we are interested in:

\|N\| = n

\|M\| = 2

\|C\| = 12

 

Our Predictors/Regressors/Features are:

 

Discrete:

| Feature  | Range of Feature                                                                                                              | dim |
|----------|-------------------------------------------------------------------------------------------------------------------------------|-----|
| evtYears | {2017, 2018, 2019}                                                                                                            | 3   |
| evtOrig  | {AUS, BRA, FRA, IDN, PRT, USA, ZAF}                                                                                           | 7   |
| evtName  | {Bali Pro, Bells Beach, France, Gold Coast, J-Bay Open, Margaret River, Peniche Pro, Pipe Masters, Rio Pro, Tahiti, Trestles} | 11  |
| rnd      | {0,1,2,3,4,5,6,7}                                                                                                             | 1   |
| heat     | {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}                                                                                      | 16  |
| athOrig  | {AUS, BRA, FRA, IDN, ITA, JPN, NZL, PRT,USA,ZAF}                                                                              | 10  |
| I_match  | {1,0}                                                                                                                         | 2   |

Note: Questions:

-   heat variable as categorical —\> dim=16 vs. heat as ordinal/cts —\> dim=1 ?

-   rnd variable is ordinal, so dim should be 1, right?

    -   values are very much not independent

-   I_match := 1 if athOrig in panelOrigs, else 0

    -   Do we need a panelOrigs variable?

 

Cts or Integer-Valued:

| Feature       | Range of Feature | dim |
|---------------|------------------|-----|
| currentPoints | [0, 57855]       | 1   |
| endingPoints  | [265, 62490]     | 1   |
| waveInHt      | [3,42]           | 1   |
| actualSco     | [0, 10]          | 1   |

Note 1: I should pick a more indicative name than “actualSco”, perhaps “wslSco”.
This is the score the WSL determines (it is the trimmed mean of 5 judge scores).

Note 2: we get 3 or 5 sub scores (usually 5) and while “actualSco” (or “wslSco”)
is technically an observed variable, we could still transform subScos in various
ways. Basic ideas: min( ), max( ), mean( ), and var( ). An interesting idea (I
think) is to use any of the basic ideas but applied to map( sco -\> min(sco,
10-sco), subScos) because any sort of distribution over rankings will likely be
tempered when closer to the boundaries of score range (in some sense there is
less subjectivity, but also simply less variance (is this truncation?) ).

 

 

dims are: (evtYear, evtOrig, evtName, rnd, heat, athOrig, I_match, waveInHt,
actualSco)

In conclusion we need tensor with dims = (3,7,11,1,16,10,2,1,1)

 

Our Outcome/Response variables are:

|                         | Feature    | Range of Feature                      | empirical max dim | max dim |
|-------------------------|------------|---------------------------------------|-------------------|---------|
| Restriction to Subspace | panelOrigs | subset of C with less than 5 elements | 5                 | 12      |
| Point Dist. on Subspace |            | subsets of S_panelOrigs               | 5x5               | 12x12   |
| Point Dist.             |            | subsets of S_C                        | n/a               | 12x12   |

 

We can construct alternative Outcome/Response variables:

Given Point Dist. for Data ( P_i )_i=1, we may be interested in how point
estimates change over an interval (MOVING FORWARD)

i.e. we are curious about F, where F P_i= P_{i+1}

Lucky for us, F_i are invertible, so: F =P_{i+1}P_i\^{-1}

 

We can also move backwards (Get from “A model for a data sequence”)

 

Furthermore, if we are interested in how this changes over larger intervals, say
of size k, let:

F_k P_i = P_{i+k} —\> F_k = P_{i+k}P_i\^{-1}

B_k P_{i+k} = P_i —\> B_k = P_i P_{i+k}\^{-1}

 

Hypothesis Test:

H0,0: EmpDist on Subspace = Unif

H1,0: EmpDist on Subspace != Unif

 

H0,1: EmpDist = Unif

H1,1: EmpDist != Unif

 

P(D\|Model)

 

 

 

 
