Statistics Using The Symmetric Group
====================================

 

Data
----

We have some data from WSL Surfing Competitions comprised of 21,615 waves surfed
in the 2017, 2018, and 2019 World Surf League Championship Tour.

 

The rules of the World Surf League Championship Tour were modified between each
season. The 2019 rules are as follows:

Chapter 1 Article 3.01: "CT events are limited to 13 per year and a limited
number of events in any one country as decided by the Office of Tours and
Competition. For this purpose Tahiti, Reunion, Hawaii, and other places as
determined by the Office of Tours and Competiion are classified as countries."

Chapter 1 Article 7.01: Mens CT events, consist of 36 surfers, in the following
format:

-   Round 1 has 12 Heats of 3 Surfers where

    -   1st & 2nd Place —\> Round 3

    -   3rd Place —\> Round 2

-   Round 2 has 4 Heats of 3 Surfers where

    -   1st & 2nd Place of each heat —\> Round 3

    -   3rd Place —\> exit

-   Round 3 has 16 Heats of 2 Surfers

    -   1st Place —\> Round 4

    -   2nd Place —\> exit

-   Round 4 has 8 heats of 2 surfers where

    -   1st Place —\> Quarter Finals

    -   2nd Place —\> exit

-   Quarter Finals has 4 heats of 2 Surfers where

    -   1st place —\> Semi Finals

    -   2nd place —\> exit

-   Semi Finals has 2 Heats of 2 Surfers where

    -   1st place —\> Final

    -   2nd place —\> exit

-   Final has 1 Heat of 2 Surfers where

    -   1st place —\> Final

 

Chapter 13 Article 179: Judging Panel Composition

179.01 For mens CT events, there will be 1 International Head Judge, 7
International Judges, and 1 International Priority Judge.

179.13 The WSL Head Judge is responsible for assuring that a minimum of 5 judges
sit on the panel for every heat of all CT Events. These 5 judges must be subset
of the 7 International Judges and the 1 WSL Head Judge.

179.17 "At CT Events, the number of Judges from any one Regional area is limited
to 3"

 

We are interested in some finite sets:

Binary Example:

$$
B = \{ Match, Non-Match} = { M, \neg M }
$$

Country Example:

C = {AUS, BRA, ESP, FJI, FRA, IDN, ITA, JPN, NZL, PRT, USA, ZAF}

General Case: Dimensional = D = { 1, 2, … , d-1, d}

 

Number of things we are interested in:

Binary Example: \|B\| = 2

Country Example: \|C\| = 12

General Case: \|D\| =

 

Information is Discrete:

| Feature    | Range of Feature                                                                                                              | dim | Notes                               |
|------------|-------------------------------------------------------------------------------------------------------------------------------|-----|-------------------------------------|
| evtYears   | {2018, 2019}                                                                                                                  | 3   |                                     |
| evtOrig    | {AUS, BRA, FRA, IDN, PRT, USA, ZAF}                                                                                           | 7   |                                     |
| evtName    | {Bali Pro, Bells Beach, France, Gold Coast, J-Bay Open, Margaret River, Peniche Pro, Pipe Masters, Rio Pro, Tahiti, Trestles} | 11  |                                     |
| rnd        | {1,2,3,4,5,6,7,8}                                                                                                             | 8   | 1 is final, everything else follows |
| heat       | {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}                                                                                      | 16  |                                     |
| athOrig    | {AUS, BRA, FRA, IDN, ITA, JPN, NZL, PRT,USA,ZAF}                                                                              | 10  |                                     |
| panelOrigs | subsets of {AUS, BRA, ESP, FRA, PRT, USA, ZAF}                                                                                |     |                                     |
| I_match    | {1,2}                                                                                                                         | 2   | 1 if athOrig in PanelOrigs, else 2  |
|            |                                                                                                                               |     |                                     |
|            |                                                                                                                               |     |                                     |

Within Each year, there exists:

-   10 Events:

-   1 WSL Head Judge

For each Event,

-   The WSL provides 7 Judges (call it EvtJudges)

-   There are 6 to 8 Rounds

For each Round

-   There are some Number of Heats

For each Heat

 

EvtJudgesOrigs = union over j in EvtPanel Origs(j)

Range of \#(EvtJudgesOrigs) is {3,4,5,6}

```
[ union([x.panelOrigs for x in filter(w-> w.evtId == evt, waves)]...) for evt in EvtIds ]


union(A::Array{Multiset{T}}) where {T} = 
keyset = union(map(x->keys(x.data),A)...)
return Multiset([repeat(])




reduce(union, Set(A) )
```

 

which is hosted in someplace in some country, sometimes even two countries for
one event

there are up to 8 rounds

and up to 16 heats in each round

usually at some point during a heat, a surfer will ride a wave

there are five judges on a panel that rate the surfers ride

 

For each heat, there is a set judges, each of which has a nationality (ie.
country of origin).

Panel

 

We have some additional information:

The nationality of a surfer

that batuibakutt of each judge

 

 

Panels = { S \\in 2\^({:AUS, :BRA, :ESP, :FRA, :PRT, :USA, :ZAF}) \|
\#(Union(S)) =5 }

Panels

([0:0.1:10]\^(Ctry) )\^5

([0:0.1:10]\^(Ctry) ) x ([0:0.1:10]\^(Ctry)) x ([0:0.1:10]\^(Ctry)) x
([0:0.1:10]\^(Ctry)) x ([0:0.1:10]\^(Ctry))

([0:0.1:10] x [0:0.1:10] x [0:0.1:10] x [0:0.1:10] x [0:0.1:10] ) \^ ( Ctry x
Ctry x Ctry x Ctry x Ctry )

 

(0:0.1:10\^Ctry same as maps Ctry —\> 0:0.1:10 )

Note that 0:0.1:10 is a totally ordered set, ie. we have (\<,0: 0.1:10)

 

(\<,0: 0.1:10)\^Panels

```
for k in 1:20 print("hello") end
```

 

Some Stats:

histogram(map(w-\>length(w.panelOrigs), waves))

 

A variable: \|PanelOrigs\|

Mean: 3.907379612257661

Var: 0.5189386620233642

SD: 0.7203739737270942

 

```
sum(map(w->(length(w.panelOrigs)-3.907379612257661)^2, waves))/(length(waves))
```

 

videofile("ExponetiationOfError.mp4") do

       for k in 0:10 

       histogram(map(w-\>abs(length(w.panelOrigs)-3.907379612257661)\^k /
factorial(k), waves) )

       title("Exponentiation of Error (m≈3.9074): (\|{Nationality of Judge s.t.
Judge on Panel}\| - m)\^ \$k / k!")

       sleep(1)

       for wait in 1:10 draw(gcf()) end

       end

       end

 

for k in 0:10 

       histogram(map(w-\>(length(w.panelOrigs)-3.907379612257661)\^k /
factorial(k), waves) )

       title("Exp of Error (m≈3.9074): (\|{Nationality of Judge s.t. Judge on
Panel}\| - m)\^ \$k / \$(k)! ")

       hold(true)

       for wait in 1:10 draw(gcf()) end

       end

 

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

 

 

 

 
