*SAS program to estimate heritability from random lines evaluated in an incomplete block
(lattice) design with 3 replications at one environment, with some missing data;
data one;
filename rust DDE "EXCEL|[OatYield.xls]yielddata!R2C1:R2161C6";
Infile rust NOTAB DLM= '09'X DSD MISSOVER lrecl = 10240;
input env $ rep block entry line $ yield;

*eliminate checks and lines from different population from data set;
data two; set one; if entry < 96;

*eliminate data from all environments except one, AME96;
data three; set two; if env = "96AME";

proc print; run;

proc mixed asycov; class rep block entry; model yield = ; 
random rep block(rep) entry;
ods listing exclude AsyCov CovParm; ods output asycov = covmat covparms = estmat;
proc iml;
start seh(V, C, LG, LP, H, SE);
Vp = LP`*V;
Vg = LG`*V;
H = VG/Vp;
d = (1/Vp)*(LG - (LP*H));
VH = d`*C*d;
SE = sqrt(VH);
finish seh;

use estmat; read all into v; use covmat; read all into c;
* Note that SAS introduces an extra first column into the matrix which must be removed;
C = C(|1:nrow(C), 2:ncol(C)|);
*order of variance components in v and c matrices is V(R), V(B), V(G), V(error);
LG = {0, 0, 1, 0};
LP = {0, 0, 1, 1};
call seh(V, C, LG, LP, H, SE);
print "Heritability on a Plot Basis", H, SE;
*based on the data structure, we computed that the harmonic mean of reps in which
lines were tested was 2.9;
r = 2.9;
LP = 0//0//1//(1/r);
call seh(V, C, LG, LP, H, SE);
print "Heritability on a Family Mean Basis", H, SE;
quit; run;
