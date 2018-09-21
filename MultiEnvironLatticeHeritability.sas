*SAS program to estimate heritability from random lines evaluated in an incomplete block
(lattice) design with 3 replications at each of 5 environments, with some missing data;
data one;
filename rust DDE "EXCEL|[OatYield.xls]yielddata!R2C1:R2161C6";
Infile rust NOTAB DLM= '09'X DSD MISSOVER lrecl = 10240;
input env $ rep block entry line $ yield;

*eliminate checks and lines from different population from data set;
data two; set one; if entry < 96;

proc sort; by env;
proc print; run;

proc mixed asycov; class env rep block entry; model yield = ; 
random env rep(env) block(rep*env) entry entry*env;
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
*order of variance components in v and c matrices is V(E), V(R), V(B), V(G), V(GE), V(error);
LG = {0, 0, 0, 1, 0, 0};
LP = {0, 0, 0, 1, 1, 1};
call seh(V, C, LG, LP, H, SE);
print "Heritability on a Plot Basis", H, SE;
*based on the data structure, we computed that the harmonic mean of environments in which
lines were tested was 4.9, and the harmonic mean of total replications across all environments
per line was 14.7;
e = 4.9;
r = 14.7;
LP = 0//0//0//1//(1/e)//(1/r);
call seh(V, C, LG, LP, H, SE);
print "Heritability on a Family Mean Basis", H, SE;
quit; run;
