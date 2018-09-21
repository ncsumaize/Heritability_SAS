* a sas program to implement estimate parent-offspring regression of F3 progeny means 
on individual F2 parents and simultaneously estimate genotypic and phenotypic variance
components of F2:3 families - data from Holland et al. (1998) Theoretical and Applied
Genetics 96:232-241;

*program written by Jim Holland, USDA-ARS, NCSU.  I give big props to an anonymous reviewer
at Crop Science of a previous manuscript who suggested the use of the iml module and a 
simpler way to code the standard error estimation routine! Thanks, person - whoever you are!;


options nocenter ls = 76;

*input data from excel spreadsheet available at www4.ncsu.edu/~jholland/heritability.html;
data one;
filename rust DDE "EXCEL|[CornRustData.xls]CornRustData!R2C1:R525C5";
Infile rust NOTAB DLM= '09'X DSD MISSOVER lrecl = 10240;
input env $ rep gen $ family rust;

proc print; run;

*first conduct standard ANOVA and regression analysis to obtain initial estimates
of variance and covariance parameters to input into the mixed models analysis;

* do regression analysis of F3 family means on F2 parents;
data f91; set one; if env = "F91";
rustF2 = rust; keep family rustf2;proc sort; by family;
data c91; set one; if env = "C91";
rustF3 = rust; keep family rustF3;proc sort; by family;
data poreg; merge f91 c91; by family;
proc reg; model rustF3 = rustF2;

*estimate total variance among F2 plants;
data f2; set one; if gen = "F2";
proc univariate; var rust;

* estimate genetic variance among F3 families using both ANOVA (proc glm) and
mixed models (proc mixed) as a comparison;
data f3; set one; if gen = "F3";
proc glm; class env rep family;
model rust = env rep(env) family family*env;

proc mixed; class env rep family;
model rust = ;
random env rep(env) family family*env;
run;

* use proc mixed to estimate covariance between F2 & F3, variance among F2, 
and genetic and phenotypic variances of F3 families;
proc mixed data = one asycov;class env rep family gen;
model rust = ;

*the macro-environment and block within environment effects are treated as random variables 
each with a single variance component;
random env rep(env) ;

*unique family - by - environment interaction variances are modeled for the two different 
generations by using the group option;
random family*env/group = gen;

*unique family variances are modeled for the two different generations, and the covariance 
between families with the same code is modeled by using the subject option and specifying an 
unstructured covariance matrix and specifying that each family has the same variance-covariance
structure;
random gen/subject=family type = un;

*unique error variances are modeled for the two different generations with the repeated 
command and group option;
repeated /group = gen;

* initial values of variance-covariance parameters based on the preliminary ANOVAs are 
introduced with the parms command - in the order that effects are specified in the random 
statements (VE, VR(E), VFE(S0),  VFE(S1), VF(S0), VG(S0,S1), VF(S1), Ve(S0), Ve(S1)).
The variance components for GE and Error within the F2 generation (VFE(S0) Ve(S1)) are forced
to be zero with the hold option - this is necessary because only one variance component is 
estimable in the F2 generation, as data were taken on individual plants, so the component 
VF(S0) is actually the phenotypic variance in the S0 generation;
parms (0.1451) (0.0057) (0) (0.146) (4.274) (2.5741) 
(3.0086) (0) (0.7515)/ hold=3,8;

*output estimates of variance and covariance component into data set "estmat". 
output variance-covariance matrix of the estimates into data set "covmat";
* to suppress the printing of estmat and covmat, use the command "
ods listing exclude asycov covparms";
ods output asycov = covmat covparms = estmat;

*start iml routine and create a module called "seh" that will return an estimate of heritability
and the Standard Error of Heritability estimate when given C and V matrices from
the proc mixed output and Lp and Lg vectors as defined in Holland and Cervantes-Martinez, 2001;
proc iml;
start seh(V, C, LG, LP, H, SE);
Vp = LP`*V;
Vg = LG`*V;
H = VG/Vp;
d = (1/Vp)*(LG - (LP*H));
VH = d`*C*d;
SE = sqrt(VH);
finish seh;

*make matrices v and c, corresponding to estmat and covmat data sets;
use estmat; read all into v; use covmat; read all into c;
* Note that SAS introduces an extra first column into the matrix which must be removed;
C = C(|1:nrow(C), 2:ncol(C)|);
*Note carefully the order of variance components in v and c matrices:VE, VR(E), VFE(S0), VFE(S1),
VF(S0), VG(S0,S1), VF(S1), Ve(S0), Ve(S1).  The vector of variance components estimates that is
of interest for estimating heritability based on S0:1 line variances includes only VF(S1), 
VFE(S1), and Ve(S1) - these are the 7th, 4th, and 9th components of the V and C matrices, 
respectively;
v = v(|7|)//v(|4|)//v(|9|);
c = c(|{7 4 9}, {7 4 9}|);
LG = {1, 0, 0};
LP = {1, 1, 1};
call seh(V, C, LG, LP, H, SE);
print "Heritability on a Plot Basis", H, SE;

*the harmonic mean of the number of plots per S0:1 family is 2.42 and the number of 
environments in which each family was tested is 1.23 -these are the divisors of error and
GE variance components, respectively, in the denominator of the heritability estimate 
(see Holland and Cervantes-Martinez, 2001);
eh = 1.23;
ph = 2.42;
lp = 1//(1/eh)//(1/ph);
call seh(V, C, LG, LP, H, SE);
print "Heritability on a Family Mean Basis", H, SE;

*now create a new pair of v and c matrices to estimate heritability from parent offspring 
regression.  In this case the variance components of interest are the 6th and 5th, respectively: 
VG(S0,S1) and V2F(S0);
use estmat; read all into v; use covmat; read all into c;
v = v(|6|)//v(|5|);
C = C(|1:nrow(C), 2:ncol(C)|);
c = c(|{6 5}, {6 5}|);
LG = {1, 0};
LP = {1, 1};
call seh(V, C, LG, LP, H, SE);
print "Heritability from regression of S1 offspring on individual parents", H, SE;
quit;
run;