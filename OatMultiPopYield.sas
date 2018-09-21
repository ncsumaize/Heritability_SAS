options ps=5000 ls=100  nonumber nodate nocenter;
DATA A;
FILENAME DATA1 DDE "EXCEL|[OatMultiPopYield.xls]YieldData!R3C1:R1630C10";
INFILE DATA1 NOTAB DLM= '09'X DSD MISSOVER lrecl = 10240;
INPUT ENV $ REP	SET	BLOCK ENT	CROSS $ LINE CYCLE $ BIOMASS YIELD;

proc sort; by cycle;run;


*estimate variance components for each cycle indpendently, ignoring incomplete blocks;

PROC MIXED COVTEST data = a;by cycle;
CLASSES ENV SET REP ENT ;
MODEL YIELD = ;
RANDOM ENV SET SET*ENV REP(SET ENV) ENT(SET) ENV*ENT(SET);
run;


data c; set a; if cycle = "C0" or cycle = "C1";
proc print; run;
PROC MIXED asycov data = c;
CLASSES ENV REP SET ENT CYCLE ;
MODEL YIELD = CYCLE/ddfm = satterth;
RANDOM ENV SET REP(ENV*SET) CYCLE*ENV SET*ENV; 
random ENT(SET) ENV*ENT(SET)/group = cycle;
LSMEANS CYCLE/PDIFF;
ods output asycov = covmat covparms = estmat;
run;


proc iml;
start seh(V,C,LG,LP,H,SE);
Vp = LP`*V;
Vg = Lg`*V;
H = Vg/Vp;
d = (1/Vp)*(LG - (Lp*H));
VH = d`*C*d;
SE = sqrt(VH);
finish seh;

*get heritability for C0 population;
use estmat; read all into v;use covmat; read all into c;
*Note that SAS introduces an extra first column into the C matrix which must be removed;

C = C(|1:nrow(C), 2:ncol(C)|);

*order of variance components in v and c matrices is ENV, SET, REP, ENV*CYCLE,
ENV*SET, ENT(SET) - CYCLE 0, ENT(SET) - CYCLE 1, ENV*ENT(SET) - C0, ENV*ENT(SET) - C1, Error;

*LG and LP vectors for cycle 0;
LG = {0,0,0,0,0,1,0,0,0,0};
LP = {0,0,0,0,0,1,0,1,0,1};

call seh(V,C,LG,LP,H,SE);
print "Heritability on a Plot Basis - Cycle 0", H, SE;
e = 4;
r = 3;
LP = 0//0//0//0//0//1//0//(1/e)//0//(1/(e*r));print LP;
call seh(V,C,LG,LP,H,SE);
print "Heritability on a Family-Mean Basis - Cycle 0", H, SE;

*LG and LP vectors for cycle 1;
LG = {0,0,0,0,0,0,1,0,0,0};
LP = {0,0,0,0,0,0,1,0,1,1};

call seh(V,C,LG,LP,H,SE);
print "Heritability on a Plot Basis - Cycle 1", H, SE;
e = 4;
r = 3;
LP = 0//0//0//0//0//0//1//0//(1/e)//(1/(e*r));
call seh(V,C,LG,LP,H,SE);
print "Heritability on a Family-Mean Basis - Cycle 1", H, SE;
quit; 
run;




