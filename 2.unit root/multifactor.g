new;
cls;

t=28;      //number of time period
n=9;       //number of countries

load ymat[t,n]  = lnfd.txt;
load x1[t,n]    = lninc.txt; /*factor1*/
load x2[t,n]    = lnrem.txt; /*factor2*/
xmat            = x1|x2;
xbar_mat        = meanc((x1~x2)');
case            = 2;         /* case=2: with intercept*/
                             /* case=3: with intercept and trend*/                           
maxp            = 3;         /* maksimum number of lags*/       

{cips_p_mat,cadf_p_mat,res_p_mat}=cipsm(ymat,xmat,maxp,case);
{CSB_p_mat,CSBi_p_mat}=CSB(ymat,xbar_mat,maxp,case);

format 8,3; 
print "    lag     CIPS      CSB";;seqa(0,1,maxp+1)~cips_p_mat[.,1]~CSB_p_mat;
"The lag  in Pesaran et al (2013)= ";;int(4*(T/100)^(1/4));



_dxmiss=1;
proc(3)=cipsm(y_mato,x_mato,maxp,case);
local outpt,n,t,tlag,y_mat,x_mat,y_mat_1,Dy_mat_temp,Dy_mat,Dvar_CSM,var_1_CSM;
local k, c_k,x_mato_k,xbar_mat,xbar_mat_1,Dxbar_mat_temp,Dxbar_mat;
local tt,y,z,hi,x,h,k1,k2,bvec,sevec,t_vec_cce,t_bar,i,c_p,temp;
local trncl,truncu,trnc1,cadf_p,cadf_ps,cips_p_mat,cips_p,cips_ps;
local cadf_p_mat,cadfs_p_mat,idvec,lgorder,res_p_mat,res_p,ckmat;


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
outpt=1;/* 1: reports ourput, 0:supress the output */
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
outwidth 256;

t=rows(y_mato);
n=cols(y_mato);



if x_mato==0;K=0;
else;
K=rows(x_mato)/t;
endif;


/** diviation from the initial cross section means **/
y_mat=y_mato-meanc(y_mato[1,.]');
x_mat={};
c_k=1;
do while c_k<=k;
  x_mato_k=x_mato[1+(c_k-1)*T:c_k*T,.];
  x_mat=x_mat|(x_mato_k-meanc(x_mato_k[1,.]'));
c_k=c_k+1;
endo;

tlag=maxp+1;

y_mat_1=lag(y_mat);
Dy_mat_temp=y_mat-y_mat_1;
    y_mat_1=y_mat_1[tlag+1:t,.];
    Dy_mat=Dy_mat_temp[tlag+1:t,.];

/* Txk matrix of Xbar */
if k>0;
xbar_mat=reshape(meanc(x_mat'),K,T)';
xbar_mat_1=lag(xbar_mat);
Dxbar_mat_temp=xbar_mat-xbar_mat_1;
    xbar_mat_1=xbar_mat_1[tlag+1:t,.];
    Dxbar_mat=Dxbar_mat_temp[tlag+1:t,.];

Dvar_CSM=meanc(Dy_mat')~Dxbar_mat;
var_1_CSM=meanc(y_mat_1')~xbar_mat_1;

elseif k==0;
Dvar_CSM=meanc(Dy_mat');
var_1_CSM=meanc(y_mat_1');

endif;


tt=t-tlag;

/*************/
y=vec(Dy_mat);
Hi=Dvar_CSM~var_1_CSM;
x=vec(y_mat_1);

/** Choose K1 and K2 (truncation threshold) */
/** NOT USED IN THIS VERSION!!   **/

ckmat={
6.12	4.16	6.19	2.61	6.42	1.70,
6.53	3.75	6.65	2.57	6.90	1.82,
6.78	3.34	6.86	2.36	7.14	1.73,
7.04	3.10	7.20	2.30	7.48	1.74,
7.28	2.87	7.41	2.16	7.63	1.61,
7.52	2.67	7.66	2.04	7.86	1.55,
7.59	2.42	7.75	1.87	7.97	1.43,
7.83	2.28	7.99	1.78	8.26	1.42};

if case==1;k1=-ckmat[k+1,1];k2=ckmat[k+1,2];
elseif case==2;Z=ones(n*tt,1);k1=-ckmat[k+1,3];k2=ckmat[k+1,4];
elseif case==3;Z=ones(n*tt,1)~(ONES(N,1).*.seqa(1,1,tt));k1=-ckmat[k+1,5];k2=ckmat[k+1,6];
else;"please choose the 'case' only from 1,2, or 3";end;
endif;

/**************** ESTIMATION ****************/
c_p=0;
cips_p_mat={};
cadf_p_mat={};
cadfs_p_mat={};
res_p_mat={};
do while c_p<=maxp;

if c_p>0;
    temp=(Dy_mat_temp[tlag+1-c_p:t-c_p,.]);

        if k>0;    
            Hi=Hi~meanc(temp')~Dxbar_mat_temp[tlag+1-c_p:t-c_p,.];
        elseif k==0;
            Hi=Hi~meanc(temp'); 
        endif;

    x=x~vec(temp);

endif;
/*************/
if case==1;H=(ONES(N,1).*.Hi);
else;H=Z~(ONES(N,1).*.Hi);
endif;

/*** MGCCE ****/
{bvec,sevec,t_vec_cce,res_p}=mgsimple_cipsm(y,x~h,n,tt);

trncl=k1*(t_vec_cce[.,1].<k1);truncu=k2*(t_vec_cce[.,1].>k2);
trnc1=(t_vec_cce[.,1].>k1).*(t_vec_cce[.,1].<k2);

cadf_p = t_vec_cce[.,1];
cadf_ps = (t_vec_cce[.,1].*trnc1+trncl+truncu);

cips_p=meanc(cadf_p);
cips_ps=meanc(cadf_ps);

if outpt==1;
format /rd 2,0;
/*
"@@@@@@@@@@@@@@@@@@@ CADF( p=";;C_P;;") @@@@@@@@@@@@@@@@@@@@";
format /rd 2,0;
"CIPSM TEST: CASE";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT AND TREND";
endif;
format /rd 5,0;
"N=";;n;;", T=";;tt;;", (";;t;;"data points used)";
"Number of X, k=";;k;
format /rd 8,3;
"   CIPSM            :";;cips_p;
"   CIPSM*(truncated):";;cips_ps;
"* Truncation is done for CADFM_i in such a way that" ;
format /rd 1,2;
"when CADF_i<k1, CADF_i=k1 and when CADF_i>k2, CADF_i=k2,";
"where ";;"k1=";;k1;;" and k2=";;k2;;".";
"Appropriate critical values are in Pesaran et al (2007).";
"-------------------------------------------------------";*/
"";
endif;
cips_p_mat=cips_p_mat|(cips_p~cips_ps);
cadf_p_mat=cadf_p_mat~cadf_p;
cadfs_p_mat=cadfs_p_mat~cadf_ps;

res_p_mat=res_p_mat|res_p;

    c_p=c_p+1;
endo;

if outpt==1;
format /rd 8,0;
idvec=ftocv(seqa(1,1,n),1,0);
lgorder="id/p"~ftocv(seqa(0,1,maxp+1)',1,0);
/*
format /rd 2,0;
"";
"*****************************************************";
"CADF_i(p) Statistics: Case";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT AND TREND";
endif;
format /rd 8,0;
$ (lgorder|(idvec~ftocv(cadf_p_mat,1,3)));
"";

"*****************************************************";
format /rd 2,0;
"CADF*_i(p) Statistics (truncated): Case";;case;;
if case==1;", NO INTERCEPT";
elseif case==2;", WITH INTERCEPT";
elseif case==3;", WITH INTERCEPT and TREND";
endif;
format /rd 8,0;
$ (lgorder|(idvec~ftocv(cadfs_p_mat,1,3)));
"* Truncation is done for CADF_i in such a way that" ;
format /rd 1,2;
"when CADF_i<k1, CADF_i=k1 and when CADF_i>k2, CADF_i=k2,";
"where ";;"k1=";;k1;;" and k2=";;k2;;".";
"Appropriate critical values are in Table 1, Pesaran (2006).";*/

endif;

format /rd 16,8;
retp(cips_p_mat,cadf_p_mat,res_p_mat);
endp;

/*** OLS regression for each cross section unit ****/
proc(4)=mgsimple_cipsm(y,x,n,t);
local bvec,sevec,tvec,i,y_i,x_i,beta_i,e_i,zig2,se_i,t_i,res_mat;
bvec=zeros(n,cols(x));
sevec=zeros(n,cols(x));
tvec=zeros(n,cols(x));
res_mat={};

if t<cols(x);
bvec=_dxmiss*ones(n,cols(x));
sevec=_dxmiss*ones(n,cols(x));
tvec=_dxmiss*ones(n,cols(x));
res_mat=_dxmiss*ones(t,n);
goto tln;
endif;

i=1;
do while i<=n;
    y_i=y[1+(i-1)*t:(i)*t,.];
    x_i=x[1+(i-1)*t:(i)*t,.];
    beta_i=pinv(x_i'x_i)*x_i'y_i;
    e_i=y_i-x_i*beta_i;
    zig2=e_i'e_i/(T-cols(x));
    se_i=sqrt(diag(zig2*pinv(x_i'x_i)));
    t_i=beta_i./se_i;

    bvec[i,.]=beta_i';
    sevec[i,.]=se_i';
    tvec[i,.]=t_i';
    res_mat=res_mat~e_i;

i=i+1;
endo;

tln:
retp(bvec,sevec,tvec,res_mat);
endp;

proc(2)=CSB(y_mat,xbar_mat,maxp,case);
local i,n,y,y_1,dy,c_p,dyp,s,s2,tstar,x,xx,phi,e,se,SBstat,xr,phir,er,
CSB_p_mat,CSBi_p_mat,sr2,xbar,zbar, dzbar, dzbarp,k;

n=cols(y_mat);
if xbar_mat==0;k=0;
zbar=meanc(y_mat');
zbar=zbar-zbar[1,.];
dzbar=zbar-lag(zbar);
dzbar=trim(dzbar,1,0);
dzbarp=plag(dzbar,maxp);

else;

xbar=xbar_mat;
k=cols(xbar);
zbar=meanc(y_mat')~xbar;
zbar=zbar-zbar[1,.];
dzbar=zbar-lag(zbar);
dzbar=trim(dzbar,1,0);
dzbarp=plag(dzbar,maxp);

endif;

CSBi_p_mat={};

i=1;
do while i<=n;
y=y_mat[.,i]-y_mat[1,i];/** deviation from the initial values **/

y_1=lag(y);
dy=y-y_1;

c_p=1;dyp=dy;
do while c_p<=maxp;
if c_p==1;dyp=lag(dyp);
else; dyp=dyp~lag(dyp[.,c_p-1]);
endif;
c_p=c_p+1;
endo;
dyp=dyp[maxp+2:rows(dy),.];


c_p=0;
do while c_p<=maxp;

s=dy[maxp+2:rows(dy)];
tstar=rows(s);/**effective # of obs. */

if c_p==0;

    if case==2;

        xr=dzbarp[.,1:(k+1)];
            if tstar<=cols(xr);SBstat=SBstat|_dxmiss;goto label1;
            else;
            phir=s/xr;er=s-xr*phir;sr2=er'er/(tstar-cols(xr));
            endif;

    elseif case==3;

        
        xr=ones(tstar,1)~dzbarp[.,1:(k+1)];
           if tstar<=cols(xr);SBstat=SBstat|_dxmiss;goto label1;
            else;phir=s/xr;er=s-xr*phir;sr2=er'er/(tstar-cols(xr));
		   endif;
	

    elseif (case/=2 or case/=3);
	print "ERROR MESSAGE: Please use case 2 or case 3 only";
    end;	
		endif;

SBstat=sumc(cumsumc(er)^2)/(sr2*rows(er)^2);

label1:

else;

xr=xr~dyp[.,c_p]~dzbarp[.,1+c_p*(k+1):(c_p+1)*(k+1)];

if tstar<=cols(xr);
SBstat=SBstat|_dxmiss;
else;
phir=s/xr;er=s-xr*phir;sr2=er'er/(tstar-cols(xr));
SBstat=SBstat|(sumc(cumsumc(er)^2)/(sr2*rows(er)^2));
endif;

endif;

c_p=c_p+1;
endo;

CSBi_p_mat=CSBi_p_mat|SBstat';
i=i+1;
endo;

CSB_p_mat=meanc(CSBi_p_mat);

retp(CSB_p_mat,CSBi_p_mat);
endp;




proc(1)=plag(y,p);
local yp,c_p;
yp=y[p+1:rows(y),.];
c_p=1;
do while c_p<=p;
yp=yp~y[p+1-c_p:rows(y)-c_p,.];
c_p=c_p+1;
endo;
retp(yp);
endp;


