/* The original version: October 25, 2010 by Eiji Kurozumi.*/

 /* The original version is modified to select optimal lags
  for each individuals as well as for unit root tests for both level and 1.difference
  by Saban Nazlioglu January 19, 2013.  */

/* -----------------------------------------------------------------------
**      Obtain the Panel Stationarity Tests with Cross-Sectional Dependence
**      Proposed by Hadri and Kurozumi (2008), "A Simple Panel Stationarity
**      Tests in the Presence of Cross-Sectional Dependence"
**
**      format   {ZA_spc,ZA_la}=panelst_HK(y,const,p_all);
**
**      input   y:      Panel Data (N by T)
**              const:  const=1 if only a constant is included
**                      const=2 if a constant and a linear trend is included
**              p_all:  N by 1 vector containing lag orders for individuals
**                      if p_all is scalar, then the same lag order (p_all) is
**                      used for all the individuals
**
**      output  ZA_spc: the augmented panel KPSS test statistic with long-run
**                      variance correced by the SPC method
**              ZA_la:  the augmented panel KPSS test statistic with lonsg-run
**                      variance correced by the LA method
** -----------------------------------------------------------------------*/
new;
cls;

t    =28;                  // time dimension
n    =9;                   // cross-section dimension
ls   =4;                   // maximum lags
level=0; 		           // level=0 seviye, level=1 birinci fark
const=1;                   // const=1 constant; const=2 constant and trend

load data[t,n] = lnfd.txt;        // Panel Data (T by N)


if level==1; y=data[2:T,.]-data[1:T-1,.];
else;       y=data;
endif;

//selection of lags (p_all is automatically selected by Schwarz information criterion
iter=1;
{yt1}=varlags(y,1);
dy=trimr(y,1,0)-yt1;  
N=cols(dy);
lagmat=zeros(1,N);
{dyj}=lagsj(dy,ls);
yt1=trimr(yt1,ls,0);
dy=trimr(dy,ls,0);

do until iter>N;
ymat=dy[.,iter];     
yt1mat=yt1[.,iter];   

p=1;
SICMAT=zeros(ls,1);
   do until p>ls;
xlags=dyj[.,1:p];
T=rows(dy);
Constant=ONES(T,1);
dt=const;
if dt==1;
 xtrend=seqa(1,1,T);
 xmat=yt1mat~Constant~xtrend~xlags;
else;
xmat=constant~yt1mat~xlags;
endif;
k=cols(xmat);

b=ymat/xmat;
e=ymat-xmat*b;
rss=e'e;   
SIC=T*ln(RSS)+k*ln(T);
SICMAT[p,.]= SIC;
if p==ls;
  maklag=minindc(SICMAT);
  lagmat[.,iter]=maklag;
endif;
p=p+1;
    endo;
iter=iter+1;
endo;
proc(1)=varlags(var,lags);
    local yt1;
    yt1=trimr(lag1(var),1,0);
    retp(yt1);
endp;

proc(1)=lagsj(varx,j);
    local N, T, i,l,cs,p,m,dyj;
    N=cols(varx)*j;
    T=rows(varx)-j;
    dyj=zeros(T,N);
    cs=cols(varx);
    l=1;
    p=1; m=cs;
    do until m>N;
    dyj[.,l:m]=trimr(lagn(varx,p),j,0);
    p=p+1;
    l=l+cs ;
    m=m+cs;
    endo;
  retp(dyj);
endp;


p_all=lagmat';

y=y';                               //do not change here         // Panel Data (N by T)

{ZA_spc,ZA_la,ybar}=panelst_HK(y,const,p_all);


format /m1/rd 8,4;
print "ZA_spac and p-val =";;ZA_spc;; "[";;cdfnc(za_spc);;"]";
print "ZA_la   and p-val =";;ZA_la;;; "[";;cdfnc(za_la);;"]";
"";
//print "Optimal Lags for Individuals";;p_all;

proc (3)=panelst_HK(y,const,p_all);

local N_dim,T_dim,yn,ybar,sqrn,meanv2,sigmav2,z,p;
local yx_all,yn1,mx,epsilon1,pars,num,lrv1,lrv2,j1,LM1,LM1bar,LM2,LM2bar;
local ZA_la,ZA_spc;

N_dim=rows(y);
T_dim=cols(y);

if rows(p_all) == 1;
    p_all=p_all+zeros(N_dim,1);
endif;

yn=y';
ybar = sumc(y);

if const == 1;
    sqrn=sqrt(N_dim); meanv2=1/6;  sigmav2=sqrt(1/45);
    z=ones(T_dim,1);
elseif const == 2;
    sqrn=sqrt(N_dim); meanv2=1/15;  sigmav2=sqrt(11/6300);
    z=ones(T_dim,1)~seqa(1,1,T_dim);
endif;

LM1=0;
LM2=0;
j1=1;
do until j1 > N_dim;
        p        = p_all[j1,1];
        yx_all   = packr(yn[.,j1]~z~lagn(ybar,seqa(0,1,p+1)) );
        yn1      = yx_all[.,1];
        mx       = yx_all[.,2:cols(yx_all)];
        epsilon1 = yn1-mx*invpd(mx'mx)*mx'*yn1;
        pars     = cumsumc(epsilon1);	/*partial sum of residuals*/
        num      = sumc(pars.*pars);
		lrv1     = lrv_Sul(yn[.,j1],ybar,const,p);
        lrv2     = lrv_LA(yn[.,j1],ybar,const,p);
        LM1      = LM1+T_dim^(-2)*num./lrv1;
        LM2      = LM2+T_dim^(-2)*num./lrv2;
        j1=j1+1;
endo;

LM1bar=LM1/N_dim;
LM2bar=LM2/N_dim;

ZA_spc  = sqrn*(LM1bar-meanv2)/sigmav2;
ZA_la   = sqrn*(LM2bar-meanv2)/sigmav2;

retp(ZA_spc,ZA_la,ybar);
endp;




/* -----------------------------------------------------------------------
**      calculate the long-run variance using Sul, Phillips and Choi's method
**      augmented by ybar
**
**
**      format   {lrv1}=lrv_Sul(y,ybar,const,p);
**
**      input   y:      time series (T by 1)
**              ybar:   cross sectional average of y
**              const:  const=1 if only a constant is included
**                      const=2 if a constant and a linear trend is included
**              p:      the order of the AR part
**
**      output  lrv1:   the long-run variance of y
** -----------------------------------------------------------------------*/

proc(1)=lrv_Sul(y,ybar,const,p);

local y_all,y_0,y_l,z,et,var_e,coef_all,coef_ar,phi1,bd1,lrv1;
@@
y_all=packr( (lagn(y,seqa(0,1,p+1 )))~(lagn(ybar,seqa(0,1,p+1 ))) );
y_0=y_all[.,1];
y_l=y_all[.,2:cols(y_all)];
if const == 1;
    z=y_l~ones(rows(y_0),1);
elseif const == 2;
    z=y_l~ones(rows(y_0),1)~seqa(1,1,rows(y_0));
endif;
coef_all=invpd(z'z)*z'y_0;
et=y_0-z*coef_all;
var_e=meanc(et.^2);
coef_ar =coef_all[1:p,1];
phi1=sumc(coef_ar);
bd1=1-1/sqrt(rows(y_0));
phi1=minc(phi1|bd1);
lrv1=var_e/(1-phi1)^2;
@@
retp(lrv1);
endp;





/* -----------------------------------------------------------------------
**      calculate the long-run variance using the lag augmenting method
**      augmented by ybar and its lags
**
**      y(t)=const+a(1)*y(t-1)+...+a(p)*y(t-p)+e(t)
**
**      lrv=var(e)/(1-a(1)-...-a(p))^2
**
**
**
**      format   {lrv1}=lrv_LA(y,const,p);
**
**      input   y:      time series (T by 1)
**              ybar:   cross sectional average of y
**              const:  const=1 if only a constant is included
**                      const=2 if a constant and a linear trend is included
**              p:      the order of the AR part
**
**      output  lrv1:   the long-run variance of y
** -----------------------------------------------------------------------*/

proc(1)=lrv_LA(y,ybar,const,p);

local y_all,y_0,y_l,z,et,var_e,coef_all,coef_ar,lrv1;
@@
@ estimation of Var(e(t)) @
@@
y_all=packr( (lagn(y,seqa(0,1,p+1 )))~(lagn(ybar,seqa(0,1,p+1 ))) );
y_0=y_all[.,1];
y_l=y_all[.,2:cols(y_all)];
if const == 1;
    z=y_l~ones(rows(y_0),1);
elseif const == 2;
    z=y_l~ones(rows(y_0),1)~seqa(1,1,rows(y_0));
endif;
et=y_0-z*invpd(z'z)*z'y_0;
var_e=meanc(et.^2);
@@
@ estimation of the AR coefficeints @
@@
y_all=packr(lagn(y,seqa(0,1,p+2 ))~lagn(ybar,seqa(0,1,p+1 )));
y_0=y_all[.,1];
y_l=y_all[.,2:cols(y_all)];
if const == 1;
    z=y_l~ones(rows(y_0),1);
elseif const == 2;
    z=y_l~ones(rows(y_0),1)~seqa(1,1,rows(y_0));
endif;
coef_all=invpd(z'z)*z'y_0;
coef_ar =coef_all[1:p,1];

lrv1=var_e/(1-sumc(coef_ar))^2;
@@
retp(lrv1);
endp;
