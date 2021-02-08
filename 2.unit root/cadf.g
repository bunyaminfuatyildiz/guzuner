new;
cls;

t    =28;       //number of time period
n    =9;        //number of countries
ls   =4;        //maxiumum number of lags@
level=0; 	    //level=0 seviye; level=1 birinci fark @
DT   =0;        //DT=0 constant; DT=1 constant+trend @ 
load data[t,n] = lnfd.txt;    


//Do not modify after here
tstart=date;
T=rows(data);
N=cols(data);
y=data;

if level==1;	
y=y[2:T,.]-y[1:T-1,.];
else;
y=y;
endif;


my=meanc(y');       @compute means@

@selection of lags@
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


lagmatt=lagmat';
maxls=maxc(lagmatt);


@compute t statistics@ 
t_stat=zeros(1,N);
resid=zeros(T,N);
{myt1}=varlags(my,1);
mdy=trimr(my,1,0)-myt1;  

{yt1}=varlags(y,1);
dy=trimr(y,1,0)-yt1; 
T=rows(dy);

iter=1;
tt=rows(dy)-maxls;
resid=zeros(tt,N);
  do until iter>N;

dymat=dy[.,iter];     
yt1mat=yt1[.,iter];
ls=lagmat[.,iter];
{dyj}=lagsj(dymat,ls);
{mdyj}=lagsj(mdy,ls);
ymat=trimr(dymat,ls,0);
yt1mat2=trimr(yt1mat,ls,0);
matmdy=trimr(mdy,ls,0);
matmyt1=trimr(myt1,ls,0);
rs1=rows(yt1mat);
rs2=rows(myt1);


xlags=yt1mat2~matmyt1~dyj~matmdy~mdyj;
T=rows(xlags);
CONSTANT=ones(T,1);
      if dt==1;
         xtrend=seqa(1,1,T);
         xmat=xlags~CONSTANT~xtrend;
     else;
     xmat=xlags~CONSTANT;
    endif;
k=cols(xmat);
b=ymat/xmat;
e=ymat-xmat*b;
rss=e'e;        
sigma=rss/(T-k);
var_b=sigma*inv(xmat'xmat);
s_b=sqrt(diag(var_b));
tb=b./s_b;
t_stat[.,iter]=tb[1,1];
ts=maxls-ls;
e2=trimr(e,ts,0);
resid[.,iter]=e2[.,1];
iter=iter+1;

endo;
T=rows(resid);

format 12,3; 
"";
"     Lags     CADF-stat     ";;lagmat'~t_stat';

t_stat=t_stat';
CIPS=meanc(t_stat);
"";
print "CIPS-stat=" CIPS;



proc(1)=varlags(var,lags);
    local yt1;
    yt1=trimr(lag1(var),1,0);
    retp(yt1);
endp;

cx=corrx(resid);

//print cx;
format /rdn 6, 2;
cd=vech(cx);
temp=sumc(cd^2);
@bu ikisinden birisi kullan�labilir@
cd_LM1=T*((cd'cd)-N);
cd_LM12=T*(temp-N);
DF_CDLM1=(N*(N-1))/2;
prob_cdlm1=cdfchic(cd_LM1,DF_CDLM1); 
format /rdn 7,8;

@bu ikisinden birisi kullan�labilir@
cd_LM2=sqrt((1/(N^2-N)))*(CD_LM1-(N^2-N)/2);
cD_LM22=sqrt((1/(N^2-N)))*(sumc(T*cd^2-1)-(N*(T-1)));
prob_cdlm2=1-cdfn(abs(cd_LM2));

cd_LM=(sumc(cd)-N)*sqrt((2*T)/(N^2-N)); 
prob_cdlm=1-cdfn(abs(cd_LM));
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
/*
format 7,4; 
print "Time period=" T;
print "cd Lm1=" cd_LM1;; "  ";; "pvalue=" prob_cdlm1;
print "cd LM2=" cd_LM2;; "  ";; "pvalue=" prob_cdlm2;
print "cd LM="  cd_LM;;  "  ";; "pvalue=" prob_cdlm;
"program is written and checked by Dr. Bulent Guloglu 27 Ocak 2007" ;
//print "correlation matrix of residulas" cx;
*/
