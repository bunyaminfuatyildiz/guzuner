new;


t=28;   //time period
n=9;   //cross-section dimension

load y[t,n]=lnfd.txt;

ls=1;                   @number of lags@
DT=0;   		        @DT=0 constant, DT=1 constant and trend@ 
level=0; 		        @level=0 duzey, level=1 birinci fark@
T=rows(y);
N=cols(y);



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

//print "number of lags for each country " lagmat;
lagmatt=lagmat';
maxls=maxc(lagmatt);
//print "maximum lag" maxls;

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

//print "dymat" dymat;
//print "CADF_statistics for each country";
//print t_stat;
/*
t_stat=t_stat';
CIPS=meanc(t_stat);
print "CIPS statistics";;CIPS;
*/
cx=corrx(resid);
cd=vech(cx);
temp=sumc(cd^2);

cd_LM1=T*((cd'cd)-N);                           //Breusch-Pagan (1980) LM test
DF_CDLM1=(N*(N-1))/2;
prob_cdlm1=cdfchic(cd_LM1,DF_CDLM1); 
format /rdn 7,3;

cd_LM2=sqrt((1/(N^2-N)))*(CD_LM1-(N^2-N)/2);    //Pesaran (2004) CDlm test
prob_cdlm2=1-cdfn(abs(cd_LM2));

cd_LM=(sumc(cd)-N)*sqrt((2*T)/(N^2-N));         //Pesaran (2004) CD test
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
proc(1)=varlags(var,lags);
    local yt1;
    yt1=trimr(lag1(var),1,0);
    retp(yt1);
endp;


/** Compute bias-adjusted LM test statistic proposed by
Pesaran, M.H., Ullah, A., and Yamagata, T.,
"A bias-adjusted LM test for error cross-section independence",
 Econometrics Journal (2008), vol 11 pp. 105-127.
 
 This is for balanced panel model.

y_it = x_it'beta_i + u_it, i=1,2,...,N;t=1,2,...,T
where x_it is a kx1 exogenous regressor vector, including deterministics (e.g. intercept),
u_it~IIDN(0,sigma_i^2)

The hull hypothesis to be tested is
H0: E(u_it,u_jt)=0 for all i,j,t.

{lma,p_lma}=lmadj(y,X,T,N);

INPUT: 
y = (y_1',y_2',...,y_N')' with y_i=(y_i1,y_i2,...,y_iT)';: a TNx1 vector of dependent variable
X = (X_1',X_2',...,X_N')' with X_i=(x_i1,x_i2,...,x_iT)': a TNxk vector of explanatory variable
T: time series dimension (to be the same for all cross section units)
N: cross section dimension

OUTPUT:
lma: bias-adjuested LM test statistic
p_lma: p-value of lma

Note that under the null hypothesis, lma~N(0,1) approximately for large T and N.
The test is upper-one-tailed test.

TY Sept 2011

*/

y=reshape(y,T*N,1);
//print "y"; y;
x=xmat;
x=reshape(x,T*N,k);
//print "x"; x;

{lma,p_lma}=lmadj(y,x,T,N);
//$"LMba"~"p-value";
//lma;;"[";;p_lma;;"]";



proc(2)=lmadj(y,x,T,N);
local k,mmat,emat,x_i,px_i,stat_ij,i,j,e_i,v_i,e_j,v_j,cov_ij,m_j,m_i,rho_ij
,mu_ij,a1_ij,a2_ij,vsq_ij,sij,lma;

k=cols(x);
mmat={};
emat={};
i=1;
do while i<=n;
    x_i=x[1+(i-1)*T:i*T,.];
    px_i=x_i*inv(x_i'x_i)*x_i';
    mmat=mmat~(eye(t)-px_i);
    emat=emat~((eye(t)-px_i)*y[1+(i-1)*T:i*T]);
    i=i+1;
endo;

stat_ij=0;
i=1;
do while i<=n-1;
m_i=mmat[.,1+(i-1)*t:i*t];
    e_i=(emat[.,i]-meanc(emat[.,i]));
    v_i=sumc(e_i^2);
j=i+1;
do while j<=n;
m_j=mmat[.,1+(j-1)*t:j*t];

    e_j=(emat[.,j]-meanc(emat[.,j]));
    v_j=sumc(e_j^2);
    cov_ij=(e_i'e_j);
    rho_ij=cov_ij/(sqrt(v_i)*sqrt(v_j));
    sij=sumc(diag(m_i*m_j));
    mu_ij= sij/(t-k);
    a2_ij=3*(((t-k-8)*(t-k+2)+24)/((t-k+2)*(t-k-2)*(t-k-4)))^2;
    a1_ij=a2_ij-1/(t-k)^2;

    vsq_ij=(sij^2)*a1_ij + 2*sumc(diag(m_i*m_j*m_i*m_j))*a2_ij;

    stat_ij=stat_ij+(((t-k)*(rho_ij^2) - mu_ij)/sqrt(vsq_ij));
j=j+1;
endo;
i=i+1;
endo;

lma=sqrt(2/(n*(n-1)))*(stat_ij);
p_lma=cdfnc(lma);
retp(lma,p_lma);
endp;



/*printing results*/
print "CD Tests                        Stat     prob     ";
print "LM    (Breusch,Pagan 1980)    " cd_LM1;; prob_cdlm1;
print "CDlm  (Pesaran 2004)          " cd_LM2;; prob_cdlm2;
print "CD    (Pesaran 2004)          " cd_LM;;  prob_cdlm;
print "LMadj (PUY, 2008)             " lma;;  p_lma;
"";
"Lags";
lagmat';
