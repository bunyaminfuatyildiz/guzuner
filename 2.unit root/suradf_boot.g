new;cls;

t       =28;       //number of time period
n       =9;        //number of countries
ls      =4;        // maximum lags
maxboot =1000;     //bootstrap number 10000
level   =0;        @ level=0 seviye; 
                     level=1 birinci fark@
mod     =1;        @ mod=0 constant; 
                     mod=1 constant+trend @

load data[t,n] =lninc.txt;

//do not modify the following lines//
y=data;
if level==1; y=y[2:T,.]-y[1:T-1,.];
else; y=y;  endif;

//selection of lags by Schwarz

{yt1}=varlags(y,1);
dy=trimr(y,1,0)-yt1;  
N=cols(dy);
lagmat=zeros(1,N);
{dyj}=lagsj(dy,ls);
yt1=trimr(yt1,ls,0);
dy=trimr(dy,ls,0);

i=1;
do until i>N;
ymat=dy[.,i];     
yt1mat=yt1[.,i];   

p=1;
SICMAT=zeros(ls,1);
   do until p>ls;
xlags=dyj[.,1:p];
T=rows(dy);
Constant=ONES(T,1);
//dt=constant;
if mod==1; xmat=yt1mat~Constant~seqa(1,1,T)~xlags;endif;
if mod==0; xmat=constant~yt1mat~xlags; endif;
k=cols(xmat);

b=ymat/xmat;
e=ymat-xmat*b;
rss=e'e;   
SIC=T*ln(RSS)+k*ln(T);
SICMAT[p,.]= SIC;
if p==ls;
   maklag=minindc(SICMAT);
   lagmat[.,i]=maklag;
endif;
p=p+1;
    endo;
i=i+1;
endo;
selectmat=lagmat';


N=cols(data);   @ N is the number of countries @
g=data;
dg=g[2:rows(g),.]-g[1:rows(g)-1,.];
g1=g[1:rows(g)-1,.];
g=g[2:rows(g),.];
T=rows(g);

p=maxc(selectmat);
tef=T-p;

@print "dg=" dg;@
@print "g1=" g1;@
@print "T " T; @ 

cs=sumc(selectmat);
{beta,vbeta,varresid}=panelSURADF(dg,g1,p);

suradf_t=beta[N+1:2*N,1]./sqrt(vbeta[N+1:2*N,1]);
suradf_tboot=zeros(maxboot,N);
b_dg=zeros(N,p);
beta=beta';
itc=1;
do while itc<=N;
gec=sumc(selectmat[1:itc]); 

b_dg[itc,1:selectmat[itc]]=beta[1,2*N+gec-selectmat[itc]+1:2*N+gec];
  itc=itc+1;
endo;
ben=beta';
level_lags=ones(1,N)+selectmat';
Tint=T-1;
@"compute critical values"@
boot=1;
do while boot<=maxboot;


ut1=rndn(T+50,N+50);
ut2=ut1[51:T+50,51:N+50];
et=chol(varresid)*ut2';
et=et';

dggen=zeros(T,N);
   ggen=zeros(T+1,N);
   ggen[1,.]=g[1,.];
   dggen[1:p,.]=dg[1:p,.];
   itc=1;
  do while itc<=N;
      gec=sumc(selectmat[1:itc]); 
      itf=p+1;
      do while itf<=T;
         dggen1=zeros(rows(ben),1);
         dggen1[itc]=0;
         dggen1[N+itc]=0;
         dggen1[2*N+gec-selectmat[itc]+1:2*N+gec]=
             rev(dggen[itf-p:itf-p+selectmat[itc]-1,itc]);   
            dggen[itf,itc]=dggen1'ben+et[itf-p,itc];
         itf=itf+1;
      endo;
      it=2;      
    do while it<=T+1;
         ggen[it,itc]=dggen[it-1,itc]+ggen[it-1,itc];
         it=it+1;
      endo;
      itc=itc+1;
   endo;
   g1gen=ggen[1:T,.];
 
{betagen,vbetagen,varresidgen}=panelSURADF(dggen,g1gen,p);

suradf_tboot[boot,.]=(betagen[N+1:2*N,1]./sqrt(vbetagen[N+1:2*N,1]))';
boot=boot+1;

endo;

suradf_cv=zeros(3,N); 
i=1;
 do while i<=N;
 t1=suradf_tboot[.,i];

t_ro1sort=sortc(t1,1);  
kant10=0.1;
kant5 =0.05;
kant1 =0.01;

suradf_cv[1,i]=quantile(t1,kant10);
suradf_cv[2,i]=quantile(t1,kant5);
suradf_cv[3,i]=quantile(t1,kant1);
i=i+1;
endo;

      
proc(3)=panelSURADF(dzm,z1m,pmax); 
   local rdzm,y,x,gec,unos,x1,x2,trend,itc,itp,be,en,em,
         vari1,vari2,ivari1,ivari2,bes,ens,ems,vari1s,xxi;
   rdzm=rows(dzm);
   y=zeros(tef*N,1);
   unos=zeros(tef*N,N);
   x1=zeros(tef*N,N);   
   x2=zeros(tef*N,cs);
   trend=zeros(tef*N,N); 
   itc=1;
   do while itc<=N;
      y[(itc-1)*tef+1:itc*tef]=dzm[pmax+1:rdzm,itc];
      unos[(itc-1)*tef+1:itc*tef,itc]=ones(tef,1);
      x1[(itc-1)*tef+1:itc*tef,itc]=z1m[pmax+1:rdzm,itc];
      trend[(itc-1)*tef+1:itc*tef,itc]=seqa(1,1,tef);
      itp=1;
      do while itp<=selectmat[itc];
          gec=sumc(selectmat[1:itc]); 
          x2[(itc-1)*tef+1:itc*tef,(gec-selectmat[itc]+itp)]=dzm[p-itp+1:T-itp,itc];
           
  itp=itp+1;
    endo;
      itc=itc+1;
    endo;
  if mod==1; x=unos~x1~x2~trend;endif;
  if mod==0; x=unos~x1~x2;endif;         
  
   be=inv(x'x)*x'y;
   en=y-x*be;
   em=en[1:tef];
   itc=2;
   do while itc<=N;
      em=em~en[(itc-1)*tef+1:itc*tef];
      itc=itc+1;
   endo;
   vari1=em'em/tef;
   ivari1=inv(vari1);
   vari2=vari1.*.eye(tef);
   ivari2=ivari1.*.eye(tef);
   xxi=inv(x'ivari2*x);
   bes=xxi*x'ivari2*y;
   ens=y-x*bes;
   ems=ens[1:tef];
   itc=2;
   do while itc<=N;
      ems=ems~ens[(itc-1)*tef+1:itc*tef];
      itc=itc+1;
   endo;
   vari1s=ems'ems/tef;
   retp(bes,diag(xxi),vari1s);
endp;


 format 8,4;
"    id     Lags    SURADF     10%     5%    1%";;seqa(1,1,n)~selectmat~suradf_t~suradf_cv';


//  " number of periods  " tef;
//  " maximum lag for level" p+1;
 // "   level lags " level_lags;
 // " birler matrisi" ones(1,N);
 // " sifirlar matrisi" zeros(1,N);

 // " coefficients on lagged differences"  b_dg; 
 // " SUR variance covariance matrisi" varresid;


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
