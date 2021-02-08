new;


t=28;       /*time period*/
n=9;       /*cross-section dimension*/
k=2;       /*number of regressors*/ 

load y[t,n] =lnfd.txt;     /*dependent variable*/
load x1[t,n]=lninc.txt;    /*1. regressor*/
load x2[t,n]=lnrem.txt;    /*2. regressor*/

data=x1~x2;



/*do not change here after if you are not familiar with GAUSS*/

x=zeros(t*n,k);
j=1;
do while j<=k;
    if j==1;    x[.,j]=reshape(data[.,j:n]',T*N,1); endif;
    if j>1;     x[.,j]=reshape(data[.,(j-1)*n+1:j*n]',T*N,1); endif;   
    y=reshape(y',T*N,1);
j=j+1;
endo;

i=1;
resid=zeros(t,N);
  do until i>N;
CONSTANT=ones(T,1);
ymat=y[1+(i-1)*T:i*T,.];
xmat=x[1+(i-1)*T:i*T,.]~CONSTANT;
b=ymat/xmat;
e=ymat-xmat*b;
resid[.,i]=e[.,1];
i=i+1;
endo;

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

print "CD Tests                          Stat     prob     ";
print "LM   (Breusch & Pagan 1980)    " cd_LM1;; prob_cdlm1;
print "CDlm (Pesaran 2004)            " cd_LM2;; prob_cdlm2;
print "CD   (Pesaran 2004)            " cd_LM;;  prob_cdlm;


{lma,p_lma}=lmadj(y,x,T,N);

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

print "LMadj                          " lma;;  p_lma;
print "";
print "";
