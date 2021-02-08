new;
cls;

t=46;
n=7;
k=2;    /*dependent+explanatory variables*/

load data[t*n,k]=data_cad.txt; 
it=zeros(t*n,2);
for i(1,n,1);
    it[(i-1)*t+1:i*t,1]=zeros(t,1)+i;
    it[(i-1)*t+1:i*t,2]=seqa(1,1,t);
endfor;
 
y=data[.,1]; 
x=data[.,2:cols(data)];

d=1;         /*1=constant; 2=constant and trend*/
p=2;         /*lags in VAR model*/

{b2step,se}= pan2step(it,y,x,d,p);

/*==========================================================
  Computes Panel Two-Step Estimator of Cointegration Vectors             

  Input:               it       NT x 2  matrix. First column: Cross section index   
                                                Second column: Time index               
                        y       NT x 1 vector of dependent variables                   
                        x       NT x k matrix of explanatory variables  
                        d       Deterministic term. d=1: individual specific intercept  
                                d=2: individual specific trend                              
                        p       Lag order of VAR model. p=1: No lagged differences    

  Output:               b2step       k x 1 vector of estimated cointegration coefficients      
                        se      k x 1 vector of standard errors                                      
  ===========================================================          
                        Joerg Breitung 02/02/05         
  =========================================================== */
proc(2)=pan2step(it,y,x,d,p);
local NT,k,xx,xy,yy,N,NTges,pos,z,ind,yi,Ti,ct,yt,xt,b,e,e1,dx,dy,dz,
        z1,lags,j,dztilde,z1tilde,e1tilde,ages,a,u,S,Si,gam,dzstar,b2step,sigu,se;

  NT=rows(x);
  k=cols(x);  
 
  p=p-1;
  xx=0;
  xy=0;
  yy=0;
  N=0;
  NTges=0;
  pos=1;
  z=y~x;
  z=z|(zeros(1,k+1));
  it=it|zeros(1,cols(it));

  do while pos < NT-1;

    ind=it[pos,1];
    yi=z[pos,1:k+1];
    pos=pos+1;
    do while it[pos,1]==ind;
      yi=yi|z[pos,1:k+1];
      pos=pos+1;
    endo;
    Ti=rows(yi);
    y=yi[.,1];
    x=yi[.,2:k+1];                            
    ct=ones(Ti,1);
    if d==2;
      ct=ct~seqa(1,1,Ti);
    endif;
    if d>0;
      yt=y-ct*invpd(ct'ct)*ct'y;
      xt=x-ct*inv(ct'ct)*ct'x;
    endif;

    b=inv(xt'xt)*xt'yt;
  
    e=y-x*b;
    e1=e[1:Ti-1];
    dx=x[2:Ti,.]-x[1:Ti-1,.];
    dy=y[2:Ti]-y[1:Ti-1];
    dz=dx~dy;
    z1=x[1:Ti-1,.]~y[1:Ti-1];

    lags=ct[p+2:Ti,.];

    if p>0;
      lags=lags~dz[p:Ti-2,.];
      j=2;
      do while j<=p;
        lags=lags~dz[p+1-j:Ti-j-1,.];
      j=j+1;
      endo;
      dx=dx[p+1:Ti-1,.];
      dy=dy[p+1:Ti-1,.];
      dz=dz[p+1:Ti-1,.];
      e1=e1[p+1:Ti-1,.];
      z1=z1[p+1:Ti-1,.];
   endif;
      dztilde=dz-lags*invpd(lags'lags)*lags'dz;               
      z1tilde=z1-lags*invpd(lags'lags)*lags'z1;
      e1tilde=e1-lags*invpd(lags'lags)*lags'e1;
 
    ages=inv(e1tilde'e1tilde)*e1tilde'dztilde;
    a=ages[1,.]';
                                                                   
    u=dztilde-e1tilde*ages;
    S=u'u/(Ti-1-p);
    Si=inv(S);
    gam=Si*a;
    dzstar=z1tilde[.,k+1]-dztilde*gam*inv(a'gam);
    z1tilde=z1tilde[.,1:k];
    xy=xy+dzstar'z1tilde;
    xx=xx+z1tilde'z1tilde;
    yy=yy+dzstar'dzstar;

    N=N+1;
    NTges=NTges+(Ti-p-1);
  endo;
 
  b2step=inv(xx)*xy';
  sigu=yy-b2step'*xy';
  sigu=sigu/NTges;
  se=sqrt(sigu*diag(inv(xx[1:k,1:k])));

format 8,4;
"beta     std.error   t-stat";b2step~se~b2step./se;

retp (b2step[1:k],se);

endp;


@ ============== end of Programm ================ @

