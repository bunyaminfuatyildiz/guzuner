
new;
cls;
  format /ld 12,3;
  library coint pgraph; graphset;
_olsres  = 1;
 __output= 0;
  __miss = 1; 
rndseed 83478;

d1 = date;

/* Set-up for parameters */
  t   = 28;    @number of time period@
  n   = 9;     @number of countries@
  nb  = 5000;  @Bootstrap Replications@
  pb  = 11;     @Maximum lag length @
  m   = 100;   @Block size @
  d   = 0;     @0=constant only, 1=constant and time trend@ 
  level=0; 	   @level=0 seviye; level=1 birinci fark @
/*======================================================================*/ 

load y[t,n]=lnfd.txt;

if level==1; y=y[2:T,.]-y[1:T-1,.]; else; y=y; endif;

  T = rows(y);
  N = cols(y);

/*  output file = bt.out reset; */
 

/* determine the number of lags for each country */
  c  = ones(T,1);
  Lv = -999999*ones(N,1);
  i = 1;
  do while i <= N;
    y0 = y[.,i];
    xv = -999999*ones(T,pb+2);
    xi = pb;
    do while xi > 0;
      xv[.,1] = c;
      xv[.,2] = lagn(y0,1);
      dy = y0 - lagn(y0,1);
      ip = 1;
      do while ip <= xi;
        xv[.,ip+2] = lagn(dy,ip);
        ip = ip + 1;
      endo;
      yv = y0[xi+2:T];
      xv = xv[xi+2:T,.];

       __con=0;
      {vnam,mm,beta,stb,vc,stderr,sigma,cx,rsq,resid,dwstat}= ols(0,yv,xv);
      np = rows(beta);
      stat = beta[np]/stderr[np];

      if abs(stat) > 1.65;
        break;
      else;
        xi = xi - 1;
      endif;
      xv  = -999999*ones(T,xi+2);
    endo;
    Lv[i] = xi;
  i = i + 1;
  endo;



/* calculate test statistics */
  tem1 = 0;
  tem2 = 0;
  tem3 = 0;
  tem4 = 0;
  tem5 = 0;
  tem6 = 0;
  vv=detrend(y,d);
  fr = zeros(N,2);
  i = 1;
  do while i <= N;
    {m1,m2,m3} = adf(y[.,i],d,Lv[i]);
    tem1 = tem1 + m2;
    fr[i,1] = m2;
    {n1,n2,n3} = adf(rev(y[.,i]),d,Lv[i]);
    tem2 = tem2 + n2;
    fr[i,2] = n2;
    tem3 = tem3 + maxc((m2~n2)');
    lq1 = lmt(T,y[.,i],Lv[i],d);
    tem4 = tem4 + lq1;
    lq2 = lmt(T,rev(y[.,i]),Lv[i],d);
    tem5 = tem5 + minc((lq1~lq2)');
    wstest=ws(T,vv[.,i],Lv[i]); 
    tem6 = tem6 + wstest;
    i = i + 1;
  endo;
  rdfa = tem1/N;
  rdfra = tem2/N;
  rmxa = tem3/N;
  rlma = tem4/N;
  rmna = tem5/N;
  rws = tem6/N;

/* compute residuals for bootstrap */
  pmax = maxc(Lv);


  dy  = y - lag(y);
  eta = zeros(pmax,N);
  e0  = zeros(T,N);

  i = 1;
  do while i <= N;
    if Lv[i] > 0;
      yi = dy[.,i];
      xi = zeros(T,Lv[i]);
      pk = 1;
      do while pk <= Lv[i];
        xi[.,pk] = lagn(yi,pk);
        pk = pk + 1;
      endo;
      __con=0;
      {vnam,mm,b,stb,vc,stderr,sigma,cx,rsq,resid,dwstat} = ols(0,yi,xi);  
   eta[1:Lv[i],i] = b;
      e0[.,i] = yi - xi*b;
    else;
      e0[.,i] = dy[.,i];
    endif;
    i = i + 1;
  endo;
  e0 = e0[pmax+2:T,.];



/* recenter bootstrap residuals */
ee0=zeros(N,T-pmax-1);
an=e0';
j=1;
do while j<=T-pmax-1;
ss=sumc(an');
ee0[.,j]=an[.,j]-((1/(T-pmax-1))*ss);
j=j+1;
endo;

e0=ee0';


/* begin bootstrap resampling */
  e0s   = zeros(rows(e0),N);
  btrdfa = zeros(nb,1);
  btrdfra = zeros(nb,1);
  btrmxa = zeros(nb,1);
  btrlma = zeros(nb,1); 
  btrmna = zeros(nb,1);
  btrws = zeros(nb,1);

  bi = 1;
  do while bi <= nb;
  id   = ceil(rndu(T,1)*rows(e0)); 
    e0s  = e0[id,.];


/*calculation of psi's*/

psi0=zeros(pmax,N);
i=1;
do while i<=N;
if Lv[i]>0;
psi0[1,i]=eta[1,i];
k=2;
do while k<=Lv[i];
psi0[k,i]=eta[1:k-1,i]'*rev(psi0[1:k-1,i])+eta[k,i];
k=k+1;
endo;
endif;
i=i+1;
endo;


psi=zeros(m,N);
eps=zeros(m,N);
i=1;
do while i<=N;
if Lv[i]>0;
psi[.,i]=recserar(eps[.,i],psi0[1:Lv[i],i],eta[1:Lv[i],i]);
endif;
i=i+1;
endo;


/* obtain initial values */
    u0  = zeros(pmax,N);

    i = 1;
    do while i <= N;
      if Lv[i] > 0;
        ei0 = e0[.,i];
    id  = ceil(rndu(m+1,1)*rows(ei0));
      eis = ei0[id,.];
        pk  = 1;
        do while pk <= Lv[i];
     psi1=(1|psi[.,i])';    
if pk==1;
          u0[pk,i] = psi1*eis;
endif;
     id1  = ceil(rndu(1,1)*rows(ei0));
      eiss = ei0[id1,.];
eispk=eis[2:m+1]|eiss;
eis=eispk;         
 u0[pk,i] = psi1*eispk;
          pk = pk +1;
        endo;
      endif;
      i = i + 1;
    endo;

/* generate u_star */
    us = zeros(rows(e0s),N);
    ys = zeros(rows(e0s),N);
    i = 1;
    do while i <= N;
      if Lv[i] > 0;
        us[.,i] = recserar(e0s[.,i],u0[1:Lv[i],i],eta[1:Lv[i],i]);
        ys[.,i] = cumsumc(us[.,i]);
      else;
        ys[.,i] = cumsumc(e0s[.,i]);
      endif;
      i = i + 1;
    endo;


/* calculate test statistics */

vv2=detrend(ys,d);


    tem1 = 0;
    tem2 = 0;
    tem3 = 0;
    tem4 = 0;
    tem5 = 0;
    tem6 = 0;

    i = 1;
    do while i <= N;
      {m1,m2,m3} = adf(ys[.,i],d,Lv[i]);
      tem1 = tem1 + m2;
      {n1,n2,n3} = adf(rev(ys[.,i]),d,Lv[i]);
      tem2 = tem2 + n2;
      tem3 = tem3 + maxc((m2~n2)');

      lq1 = lmt(T,ys[.,i],Lv[i],d);

      tem4 = tem4 + lq1;

      lq2 = lmt(T,rev(ys[.,i]),Lv[i],d);

      tem5 = tem5 + minc((lq1~lq2)');
 
     wstest=ws(T,vv2[.,i],Lv[i]); 
 
      tem6 = tem6 + wstest;

      i = i + 1;
    endo;
    btrdfa[bi] = tem1/N;
    btrdfra[bi] = tem2/N;
    btrmxa[bi] = tem3/N;
    btrlma[bi] = tem4/N;
    btrmna[bi] = tem5/N;
    btrws[bi] = tem6/N;

//      print "Replications = " bi;

    bi = bi + 1;
  endo;

/*output file="bootemp" on;*/

  pvf  = sumc(btrdfa  .<= rdfa)/nb;
  pvr  = sumc(btrdfra .<= rdfra)/nb;
  pvmx = sumc(btrmxa  .<= rmxa)/nb;
  pvlm = sumc(btrlma  .>= rlma)/nb;
  pvmn = sumc(btrmna  .>= rmna)/nb;
  pvws = sumc(btrws  .<= rws)/nb;



  print "";
  print "Number of observations(T) : "  T;
  print "Number of individuals(N)  : "  N;
  print "Degree of detrending(d)   : "  d;
  print "Bootstrap Replications    : "  nb;
//  print "Block size                : "  m;
  print "Maximum lag length        : "  pb;
  print "t-bar statistic           : "  rdfa;
  print "p-value for t-bar         : "  pvf;
  print "WS statistic              : "  rws;
  print "P-value for WS            : "  pvws;
/*  print "Max statistic             : "  rmxa;
  print "P-value for Max           : "  pvmx;
  print "LM statistic             : "  rlma;
  print "P-value for LM           : "  pvlm;
  print "MinLM statistic             : "  rmna;
  print "P-value for MinLM           : "  pvmn;

  print "Number of lagged vars(p)  : "  Lv;
*/

  d2 = date;
  et = ethsec(d1,d2);
  et = et/100;
  print "Total elapsed time [mins] :" et/60;
  print "";



/*output file="bootemp" off;*/

proc lmt(T,x,p,d);
     local dx, c, tr, tem, f1 , f2, f3, f4, f5, f6, f7, f8, f9, r1, r2, f11;
     
    __con=1;

     tr = seqa(1,1,T);

     dx = x - lagn(x,1);

     if d == 0;

        if p == 0;

           dx = dx[2:T];

           r1 = dx - meanc(dx);

           r1 = (zeros(1,1)~r1')'; 
           
       
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1));  
   
           
        elseif p == 1;
      
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1));

           r1 = (zeros(1,2)~r1')'; 
   
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1));

        elseif p == 2;
    
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2));

           r1 = (zeros(1,3)~r1')'; 
       
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2));

        elseif p == 3;
    
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3));

           r1 = (zeros(1,4)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3));

        elseif p == 4;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4));

           r1 = (zeros(1,5)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4));

   elseif p == 5;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5));

           r1 = (zeros(1,6)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5));

elseif p == 6;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6));

           r1 = (zeros(1,7)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6));


  elseif p == 7;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7));

           r1 = (zeros(1,8)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7));


  elseif p == 8;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8));

           r1 = (zeros(1,9)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8));


  elseif p == 9;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9));

           r1 = (zeros(1,10)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9));


  elseif p == 10;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9)~lagn(dx,10));

           r1 = (zeros(1,11)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9)~lagn(dx,10));

        endif;

     endif;

     if d == 1;

        if p == 0;

          {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,tr);       

           r1 = (zeros(1,1)~r1')';   

          {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,tr~lagn(x,1));
   
           
        elseif p == 1;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,tr~lagn(dx,1));

           r1 = (zeros(1,2)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,tr~lagn(x,1)~lagn(dx,1));

        elseif p == 2;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,tr~lagn(dx,1)~lagn(dx,2));

           r1 = (zeros(1,3)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,tr~lagn(x,1)~lagn(dx,1)~lagn(dx,2));

        elseif p == 3;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,tr~lagn(dx,1)~lagn(dx,2)~lagn(dx,3));

           r1 = (zeros(1,4)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,tr~lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3));

        elseif p == 4;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,tr~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4));

           r1 = (zeros(1,5)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,tr~lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4));

  elseif p == 5;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5));

           r1 = (zeros(1,6)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5));


  elseif p == 6;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6));

           r1 = (zeros(1,7)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6));


  elseif p == 7;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7));

           r1 = (zeros(1,8)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7));


  elseif p == 8;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8));

           r1 = (zeros(1,9)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8));


  elseif p == 9;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9));

           r1 = (zeros(1,10)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9));


  elseif p == 10;

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r1,f11} = ols(0,dx,lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9)~lagn(dx,10));

           r1 = (zeros(1,11)~r1')'; 

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,r2,f11} = ols(0,r1,lagn(x,1)~lagn(dx,1)~lagn(dx,2)~lagn(dx,3)~lagn(dx,4)~lagn(dx,5)~lagn(dx,6)~lagn(dx,7)~lagn(dx,8)~lagn(dx,9)~lagn(dx,10));

        endif;

     endif;

     tem = 1 - (r2'*r2)/(r1'*r1); 

     retp(T*tem);

endp;




proc ws(n,h,p);
     local r, t, w, dep, ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, ind9, ind10, ind11, ind12, ind13,
              f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, sw;

         __con=0;

         r = p + 1;

         t = seqa(1,1,n);

         w =  (sqrt(t[1:n-r])'~rev(sqrt(t[1:n-r]))')';

         dep = (h[r+1:n]'~h[1:n-r]')';   

         ind0 = (h[r:n-1]'~h[2:n-r+1]')'; 

         if r == 1; 
          
    
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0));
                                                                 
           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1];  

         endif;

         if r == 2;

           ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 


           {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 

        endif;

        if r == 3;

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';


           {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
           
        endif;

        if r == 4; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';


           {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 

        endif;

        if r == 5; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';

	
           {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 

        endif;

 if r == 6; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';
            ind5 = ((h[r-4:n-1-4]-h[r-1-4:n-2-4])'~(h[2+4:n-r+1+4]-h[3+4:n-r+2+4])')';   

           {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4)~(w .* ind5));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
      
        endif;
        
 if r == 7; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';
            ind5 = ((h[r-4:n-1-4]-h[r-1-4:n-2-4])'~(h[2+4:n-r+1+4]-h[3+4:n-r+2+4])')';   
            ind6 = ((h[r-5:n-1-5]-h[r-1-5:n-2-5])'~(h[2+5:n-r+1+5]-h[3+5:n-r+2+5])')';   
           
             {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4)~(w .* ind5)~(w .* ind6));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
      
        endif;

 if r == 8; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';
            ind5 = ((h[r-4:n-1-4]-h[r-1-4:n-2-4])'~(h[2+4:n-r+1+4]-h[3+4:n-r+2+4])')';   
            ind6 = ((h[r-5:n-1-5]-h[r-1-5:n-2-5])'~(h[2+5:n-r+1+5]-h[3+5:n-r+2+5])')';   
            ind7 = ((h[r-6:n-1-6]-h[r-1-6:n-2-6])'~(h[2+6:n-r+1+6]-h[3+6:n-r+2+6])')';
   
             {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4)~(w .* ind5)~(w .* ind6)~(w .* ind7));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
      
        endif;

 if r == 9; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';
            ind5 = ((h[r-4:n-1-4]-h[r-1-4:n-2-4])'~(h[2+4:n-r+1+4]-h[3+4:n-r+2+4])')';   
            ind6 = ((h[r-5:n-1-5]-h[r-1-5:n-2-5])'~(h[2+5:n-r+1+5]-h[3+5:n-r+2+5])')';   
            ind7 = ((h[r-6:n-1-6]-h[r-1-6:n-2-6])'~(h[2+6:n-r+1+6]-h[3+6:n-r+2+6])')';
            ind8 = ((h[r-7:n-1-7]-h[r-1-7:n-2-7])'~(h[2+7:n-r+1+7]-h[3+7:n-r+2+7])')';

             {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4)~(w .* ind5)~(w .* ind6)~(w .* ind7)~(w .* ind8));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
      
        endif;

 if r == 10; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';
            ind5 = ((h[r-4:n-1-4]-h[r-1-4:n-2-4])'~(h[2+4:n-r+1+4]-h[3+4:n-r+2+4])')';   
            ind6 = ((h[r-5:n-1-5]-h[r-1-5:n-2-5])'~(h[2+5:n-r+1+5]-h[3+5:n-r+2+5])')';   
            ind7 = ((h[r-6:n-1-6]-h[r-1-6:n-2-6])'~(h[2+6:n-r+1+6]-h[3+6:n-r+2+6])')';
            ind8 = ((h[r-7:n-1-7]-h[r-1-7:n-2-7])'~(h[2+7:n-r+1+7]-h[3+7:n-r+2+7])')';
            ind9 = ((h[r-8:n-1-8]-h[r-1-8:n-2-8])'~(h[2+8:n-r+1+8]-h[3+8:n-r+2+8])')';

             {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4)~(w .* ind5)~(w .* ind6)~(w .* ind7)
									~(w .* ind8)~(w .* ind9));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
      
        endif;

  if r == 11; 

            ind1 = ((h[r:n-1]-h[r-1:n-2])'~(h[2:n-r+1]-h[3:n-r+2])')'; 
            ind2 = ((h[r-1:n-1-1]-h[r-1-1:n-2-1])'~(h[2+1:n-r+1+1]-h[3+1:n-r+2+1])')';
            ind3 = ((h[r-2:n-1-2]-h[r-1-2:n-2-2])'~(h[2+2:n-r+1+2]-h[3+2:n-r+2+2])')';
            ind4 = ((h[r-3:n-1-3]-h[r-1-3:n-2-3])'~(h[2+3:n-r+1+3]-h[3+3:n-r+2+3])')';
            ind5 = ((h[r-4:n-1-4]-h[r-1-4:n-2-4])'~(h[2+4:n-r+1+4]-h[3+4:n-r+2+4])')';   
            ind6 = ((h[r-5:n-1-5]-h[r-1-5:n-2-5])'~(h[2+5:n-r+1+5]-h[3+5:n-r+2+5])')';   
            ind7 = ((h[r-6:n-1-6]-h[r-1-6:n-2-6])'~(h[2+6:n-r+1+6]-h[3+6:n-r+2+6])')';
            ind8 = ((h[r-7:n-1-7]-h[r-1-7:n-2-7])'~(h[2+7:n-r+1+7]-h[3+7:n-r+2+7])')';
            ind9 = ((h[r-8:n-1-8]-h[r-1-8:n-2-8])'~(h[2+8:n-r+1+8]-h[3+8:n-r+2+8])')';
            ind10 = ((h[r-9:n-1-9]-h[r-1-9:n-2-9])'~(h[2+9:n-r+1+9]-h[3+9:n-r+2+9])')';
            
	 {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11} = ols(0,(w .* dep),(w .* ind0)
                                                                 ~(w .* ind1)~(w .* ind2)~(w .* ind3)~(w .* ind4)~(w .* ind5)~(w .* ind6)~(w .* ind7)
									~(w .* ind8)~(w .* ind9)~(w .* ind10));

           sw = sqrt((n-r-1)/(2*(n-r)-1))*(f3[1] - 1)/f6[1]; 
      
        endif;


        retp(sw);

endp;




end;


