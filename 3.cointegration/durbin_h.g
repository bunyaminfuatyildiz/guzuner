//new;

/*
** durbinh.prg was written by Joakim Westerlund, Department of Economics, Lund University.
** Contact info:  joakim.westerlund@nek.lu.se
** 
** The code can be used freely as long as proper reference is given. No performance
** guarantee is made. Bug reports are welcome.
*/


k       = 2;                       /* regressors*/
n       = 9;                       /* number of cross-section */  
t       = 28;                      /* time period*/ 
kmax    = 2;                       /* max fact */
cri     = 2;                       /* criterion
                                      2=AIC; 3=BIC*/
pen     = 2;                       /* penalty */
p       = int(4*(t/100)^(2/9));    /* bandwidth */
mod     = 2;                       /* 0=nothing, 1=const, 2=const & trend */   

load y[t,n]  = lnfd.txt;        /* t x n matrix */
load x[t,n*k]= west_x.txt;


if mod == 0; x = x;                        endif;
if mod == 1; x = ones(t,1)~x;              endif;
if mod == 2; x = ones(t,1)~seqa(1,1,t)~x;  endif;

i       = 1; 
dhg     = 0;
{e,nf}  = cum(x,y,kmax,pen,cri);
do until i > n;

dhg = dhg + gdh(e[.,i],p);

i = i + 1;
endo;

mu  = 5.5464;
var = 36.7673;
dhg = sqrt(n)*(dhg/n-mu)/sqrt(var);

mu  = 0.5005;
var = 0.3348;
dhp = sqrt(n)*(pdh(e,p)/n-mu^(-1))/sqrt(mu^(-4)*var);   




    
/* printing options */
format/m1/rd 8,3;


print " dh_g = ";; dhg;;    "p-value=";;(cdfn(dhg));
print " dh_p = ";; dhp;;    "p-value=";;(cdfn(dhp));





/* procs */


proc (2) = cum(x,y,kmax,pen,cri);
local t,n,k,i,e,de,dy,dx,dyi,dxi,nf,f,lam;

t   = rows(y);
n   = cols(y);
k   = cols(x)/n;

i  = 1;
de = zeros(t,n);
y  = zeros(1,n)|y;
x  = zeros(1,k*n)|x;
dy = diff(y,1);
dx = diff(x,1);
do until i > n;

dyi     = dy[.,i];
dxi     = dx[.,1+(i-1)*k:i*k];
de[.,i] = (eye(t) - dxi*inv(dxi'dxi)*dxi')*dyi;  

i = i + 1;
endo;

nf      = fact(de,kmax,pen,cri);   
{f,lam} = prin(de,nf);
de      = de - f*lam';
e       = cumsumc(de);

retp(e,nf);
endp;




proc (1) = gdh(w,p);
local t,w0,wl,v,m,b1,b2,e,s,io,sig,dhs;

t   = rows(w);

wl  = w[1:t-1];
w0  = w[2:t];

b1  = inv(wl'wl)*(wl'w0);
b2  = inv(wl'w0)*(w0'w0);
e   = w0 - wl*b1;
s  	= (e'e)/rows(e);
if p == 0;  io = 0;     	    
else;       io 	= fejer(e,p);   endif;
sig = s + io + io';		
s   = s^2/sig;

dhs = (b2 - b1)^2/(s*inv(wl'wl));

retp(dhs);
endp;





proc (1) = pdh(w,p);
local t,n,w0,wl,v,m,b1,b2,e,s,io,l,lnt,snt,dhs;

t   = rows(w);
n   = cols(w);

w   = w[.,2:cols(w)]; 
wl  = w[1:t-1,.];
w0  = w[2:t,.];
wl  = vec(wl);
w0  = vec(w0);

b1  = inv(wl'wl)*(wl'w0);
b2  = inv(wl'w0)*(w0'w0);

i   = 1;
snt = 0;
lnt = 0;
e   = reshape((w0-wl*b1),n,t-1)';
do until i > n;

s  	= (e[.,i]'e[.,i])/rows(e);
if p == 0;  io = 0;     	    
else;       io = fejer(e[.,i],p);   endif;

snt = snt + s;	                
lnt = lnt + (s + io + io');

i = i + 1;
endo;

dhs = ((lnt/n)*(b2 - b1)^2)/(((snt/n)^2)*inv(wl'wl));

retp(dhs);
endp;







proc lagn(x,n);
local y;
y = shiftr(x',n,(miss(0,0))');
retp(y');
endp;


proc diff(x,k);
if k == 0;
retp(x) ;
endif ;
retp(trimr(x,k,0)-trimr(lagn(x,k),k,0));
endp;





proc(1) = fact(e,nf,p,c);
local t,n,pen,cr,k,s,smax,u,f,lam;

t       = rows(e);
n       = cols(e);   
{f,lam} = prin(e,nf);
u       = e - f*lam';
smax    = sumc(sumc(u.^2))/(n*t);

k    = 1;
cr   = zeros(nf,1);  
do while k <= nf;

{f,lam} = prin(e,k);
u       = e - f*lam';
s       = sumc(sumc(u.^2))/(n*t);

if c == 1 or c == 2;
if p == 1; pen = (n+t)/(n*t)*log((n*t)/(n+t));  endif;
if p == 2; pen = (n+t)/(n*t)*log(minc(n|t));    endif;
if p == 3; pen = (log(minc(n|t))/minc(n|t));    endif;
endif;

if c == 3;
if p == 1; pen = 2*(n+t-k)/(n*t);           endif;  /* aic */
if p == 2; pen = (n+t-k)*log(n*t)/(n*t);    endif;  /* bic */
endif;

if c == 1; cr[k] = s + k*smax*pen;      endif;  /* pc */
if c == 2; cr[k] = log(s) + k*pen;      endif;  /* ic */
if c == 3; cr[k] = s + k*smax*pen;      endif;  /* aic/bic */

k = k+1;
endo;

cr = sortc(seqa(1,1,rows(cr))~cr,2);

retp(cr[1,1]);
endp;




proc (2) = prin(e,nf);
local t,n,f0,v,f,lam;

t = rows(e);
n = cols(e);

if n > t;
{f0,v,f}  = svd1(e*e');
f         = f0[.,1:nf]*sqrt(t);
lam       = (e'f)/t;
else;

{f0,v,f}  = svd1(e'e);
lam       = f0[.,1:nf]*sqrt(n);
f         = (e*lam)/n;
endif;

retp(f,lam);
endp;





proc fejer(uv,k);
local i,m,a,t1,t2,f;
i = 1;
a = 0;
do until i > k;
f = i/(k+1);
m = 1 - f;
t1 = trimr(uv,i,0);
t2 = trimr(lagn(uv,i),i,0);
a = a + m*(t1't2);
i = i + 1;
endo;
retp(a/rows(uv));
endp;
