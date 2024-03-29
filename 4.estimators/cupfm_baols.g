new;cls;

t =28;                             
n =9;
k =2;                           /*number of regressors*/

nfmax=2;                        /* max fact < n */
pen  =2;                        /* penalty */
cri  =3;                        /* 1 = pc, 2 = ic, 3 = aic/bic */
p    =int(4*(t/100)^(2/9));     /* bandwidth */

load y[t,n]  =lnfd.txt ;        /* t x n matrix */
load x[t,n*k]=west_x.txt;     /* t x kn matrix */


{b1,adjb1,fmb11,fmb21,s11,s21,nf1} = tsest(x,y,nfmax,p,pen,cri);
{b2,adjb2,fmb12,fmb22,s12,s22,nf2} = itest(x,y,nfmax,p,pen,cri);



/* printing options */
format/m1/rd 6,3;


print " ";
print " t = ";; t;
print " n = ";; n;
print " fn  ";; nf1;

"----------ols----------";
"  beta   se    t-stat";
b1'~s21'~(b1./s21)';             @ OLS estimator proposed by Bai and Kao (2006)@ 

"----------cup-fm----------";
"  beta   se    t-stat";
fmb11'~s11'~(fmb11./s11)';       @ CUP-FM estimator of Bai and Kao (2006)@

"----------ba-ols----------";
"  beta   se    t-stat";
adjb1'~s11'~(adjb1./s11)';       @ bias-adjusted estimator of Westerlund (2007)@

print " ";
/*print " iter beta/standard error: ";
print " ols ";; b2;; s22;
print " adj ";; adjb2;; s12;
print " fm  ";; fmb12;; s12;
print " fma ";; fmb22;; s22;
print " fn  ";; nf2;
*/



/* procs */


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
f         = f[.,1:nf]*sqrt(t);
lam       = (e'f)/t;
else;
{f0,v,f}  = svd1(e'e);
lam       = f[.,1:nf]*sqrt(n);
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



proc (7) = tsest(x,y,fmax,p,pen,cri);
local t,n,xb,yb,xx,b,e,f0,v,dx,f,lam,i,s,si,di,d,l,o,adjb,tr,nf;
local yp,xyp,sp,fmb1,fmb2,fms1,fms2,o11,o12,o22,o13,o31,o23,o32,o33;
local s11,s22,l22,l12,l23,l13,l33,s13,s23,k,xbv,xbl,dxv;


n   = cols(y);
t   = rows(x);
k   = cols(x)/n;
xb  = x-meanc(x)';
yb  = y-meanc(y)';
dx  = x[2:t,.]-x[1:t-1,.];

i = 1;
xbv = zeros(n*t,1);
xbl = zeros(n*(t-1),1);
dxv = zeros(n*(t-1),1);
do while i <= k;

xbv = xbv~vec(xb[.,1+n*(i-1):n*i]);
xbl = xbl~vec(xb[2:t,1+n*(i-1):n*i]);
dxv = dxv~vec(dx[.,1+n*(i-1):n*i]);

i = i+1;
endo;

xbv = xbv[.,2:cols(xbv)];
xbl = xbl[.,2:cols(xbl)];
dxv = dxv[.,2:cols(dxv)];
xx  = inv(xbv'xbv);
b   = (vec(y)'xbv)*xx;
e   = vec(yb) - xbv*b';
e   = reshape(e,n,t)';

nf      = fact(e,fmax,pen,cri);            /* est no of factors */
{f,lam} = prin(e,nf);



i  = 1;
sp = zeros(k,k);
s  = zeros(nf+k+1,nf+k+1);
d  = zeros(nf+k+1,nf+k+1);
e  = e - f*lam';
do while i <= n;

v   = f[2:t,.]~e[2:t,i]~dxv[1+(t-1)*(i-1):(t-1)*i,.];
si  = (v'v)/t;
di  = fejer(v,p);
s   = s + si;
d   = d + di;
o   = si + di + di';

o11 = o[1:nf,1:nf];
o22 = o[nf+1,nf+1];
o13 = o[1:nf,nf+2:nf+k+1];
o31 = o[nf+2:nf+k+1,1:nf];
o23 = o[nf+1,nf+2:nf+k+1];
o32 = o[nf+2:nf+k+1,nf+1];
o33 = o[nf+2:nf+k+1,nf+2:nf+k+1];

s11 = o11-o13*inv(o33)*o31;
s22 = o22-o23*inv(o33)*o32;

sp  = sp + ((lam[i,.]*s11*lam[i,.]')*o33+s22*o33);

i = i+1;
endo;

sp = sp/n;
d  = d/n;
s  = s/n;
l  = s + d;
o  = s + d + d';

o11 = o[1:nf,1:nf];
o12 = o[1:nf,nf+1];
o22 = o[nf+1,nf+1];
o13 = o[1:nf,nf+2:nf+k+1];
o31 = o[nf+2:nf+k+1,1:nf];
o23 = o[nf+1,nf+2:nf+k+1];
o32 = o[nf+2:nf+k+1,nf+1];
o33 = o[nf+2:nf+k+1,nf+2:nf+k+1];

l22 = l[nf+1,nf+1];
l12 = l[1:nf,nf+1];
l23 = l[nf+1,nf+2:nf+k+1];
l13 = l[1:nf,nf+2:nf+k+1];
l33 = l[nf+2:nf+k+1,nf+2:nf+k+1];

i  = 1;
yp = zeros(t-1,n);
do while i <= n;

yp[.,i] = y[2:t,i] - dxv[1+(t-1)*(i-1):(t-1)*i,.]*inv(o33)*(o32+o31*lam[i,.]');

i = i+1;
endo;

s13  = l13-o13*inv(o33)*l33;
s23  = l23-o23*inv(o33)*l33;
s22  = o22-o23*inv(o33)*o32;

adjb = b - ((-1/2)*meanc(lam)'o13+(-1/2)*o23+meanc(lam)'l13+l23)*inv((1/6)*o33)/t;
xyp  = (vec(yp)'xbl);
fmb1 = (xyp - n*t*(meanc(lam)'s13+s23))*xx;
fms1 = sqrt(diag(6*inv(o33)*sp*inv(o33)))';
fmb2 = (xyp - n*t*s23)*xx;
fms2 = sqrt(diag(6*inv(o33)*s22))';

retp(b,adjb,fmb1,fmb2,fms1,fms2,nf);
endp;




proc (7) = itest(x,y,fmax,p,pen,cri);
local t,n,xb,yb,xx,b,olsb,e,f0,v,dx,f,lam,i,s,si,di,d,l,o,adjb,tr,nf;
local yp,xyp,sp,fmb1,fmb2,fms1,fms2,o11,o12,o22,o13,o31,o23,o32,o33;
local s11,s22,l22,l12,l23,l13,l33,s13,s23,max,cr,j,k,xbv,xbl,dxv;


n   = cols(y);
t   = rows(x);
k   = cols(x)/n;
xb  = x-meanc(x)';
yb  = y-meanc(y)';
dx  = x[2:t,.]-x[1:t-1,.];

i = 1;
xbv = zeros(n*t,1);
xbl = zeros(n*(t-1),1);
dxv = zeros(n*(t-1),1);
do while i <= k;

xbv = xbv~vec(xb[.,1+n*(i-1):n*i]);
xbl = xbl~vec(xb[2:t,1+n*(i-1):n*i]);
dxv = dxv~vec(dx[.,1+n*(i-1):n*i]);

i = i+1;
endo;

xbv  = xbv[.,2:cols(xbv)];
xbl  = xbl[.,2:cols(xbl)];
dxv  = dxv[.,2:cols(dxv)];
xx   = inv(xbv'xbv);
olsb = (vec(y)'xbv)*xx;

j   = 1;
max = 20;
cr  = 1000;
b   = olsb;
do while maxc(cr) > 0.00001;

e       = vec(yb) - xbv*b';
e       = reshape(e,n,t)';
nf      = fact(e,fmax,pen,cri);            /* est no of factors */
{f,lam} = prin(e,nf);

i  = 1;
sp = zeros(k,k);
s  = zeros(nf+k+1,nf+k+1);
d  = zeros(nf+k+1,nf+k+1);
e  = e - f*lam';
do while i <= n;

v   = f[2:t,.]~e[2:t,i]~dxv[1+(t-1)*(i-1):(t-1)*i,.];
si  = (v'v)/t;
di  = fejer(v,p);
s   = s + si;
d   = d + di;
o   = si + di + di';

o11 = o[1:nf,1:nf];
o22 = o[nf+1,nf+1];
o13 = o[1:nf,nf+2:nf+k+1];
o31 = o[nf+2:nf+k+1,1:nf];
o23 = o[nf+1,nf+2:nf+k+1];
o32 = o[nf+2:nf+k+1,nf+1];
o33 = o[nf+2:nf+k+1,nf+2:nf+k+1];

s11 = o11-o13*inv(o33)*o31;
s22 = o22-o23*inv(o33)*o32;

sp  = sp + ((lam[i,.]*s11*lam[i,.]')*o33+s22*o33);

i = i+1;
endo;

sp = sp/n;
d  = d/n;
s  = s/n;
l  = s + d;
o  = s + d + d';

o11 = o[1:nf,1:nf];
o12 = o[1:nf,nf+1];
o22 = o[nf+1,nf+1];
o13 = o[1:nf,nf+2:nf+k+1];
o31 = o[nf+2:nf+k+1,1:nf];
o23 = o[nf+1,nf+2:nf+k+1];
o32 = o[nf+2:nf+k+1,nf+1];
o33 = o[nf+2:nf+k+1,nf+2:nf+k+1];

l22 = l[nf+1,nf+1];
l12 = l[1:nf,nf+1];
l23 = l[nf+1,nf+2:nf+k+1];
l13 = l[1:nf,nf+2:nf+k+1];
l33 = l[nf+2:nf+k+1,nf+2:nf+k+1];

i  = 1;
yp = zeros(t-1,n);
do while i <= n;

yp[.,i] = y[2:t,i] - dxv[1+(t-1)*(i-1):(t-1)*i,.]*inv(o33)*(o32+o31*lam[i,.]');

i = i+1;
endo;

s13  = l13-o13*inv(o33)*l33;
s23  = l23-o23*inv(o33)*l33;
s22  = o22-o23*inv(o33)*o32;

adjb = b - ((-1/2)*meanc(lam)'o13+(-1/2)*o23+meanc(lam)'l13+l23)*inv((1/6)*o33)/t;
xyp  = (vec(yp)'xbl);
fmb1 = (xyp - n*t*(meanc(lam)'s13+s23))*xx;
fms1 = sqrt(diag(6*inv(o33)*sp*inv(o33)))';
fmb2 = (xyp - n*t*s23)*xx;
fms2 = sqrt(diag(6*inv(o33)*s22))';
cr   = abs(fmb1 - b);

if j > max;
goto top;                   /* no of iter reached upper limit */
else;
b = fmb1;
j = j+1;
endif;
endo;

top:

retp(olsb,adjb,fmb1,fmb2,fms1,fms2,nf);
endp;











