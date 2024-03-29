new;

/*
** panellmcoint.prg was written by Joakim Westerlund, Department of Economics, Lund University.
** Contact info:  joakim.westerlund@nek.lu.se
** The code produces panel cointegration tets of
   Westerlund and Edgerton (2008).A Simple Test for Cointegration in Dependent Panels with Structural Breaks. Oxford Bulletion of Economics and Statistics 70.
** The code can be used freely as long as proper reference is given. No performance
** guarantee is made. Bug reports are welcome.
*/
k       = 2;    
t       = 28;
n       = 9;
max     = 2;                        /* max factors */
mod     = 2;                        /* 0 = no shift, 
                                       1 = level shift, 
                                       2 = regime shift */
p       = int(4*(t/100)^(2/9));     /* lags */
q       = int(4*(t/100)^(2/9));     /* bandwidth */
tr      = 0.1;                      /* trimming */

load y[t,n]   = lnfd.txt;
load x[t,n*k] = west_x.txt;

brn      = ilt_br(y,x,tr,p,mod);
{lmn,nf} = ilt_fact(y,x,brn,p,q,mod,max);


/* printing options */
format/m1/rd 8,5;

print " ";
print " tau_n   = ";; lmn[1];
print " p-value = ";; (1-cdfnc(lmn[1]));
print " phi_n   = ";; lmn[2];
print " p-value = ";; (1-cdfnc(lmn[2]));
print " ";
print "#Factors = ";; nf;
print " ";

format/m1/rd 8,0;

if mod ne 0; 
    "    country   break date";;
    seqa(1,1,n)~brn';
endif;



/* procs */

proc (1) = ilt_br(y,x,tr,p,mod);
local t,n,i,br,t1,t2,j,ssr;

t = rows(y);
n = cols(y);

i  = 1;
br = zeros(1,n);
t1 = round(tr*t);
t2 = round((1-tr)*t);
do while i <= n;

j   = t1;
ssr = zeros(t,1);
do while j <= t2;
ssr[j] = ilt_ssr(y[.,i],x[.,i],j,p,mod);
j      = j + 1;
endo;

br[i] = minpr(ssr[t1:t2]) + t1 - 1;

i = i + 1;
endo;

retp(br);
endp;



proc (2) = ilt_fact(y,x,br,p,q,mod,max);
local t,n,i,de,z,dz,dy,dx,nf,f,l,d,b,s0,sl,ds,v,u,v0,vl,zt,za;

t = rows(y);
n = cols(y);

i = 1;
de = zeros(t-1,n);
do while i <= n;

z       = dum(x[.,i],br[i],mod);
dz      = diff(z,1);
dy      = diff(y[.,i],1);
dx      = (diff(x[.,i],1)~dz);
de[.,i] = (eye(t-1) - dx*inv(dx'dx)*dx')*dy;

i = i + 1;
endo;

i     = 1;
zt    = 0;
za    = 0;
de    = (zeros(1,n)|de);
if max eq 0;
nf    = 0;
f     = zeros(t,n);
else;
nf    = fact(de,max);
{f,l} = prin(de,nf);
f     = cumsumc(f*l');
endif;
do while i <= n;

z  = dum(x[.,i],br[i],mod);
dz = diff(z,1);
dx = (diff(x[.,i],1)~dz);
dy = diff(y[.,i],1);
d  = inv(dx'dx)*(dx'dy);
s0 = y[.,i] - (y[1,i] - (x[1,i]~z[1,.])*d) - (x[.,i]~z)*d  - f[.,i];
ds = diff(s0,1);
sl = (s0[p+1:t-1,.]~lagp(ds,p)~dz[p+1:t-1,1]);
ds = ds[p+1:t-1];
d  = inv(sl'sl)*sl'ds;
u  = ds - sl*d;
v0 = lrvar(u,0);
vl = lrvar(ds,q);
v  = sqrt(diag(v0*inv(sl'sl)));

zt = zt + (d[1]/v[1]);
za = za + ((t-p-1)*d[1])*sqrt(vl/v0);

i = i + 1;
endo;

zt = sqrt(n)*(zt/n+1.9675)/sqrt(0.3301);
za = sqrt(n)*(za/n+8.4376)/sqrt(25.8964);

retp((zt~za),nf);
endp;



proc(1) = fact(e,nf);
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
pen     = (n+t)/(n*t)*log(minc(n|t));
cr[k]   = log(s) + k*pen;

k = k + 1;
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


proc (1) = detm(x,br,mod);
local t,d,z,dz;

t = rows(x);

if br ne 0;
d = (zeros(br,1)|ones(t-br,1));
endif;

if mod eq 0 or mod eq 1;
z = seqa(1,1,t);
endif;

if mod eq 2;
if br eq 0; z = seqa(1,1,t);        endif;
if br ne 0; z = seqa(1,1,t)~(d.*x); endif;
endif;

retp(z);
endp;



proc (1) = dum(x,br,mod);
local t,d,z,dz;

t = rows(x);

if br ne 0;
d = (zeros(br,1)|ones(t-br,1));
endif;

if mod eq 0;
z = seqa(1,1,t);
endif;

if mod eq 1;
if br eq 0; z = seqa(1,1,t);   endif;
if br ne 0; z = seqa(1,1,t)~d; endif;
endif;

if mod eq 2;
if br eq 0; z = seqa(1,1,t);          endif;
if br ne 0; z = seqa(1,1,t)~d~(d.*x); endif;
endif;

retp(z);
endp;


proc (1) = ilt(y,x,br,p,mod);
local k,d,b,u,za,zt,t,dy,dx,tr,v,z,dz,s,s0,sl,ds,v0,vl;

t = rows(y);

z  = dum(x,br,mod);
dz = diff(z,1);
dx = (diff(x,1)~dz);
dy = diff(y,1);
d  = inv(dx'dx)*(dx'dy);
s0 = y - (y[1] - (x[1,.]~z[1,.])*d) - (x~z)*d;
ds = diff(s0,1);
sl = (s0[p+1:t-1,.]~lagp(ds,p)~dz[p+1:t-1,.]);
ds = ds[p+1:t-1];
d  = inv(sl'sl)*sl'ds;
u  = ds - sl*d;
v0 = lrvar(u,0);
vl = lrvar(ds,p);
v  = sqrt(diag(v0*inv(sl'sl)));

zt = d[1]/v[1];
za = ((t-p-1)*d[1])*sqrt(vl/v0);

retp(zt~za);
endp;




proc (1) = ilt_ssr(y,x,br,p,mod);
local d,u,t,dy,dx,z,dz,b;

t = rows(y);

z  = dum(x,br,mod);
dz = diff(z,1);
dx = (diff(x,1)~dz);
dy = diff(y,1);
d  = inv(dx'dx)*(dx'dy);
u  = dy - dx*d;

retp(u'u);
endp;


proc (1) = lagp(x,p);
local t,i,xl;

t = rows(x);

i  = 1;
xl = lagn(x,1);
do while i < p;
i  = i + 1;
xl = (xl~lagn(x,i));
endo;

retp(trimr(xl,p,0));
endp;



proc lrvar(u,k);
local sl,s0;

sl = fejer(u,k);
s0 = (u'u)/rows(u);

retp((s0 + sl + sl'));
endp;



proc fejer(u,k);
local i,w,s,u0,ul;

if k == 0; s = 0; goto out; endif;

i = 1;
s = 0;
do until i > k;
w  = 1 - i/(k+1);
u0 = trimr(u,i,0);
ul = trimr(lagn(u,i),i,0);
s  = s + w*(u0'ul);
i  = i + 1;
endo;

out:

retp(s/rows(u));
endp;


proc diff(x,k);
if k == 0;
retp(x) ;
endif ;
retp(trimr(x,k,0)-trimr(lagn(x,k),k,0));
endp;



proc (1) = minpr(x);
local d,m,i,minx;

i = 1;
m = 1;
d = x[1];
do while i <= rows(x);
if x[i] < d;  d = x[i];  endif;
i = i+1;
endo;

i = 2;
do while i <= rows(x);
if d == x[i]; m = i;  goto stops; endif;
i = i + 1;
endo;

stops:

retp(m);
endp;


proc lagn(x,n);
local y;
y = shiftr(x',n,(miss(0,0))');
retp(y');
endp;


