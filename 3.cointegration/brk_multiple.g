new;

/*
** does the bootstr lm panel coint test
*/
k       = 2;    /*number of regressors (kmax=5)*/
t       = 28;   /*time dimension*/
n       = 9;    /*number of cross-section*/
mod     = 4;    /*1= constant (no break model)
                  2= constant & trend (no break model)
                  3= break in constant (level shift)
                  4= break in const & trend (trend shift)*/
p       = int(4*(t/100)^(2/9));  /* bandwidth */
max     = 3;                     /* max no of breaks */
tr      = int(0.2*t);            /* trimming */
nboot   = 100;                  /*number of bootstrap replications*/

load y[t,n]    = lnfd.txt;
load x[t,n*k]  = west_x.txt;

{lmn,brn} = boot_lm(y,x,mod,tr,p,max,nboot);


/* printing options */
format/m1/rd 8,3;

"  LMstat  asym.pval  boot.pval";lmn;

"  country  break number  -------break dates------";
 format/m1/rd 8,0;
seqa(1,1,n)~brn';



/* procs */



proc(2) = varlags(x,p);
local xp;

xp = shiftr((ones(1,p).*.x)',seqa(1-p,1,p).*.ones(cols(x),1),miss(0,0))';

retp(trimr(x,p,0),trimr(xp,0,p));
endp;




proc(1) = lag_filt(x,rho,p);
local t,n,y,yl,i,j;

t = rows(x);
n = cols(x);

i = p+1;
y = zeros(t+p,n);
do while i <= t+p;

j  = 1;
y[i,.] = x[i-p,.];
do while j <= p;
if rows(rho) eq 1 and cols(rho) eq 1; y[i,.] = y[i,.] + y[i-j,.]*rho;                   endif;
if rows(rho) ne 1 or  cols(rho) ne 1; y[i,.] = y[i,.] + y[i-j,.]*rho[1+n*(j-1):n*j,.];  endif;
j  = j + 1;
endo;

i = i + 1;
endo;

retp(y[p+1:rows(y),.]);
endp;


proc (2) = ols(x,p);
local k,x0,xl,b,u;

k  = cols(x);
x  = x - meanc(x)';
xl = shiftr((ones(1,p).*.x)',(seqa(1-p,1,p).*.ones(k,1)),miss(0,0))';
x0 = trimr(x,p,0);
xl = trimr(xl,0,p);
b  = inv(xl'xl)*(xl'x0);
u  = x0 - xl*b;

retp(b,u);
endp;




proc(2) = boot_lm(y,x,mod,tr,p,m,nb);
local t,n,u,d,brn,br,b,bn,k,boot,i,j,e,yb,xb,rho,lmn,pr,w,r;

t = rows(y);
n = cols(y);
r = cols(x)/n;

i   = 1;
brn = zeros(m+1,n);
if mod < 3; bn = zeros(mod+r,n);             endif;
if mod > 2; bn = zeros((mod-2)*(m+1)+r,n);   endif;
rho = zeros((r+1)*p,(r+1)*n);
u   = zeros(t-p-1,(r+1)*n);
do while i <= n;

brn[.,i]                                                    = breaks(y[.,i],x[.,1+(i-1)*r:i*r],mod,tr,m);
if mod < 3; {bn[1:mod+r,i],e}                               = fm(y[.,i],x[.,1+(i-1)*r:i*r],mod,p,brn[1,i],brn[2:m+1,i]); endif;
if mod > 2; {bn[1:(brn[1,i]+1)*(mod-2)+r,i],e}              = fm(y[.,i],x[.,1+(i-1)*r:i*r],mod,p,brn[1,i],brn[2:m+1,i]); endif;
//{rho[.,1+(r+1)*(i-1):(r+1)*i],u[.,1+(r+1)*(i-1):(r+1)*i]}   = yulew((e~(x[2:t,1+(i-1)*r:i*r]-x[1:t-1,1+(i-1)*r:i*r])),p);
{rho[.,1+(r+1)*(i-1):(r+1)*i],u[.,1+(r+1)*(i-1):(r+1)*i]}   = ols((e~(x[2:t,1+(i-1)*r:i*r]-x[1:t-1,1+(i-1)*r:i*r])),p);

i = i + 1;
endo;

k    = 1;
boot = zeros(nb,1);
do while k <= nb;

i  = 1;
j  = ceil(rndu(t,1)*(rows(u)));
yb = zeros(t,n);
xb = zeros(t,n*r);
do while i <= n;

e                   = u[j,1+(r+1)*(i-1):(r+1)*i];
e                   = e - meanc(e)';
e                   = lag_filt(e,rho[.,1+(r+1)*(i-1):(r+1)*i],p);
xb[.,1+(i-1)*r:i*r] = lag_filt(e[.,2:r+1],1,1);
w                   = xb[.,1+(i-1)*r:i*r]~detc(mod,t,brn[1,i],brn[2:m+1,i]);
b                   = selif(bn[.,i],bn[.,i] .ne 0);
yb[.,i]             = w*b + e[.,1];

i = i + 1;
endo;

{boot[k],br} = lm_n(yb,xb,mod,tr,p,m);

k = k + 1;
endo;

{lmn,br} = lm_n(y,x,mod,tr,p,m);
pr       = (counts(sortc(boot,1),lmn)/nb);
retp(lmn~cdfnc(lmn)~(1-pr),br);
endp;




proc (2) = yulew(x,p);
local t,i,j,vl,v0,x0,xl,u,b,k;

t  = rows(x);
k  = cols(x);

i  = 1;
x  = x - meanc(x)';
v0 = zeros(k*p,k*p);
vl = zeros(k*p,k);
do until i > p;

j  = 1;
do until j > p;
if i eq j; v0[1+(i-1)*k:i*k,1+(j-1)*k:j*k] = (x'x);                         endif;
if i ne j; v0[1+(i-1)*k:i*k,1+(j-1)*k:j*k] = (x[j:t-i+1,.]'x[i:t-j+1,.]);   endif;
vl[1+(i-1)*k:i*k,.]                        = (x[i+1:t,.]'x[1:t-i,.]);
j  = j+1;
endo;

i  = i+1;
endo;

b  = inv(v0)*vl;
xl = shiftr((ones(1,p).*.x)',(seqa(1-p,1,p).*.ones(k,1)),miss(0,0))';
x0 = trimr(x,p,0);
xl = trimr(xl,0,p);
u  = x0 - xl*b;

retp(b,u);
endp;




proc (2) = fm(y,x,mod,p,m,tb);
local w,b,dx,e,vl,v0,s,d,t,u,k;

t  = rows(y);
k  = cols(x);

w  = x~detc(mod,t,m,tb);
b  = inv(w'w)*w'y;
u  = y - w*b;
dx = x[2:t,.]-x[1:t-1,.];
e  = (u[2:t]~dx);
vl = fejer(e,p);
v0 = (e'e)/rows(e);
s  = v0 + vl + vl';
d  = v0 + vl';
d  = t*((d[2:k+1,1]-d[2:k+1,2:k+1]*inv(s[2:k+1,2:k+1])*s[2:k+1,1])|zeros(cols(w)-k,1));

y  = y[2:t] - dx*inv(s[2:k+1,2:k+1])*s[2:k+1,1];
x  = w[2:t,.];
b  = inv(x'x)*((x'y)-d);
u  = y - x*b;

retp(b,u);
endp;



proc(2) = lm_n(y,x,mod,tr,p,m);
local n,t,k,i,lmn,brn,mon;

n = cols(y);
t = rows(y);
k = cols(x)/n;

i    = 1;
lmn  = 0;
brn  = zeros(m+1,1);
mon  = zeros(1,2);
do until i > n;

brn = brn~breaks(y[.,i],x[.,1+(i-1)*k:i*k],mod,tr,m);
lmn = lmn + lm_b(y[.,i],x[.,1+(i-1)*k:i*k],mod,brn[1,1+i],brn[2:m+1,1+i],p);
mon = mon + mom(mod,t,k,brn[1,1+i]);

i = i + 1;
endo;

lmn = sqrt(n)*(lmn/n-mon[1]/n)./sqrt(mon[2]/n);

retp(lmn,brn[.,2:cols(brn)]);
endp;



proc(1) = lm_b(y,x,mod,m,tb,p);
local lmn,yb,xb,j;

lmn = 0;

if m > 0;  tb = selif(tb,tb .> 0); endif;
if m == 0; tb = 0;                 endif;

if mod == 1 or mod == 2;
lmn = lmn + lm(y,x,mod,p);
endif;

if mod == 3 or mod == 4;
if m > 0;
yb  = y[1:tb[1]];
xb  = x[1:tb[1],.];
lmn = lmn + lm(yb,xb,mod,p);

j  = 2;
do while j <= m;
yb  = y[1+tb[j-1]:tb[j]];
xb  = x[1+tb[j-1]:tb[j],.];
lmn = lmn + lm(yb,xb,mod,p);
j   = j + 1;
endo;

yb  = y[1+tb[m]:rows(y)];
xb  = x[1+tb[m]:rows(x),.];
lmn = lmn + lm(yb,xb,mod,p);
endif;

if m == 0;
lmn = lmn + lm(y,x,mod,p);
endif;
endif;

retp(lmn);
endp;







proc(1) = lm(y,x,mod,p);
local t,z,b,e,s,io,lm;

t     = rows(y);
{b,e} = fm(y,x,mod,p,0,0);
s     = (e'e)/t;
io    = fejer(e,p);
s     = s + io + io';
lm    = (sumc(abs(cumsumc(e)).^2))/((t^2)*s);

retp(lm);
endp;



proc(1) = mom(mod,t,k,m);
local a1,a2,a3,mu,var;

if mod == 1 or mod == 3;
a2 = {  0.11708	-0.05532	0.31233	 3.14865 ,
        0.08557	-0.06494	0.24917	 6.69821 ,
        0.06730	-0.11746	0.60610	 6.52534 ,
        0.05391	-0.09801	0.49008	11.09135 ,
        0.04451	-0.09041	0.53038	13.72616 };

a3 = {  0.01130	-0.02838	-0.02232	 0.55336 ,
        0.00581	-0.02252	 0.01171	 0.30932 ,
        0.00309	-0.02008	 0.05684	-0.21322 ,
        0.00161	-0.00995	 0.01394	 0.19668 ,
        0.00095	-0.00646	 0.00484	 0.31930 };
endif;

if mod == 2 or mod == 4;
a2 = {  0.05515	 0.01680	0.21483	 4.98491 ,
        0.04749	-0.02300	0.39897	 6.46768 ,
        0.04103	-0.03676	0.45948	 9.02230 ,
        0.03558	-0.04270	0.53224	11.23719 ,
        0.03073	-0.02145	0.41501	15.37253 };

a3 = {  0.00125	-0.00527	 0.00590	0.02762 ,
        0.00085	-0.00385	 0.00168	0.10582 ,
        0.00057	-0.00232	-0.00541	0.22323 ,
        0.00036	-0.00087	-0.01211	0.33811 ,
        0.00025	-0.00023	-0.01607	0.46277 };
endif;

a1  = 1~t^(-1/2)~t^(-1)~t^(-2);
mu  = a1*a2[k,.]';
var = a1*a3[k,.]';

retp((m+1)*mu~((m+1)^2)*var);
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


proc(1) = detc(mod,t,m,tb);
local du,dt,i,z;

if mod == 1 or mod == 3 and m == 0;
z = ones(t,1);
endif;

if mod == 2 or mod == 4 and m == 0;
z = ones(t,1)~seqa(1,1,t);
endif;

if mod == 3 and m > 0;
i  = 1;
tb = selif(tb,tb .> 0);
du = zeros(t,1);
do until i > m;
du = du~(zeros(tb[i],1)|ones(t-tb[i],1));
i  = i + 1;
endo;
z = ones(t,1)~du[.,2:m+1];
endif;

if mod == 4 and m > 0;
i  = 1;
tb = selif(tb,tb .> 0);
du = zeros(t,1);
dt = zeros(t,1);
do until i > m;

du = du~(zeros(tb[i],1)|ones(t-tb[i],1));
dt = dt~(zeros(tb[i],1)|seqa(1,1,t-tb[i]));

i = i + 1;
endo;
z = ones(t,1)~seqa(1,1,t)~du[.,2:m+1]~dt[.,2:m+1];
endif;

retp(z);
endp;





proc(1) = breaks(y,x,mod,seg,m);
local t,tb,nbr,ssrv,s0;

t = rows(y);

if mod == 1 or mod == 2;
tb  = zeros(m,1);
nbr = 0;
endif;

if mod == 3 or mod == 4;
{ssrv,tb} = dat(y,x,mod,seg,m);
s0        = ssr0(y,x,mod);
nbr       = order(s0,ssrv,t,m,mod);
endif;

if nbr > 0;  tb = tb[.,nbr];  endif;
if nbr == 0; tb = zeros(m,1); endif;

retp(nbr|tb);
endp;



proc(2) = dat(y,x,mod,seg,m);
local t,i,j,z,xbar,zbar,te0,te1,de1,be1,ssr1;
local ssrn,ssr0,dt0,dt,l;

t = rows(y);

if mod == 3; z = ones(t,1);               endif;
if mod == 4; z = ones(t,1)~seqa(1,1,t);   endif;

i       = 1;
ssr0    = zeros(m,1);
dt0     = zeros(m,m);
dt      = zeros(m,m);
do while i <= m;

dt      = dati(y,(x~z),seg,i);
xbar    = pzbar(x,i,dt[1:i,i]);
zbar    = pzbar(z,i,dt[1:i,i]);
te0     = olsqr(y,(zbar~xbar));
de1     = te0[1:cols(z)*(i+1),1];
be1     = olsqr((y-zbar*de1),x);
ssr1    = (y-x*be1-zbar*de1)'(y-x*be1-zbar*de1);

j = 1;
l = 99999999;
do while l > 0.0001;

dt   = dati((y-x*be1),z,seg,i);
zbar = pzbar(z,i,dt[1:i,i]);
te1  = olsqr(y,(x~zbar));
be1  = te1[1:cols(x)];
de1  = te1[cols(x)+1:cols(x)+cols(z)*(i+1)];
ssrn = (y-(x~zbar)*te1)'(y-(x~zbar)*te1);
l    = abs(ssrn-ssr1);

if j >= 50;
goto top;
else;
j = j+1;
ssr1        = ssrn;
ssr0[i]     = ssrn;
dt0[1:i,i]  = dt[1:i,i];
endif;
endo;

i = i+1;
endo;

top:

retp(ssr0,dt0);
endp;




proc(1) = dati(y,z,seg,m);
local i,t,j,l,dt0,odt,ossr,d,ssrmin,dt;
local ssr0,ssr1,ssrv;

i    = 1;
t    = rows(y);
dt0  = zeros(m,m);
odt  = zeros(t,m);
ossr = zeros(t,m);
d    = zeros(t,1);
ssr0 = zeros(m,1);
ssrv = zeros(t*(t+1)/2,1);
do while i <= t-seg+1;
ssr1 = ssr(i,t,y,z,seg);
ssrv[(i-1)*t+i-(i-1)*i/2:i*t-(i-1)*i/2,1] = ssr1[i:t];
i = i+1;
endo;

if m == 1;
{ssrmin,dt} = parti(1,seg,t-seg,t,ssrv,t);
dt0[1,1]    = dt;
ssr0[1,1]   = ssrmin;
else;

j = 2*seg;
do while j <= t;
{ssrmin,dt} = parti(1,seg,j-seg,j,ssrv,t);
ossr[j,1]   = ssrmin;
odt[j,1]    = dt;
j = j+1;
endo;

i = 2;
ssr0[1,1] = ossr[t,1];
dt0[1,1]  = odt[t,1];
do while i <= m;

if i == m;
l  = t;
dt = i*seg;
do while dt <= l-seg;
d[dt,1] = ossr[dt,i-1]+ssrv[(dt+1)*t-dt*(dt+1)/2,1];
dt      = dt+1;
endo;

ossr[l,i] = minc(d[i*seg:l-seg,1]);
odt[l,i]  = (i*seg-1)+minindc(d[i*seg:l-seg,1]);
else;

l = (i+1)*seg;
do while l <= t;

dt = i*seg;
do while dt <= l-seg;
d[dt,1] = ossr[dt,i-1]+ssrv[dt*t-dt*(dt-1)/2+l-dt,1];
dt      = dt+1;
endo;

ossr[l,i] = minc(d[i*seg:l-seg,1]);
odt[l,i]  = (i*seg-1)+minindc(d[i*seg:l-seg,1]);

l = l+1;
endo;
endif;

j = 1;
dt0[i,i] = odt[t,i];
do while j <= i-1;
l = i-j;
dt0[l,i] = odt[dt0[l+1,i],l];
j = j+1;
endo;

ssr0[i,1] = ossr[t,i];
i = i+1;
endo;
endif;

retp(dt0);
endp;




proc(1) = ssr(d0,d1,y,z,seg);
local s,del1,del2,inv1,inv2,inv3,u,v,f,r;

s             = zeros(d1,1);
inv1          = inv(z[d0:d0+seg-1,.]'z[d0:d0+seg-1,.]);
del1          = inv1*(z[d0:d0+seg-1,.]'y[d0:d0+seg-1,1]);
u             = y[d0:d0+seg-1,1]-z[d0:d0+seg-1,.]*del1;
s[d0+seg-1,1] = u'u;

r = d0+seg;
do while r <= d1;

v       = y[r,1]-z[r,.]*del1;
inv3    = inv1*z[r,.]';
f       = 1+z[r,.]*inv3;
del2    = del1+(inv3*v)/f;
inv2    = inv1-(inv3*inv3')/f;
inv1    = inv2;
del1    = del2;
s[r]    = s[r-1]+v*v/f;

r = r+1;
endo;

retp(s);
endp;



proc(2) = parti(s0,b0,b1,s1,ssrv,t);
local k,d,j,s,i,l;

j = b0;
d = zeros(t,1);
i = (s0-1)*t-(s0-2)*(s0-1)/2+1;
do while j <= b1;
l    = j-s0;
k    = j*(t-1)-(j-1)*j/2+s1;
d[j] = ssrv[i+l]+ssrv[k];
j    = j+1;
endo;

s = minc(d[b0:b1]);
d = (b0-1)+minindc(d[b0:b1]);

retp(s,d);
endp;




proc pzbar(z,i,dt);
local t,q,zd,j;

j   = 2;
t   = rows(z);
q   = cols(z);
zd  = zeros(t,(i+1)*q);
zd[1:dt[1],1:q] = z[1:dt[1],.];
do while j <= i;
zd[dt[j-1]+1:dt[j],(j-1)*q+1:j*q] = z[dt[j-1]+1:dt[j],.];
j = j+1;
endo;

zd[dt[i]+1:t,i*q+1:(i+1)*q] = z[dt[i]+1:t,.];

retp(zd);
endp;



proc ssr0(y,x,mod);
local t,z,b,ssr;

t = rows(y);

if mod == 3; z = ones(t,1);               endif;
if mod == 4; z = ones(t,1)~seqa(1,1,t);   endif;

b   = olsqr(y,(z~x));
ssr = (y-(z~x)*b)'(y-(z~x)*b);

retp(ssr);
endp;


proc(1) = order(s0,ssr0,t,m,mod);
local bic,lwz,i,ssr,q;

i            = 0;
q            = mod-2;
ssr          = zeros(m+1,1);
ssr[1,1]     = s0;
ssr[2:m+1,1] = ssr0;
bic          = zeros(m+1,1);
lwz          = zeros(m+1,1);
do while i <= m;
bic[i+1,1] = ln(ssr[i+1,1]/t)+ln(t)*i*(q+1)/t;
lwz[i+1,1] = ln(ssr[i+1,1]/(t-(i+1)*q-i))+(i*(q+1)*(0.299)*(ln(t))^(2.1))/t;
i = i+1;
endo;

bic = minindc(bic)-1;
lwz = minindc(lwz)-1;

retp(bic);
endp;




