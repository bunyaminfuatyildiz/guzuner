new;

/*
** cointboot.prg was written by Joakim Westerlund, Department of Economics, Lund University.
** Contact info: joakim.westerlund@nek.lu.se
** 
** The code can be used freely as long as proper reference is given. No performance
** guarantee is made. Bug reports are welcome.
*/


k   = 2;                       /* regressors, k = 5 is max */
n   = 9;                       /* number of cross-section */  
t   = 28;                      /* time period*/ 
nb  = 1000;                    /* no of bootstrap replications*/
p   = int(4*(t/100)^(2/9));    /* lag length and bandwidth*/
mod = 2;                       /* 1 = constant, 2 = constant and trend*/
est = 1;                       /* sieve estimation, 1 = ols, 2 = yule-walker*/


load y[t,n] 	= lnfd.txt;        // t x n matrix //
load x[t,n*k] 	= west_x.txt;        // t x n matrix //



i = 1;
l = 0;
do while i <= n;
l = l + lm(y[.,i],x[.,1+(i-1)*k:i*k],mod,p);
//"contry-stats";;l;
i = i + 1;
endo;

{m,v} = mom(mod,k);
lmn   = sqrt(n)*(l/n-m)./sqrt(v);
bpv   = boot(y,x,mod,est,nb,p);
bpv   = counts(bpv,lmn)/nb;



/* printing options */
format/m1/rd 8,3;


print " lm statistic = ";; lmn;
print " bootst p-val = ";; (1-bpv);
print " asymp p-val  = ";; cdfnc(lmn);  





/* procs */


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




proc (1) = boot(y,x,mod,est,nb,q);
local n,t,i,b,e,p,u,j,l,nl,h,xb,yb,m,v,k;

n = cols(y);
t = rows(y);
k = cols(x)/n;

i = 1;
b = zeros(mod+k,n);
e = zeros(t-1,n);
p = zeros((k+1)*q,(k+1)*n);
u = zeros(t-q-1,(k+1)*n);
do while i <= n;
{b[.,i],e[.,i]} = fm(y[.,i],x[.,1+(i-1)*k:i*k],mod,q);
if est eq 1; {p[.,1+(k+1)*(i-1):(k+1)*i],u[.,1+(k+1)*(i-1):(k+1)*i]} = ols((e[.,i]~(x[2:t,1+(i-1)*k:i*k]-x[1:t-1,1+(i-1)*k:i*k])),q);   endif;
if est eq 2; {p[.,1+(k+1)*(i-1):(k+1)*i],u[.,1+(k+1)*(i-1):(k+1)*i]} = yulew((e[.,i]~(x[2:t,1+(i-1)*k:i*k]-x[1:t-1,1+(i-1)*k:i*k])),q); endif;
i = i + 1;
endo;

j  = 1;
nl = zeros(nb,1);
do while j <= nb;

i  = 1;
l  = 0;
h  = ceil(rndu(t,1)*(rows(u)));
do while i <= n;

e  = u[h,1+(i-1)*(k+1):i*(k+1)];
e  = e - meanc(e)';
e  = pfilter(e,p[.,1+(i-1)*(k+1):i*(k+1)],q);
xb = pfilter(e[.,2:k+1],eye(k),1);

if mod == 1; xb = xb~ones(t,1);             endif;
if mod == 2; xb = xb~ones(t,1)~seqa(1,1,t); endif;                                  

yb = xb*b[.,i] + e[.,1];	
l  = l + lm(yb,xb[.,1:k],mod,q);

i = i + 1;
endo;

{m,v} = mom(mod,k);
nl[j] = sqrt(n)*(l/n-m)./sqrt(v);
//"boostratp contry-stats""[";;"cv";;(j);;"]";;nl;
j = j + 1;
endo;

retp(sortc(nl,1));
endp;



proc(1) = pfilter(x,rho,p);
local t,n,y,yl,i,j;

t = rows(x);
n = cols(x);

i = p+1;
y = zeros(t+p,n);
do while i <= t+p;

j  = 1;
y[i,.] = x[i-p,.];
do while j <= p;
y[i,.] = y[i,.] + y[i-j,.]*rho[1+n*(j-1):n*j,.];
j  = j + 1;
endo;

i = i + 1;
endo;

retp(y[p+1:rows(y),.]);
endp;




proc(2) = mom(mod,k);
local mu,var;

if mod == 1; 
mu  = {0.11601 0.08464 0.06539 0.05295 0.04442}; 
var = {0.01151 0.00559 0.00266 0.00144 0.00086}; 
endif;

if mod == 2; 
mu  = {0.05530 0.04686 0.04063 0.03532 0.03133}; 
var = {0.00115 0.00078 0.00056 0.00036 0.00027}; 
endif;

retp(mu[k],var[k]);
endp;




proc (2) = fm(y,x,mod,q);
local w,b,dx,e,vl,v0,s,d,t,u,k;			

t = rows(y);
k = cols(x);

if mod == 1; w = x~ones(t,1);             endif;
if mod == 2; w = x~ones(t,1)~seqa(1,1,t); endif;                                  

b  = inv(w'w)*w'y;
u  = y - w*b;
dx = x[2:t,.]-x[1:t-1,.];
e  = (u[2:t]~dx);
vl = fejer(e,q);
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




proc (1) = lm(y,x,mod,q);
local w,b,dx,e,vl,v0,s,d,t,u,k;			

t = rows(y);
k = cols(x);

if mod == 1; w = x~ones(t,1);             endif;
if mod == 2; w = x~ones(t,1)~seqa(1,1,t); endif;                                  

b  = inv(w'w)*w'y;
u  = y - w*b;
dx = x[2:t,.]-x[1:t-1,.];
e  = (u[2:t]~dx);
vl = fejer(e,q);
v0 = (e'e)/rows(e);
s  = v0 + vl + vl';					
d  = v0 + vl';						
d  = t*((d[2:k+1,1]-d[2:k+1,2:k+1]*inv(s[2:k+1,2:k+1])*s[2:k+1,1])|zeros(cols(w)-k,1));

y  = y[2:t] - dx*inv(s[2:k+1,2:k+1])*s[2:k+1,1];									
x  = w[2:t,.];
b  = inv(x'x)*((x'y)-d);	
u  = y - x*b;
v0 = (u'u)/t;
vl = fejer(u,k);
s  = v0 + vl + vl';		

retp((sumc(abs(cumsumc(u)).^2))/((t^2)*s));
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

