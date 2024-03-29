/*
** cusum.prg was written by Joakim Westerlund, Department of Economics, Lund University.
** Contact info:  joakim.westerlund@nek.lu.se
** 
** The code can be used freely as long as proper reference is given. No performance
** guarantee is made. Bug reports are welcome.
*/
new;
//cls;

k       = 2;               /* regressors, k = 4 is max */
n       = 9;              /* number of cross-section */  
t       = 28;              /* time period*/ 
mod     = 2;               /* 0 = nothing, 1 = const, 2 = const and trend */	        
fm      = 1;               /* 0 = dols, 1 = fmols */
p       = 2;               /* lags and leads in dols */
q       = int(t^(1/3));    /* bandwidth */

load y[t,n]  =lnfd.txt;
load x[t,n*k]=west_x.txt;


/* test begins */

i   = 1;
cus = 0;
do until i > n;

if fm == 0; u = dols(y[.,i],x[.,1+(i-1)*k:i*k],mod,p);  endif;
if fm == 1; u = fmols(y[.,i],x[.,1+(i-1)*k:i*k],mod,q); endif;

cus = cus + cusum(u,q);

i = i + 1;
endo;

{mu,var} = mom(mod,k);
cus      = sqrt(n)*(cus/n-mu)./sqrt(var);


/* printing options */
format/m1/rd 8,3;

print " ";
print " cusum = ";; cus;
print " p-value  = ";; cdfnc(cus);






/* procedures */

proc(1) = fmols(y,x,mod,q);
local w,xd,b,dx,e,io,s,sig,del,s21,s12,del21;
local ys,xk,ixx,t,u,d,k;			

t = rows(y);
k = cols(x);

if mod == 0;
w   = x; 
xd  = x; 
else;
d   = detc(mod,t);
w   = x~d; 
xd  = x-d*inv(d'd)*d'x; 
endif;

b 	= inv(w'w)*w'y;
u   = y - w*b;
dx 	= xd[2:t,.]-xd[1:t-1,.];
e 	= u[2:t]~dx;

io 	= fejer(e,q);
s  	= (e'e)/rows(e);
sig = s + io + io';					
del = s + io';						
s21	= inv(sig[2:1+k,2:1+k])*sig[2:1+k,1];
s12 = sig[1,1] - sig[1,2:1+k]*s21;								

if mod == 0;
del21 = t*(del[2:1+k,1]-del[2:1+k,2:1+k]*s21); 
else;	
del21 = t*(del[2:1+k,1]-del[2:1+k,2:1+k]*s21)|zeros(cols(d),1);
endif;

ys  = y[2:t] - dx*s21;									
xk  = w[2:t,.];
ixx = inv(xk'xk);					
b 	= ixx*((xk'ys)-del21);	
//print "bfmols" b;
u   = ys - xk*b;

retp(u);
endp;



proc(1) = dols(y,x,mod,p);
local i,xi,yi,t,n,dx,j,z,b,u,k;

t = rows(y);
k = cols(x);

j  = 1;	
dx = zeros(t-2*p-1,k*(2*p+1)); 		
do while j <= 2*p+1;				
dx[.,1+(j-1)*k:j*k] = x[2*p+3-j:t+1-j,.] - x[2*p+2-j:t-j,.];
j = j + 1;
endo;

xi = x[1:t-2*p-1,.];
yi = y[1:t-2*p-1,.];
z  = eye(rows(dx)) - dx*inv(dx'dx)*dx';

if mod == 0;
xi = xi;
else;
xi = detc(mod,t-2*p-1)~xi;
endif;

xi = z*xi;		
yi = z*yi;
b  = invpd(xi'xi)*xi'yi;
//print "bdols" b;
u  = yi - xi*b;

retp(u);
endp;



proc(1) = detc(mod,t);
local du,dt,z;     

if mod == 1; z = ones(t,1);             endif;
if mod == 2; z = ones(t,1)~seqa(1,1,t); endif;
    
retp(z);
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



proc (1) = cusum(u,q);
local t,io,s,cus;

t  = rows(u);
io = fejer(u,q);
s  = (u'u)/rows(u);
s  = s + io + io';

cus = maxc(1/(sqrt(t*s)).*abs(cumsumc(u)));

retp(cus);
endp;



proc(2) = mom(mod,k);
local mu,var;

mu = {  1.063552 0.930258 0.843500 0.771490 0.727086 ,
        0.751281 0.676628 0.619667 0.573433 0.539597 ,
        0.584626 0.553611 0.521820 0.499409 0.477681 };

var = { 0.183683 0.126593 0.093796 0.070758 0.059383 ,
        0.050239 0.037289 0.028520 0.021072 0.017286 ,
        0.018457 0.015716 0.013679 0.012123 0.010621 };

retp(mu[mod+1,k],var[mod+1,k]);
endp;

