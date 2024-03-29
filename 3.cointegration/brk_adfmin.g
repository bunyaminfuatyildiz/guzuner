new;

/*
** adfmin.prg was written by Joakim Westerlund, Department of Economics, Lund University.
** Contact info:  joakim.westerlund@nek.lu.se
** 
** The code can be used freely as long as proper reference is given. No performance
** guarantee is made. Bug reports are welcome.
*/

t       = 28;
n       = 9;
k       = 2;                /* k = 5 is max */                    
p       = int(t^(1/3));     /* bandwidth */
t0      = 0.1;              /* trimming */
mod     = 2;                /* 1 = constant break (level shift) 
                               2 = constant break with trend 
                                   (level shift with trend)*/

load y[t,n]     = lnfd.txt;    /* t x n matrix */
load x[t,n*k]   = west_x.txt;    /* t x kn matrix */


i     = 1;
adfn1 = 0;
adfn2 = 0;
adfn3 = 0;
adfn4 = 0;
brn1  = zeros(1,n);
brn2  = zeros(1,n);
brn3  = zeros(1,n);
brn4  = zeros(1,n);
t1    = round(t0*t);
t2    = round((1-t0)*t); 
do while i <= n;

j    = t1;
adf1 = zeros(t,2);
adf2 = zeros(t,2);
do while j <= t2;
 
adf1[j,.] = adfs(y[.,i],x[.,1+k*(i-1)],j,p,mod);
adf2[j,.] = adfp(y[.,i],x[.,1+k*(i-1)],j,p,mod);

j = j + 1;
endo;

adfn1 = adfn1 + minc(adf1[t1:t2,1]);
adfn2 = adfn2 + minc(adf1[t1:t2,2]);
adfn3 = adfn3 + minc(adf2[t1:t2,1]);
adfn4 = adfn4 + minc(adf2[t1:t2,2]);

brn1[i] = minpr(adf1[t1:t2,1]) + t1 - 1;         
brn2[i] = minpr(adf1[t1:t2,2]) + t1 - 1;        
brn3[i] = minpr(adf2[t1:t2,1]) + t1 - 1;            
brn4[i] = minpr(adf2[t1:t2,2]) + t1 - 1;           
i = i + 1;
endo;

   
{mu,var} = mom(mod,k);
adfn1    = sqrt(n)*(adfn1/n-mu[1])./sqrt(var[1]);
adfn2    = sqrt(n)*(adfn2/n-mu[2])./sqrt(var[2]);
adfn3    = sqrt(n)*(adfn3/n-mu[1])./sqrt(var[1]);
adfn4    = sqrt(n)*(adfn4/n-mu[2])./sqrt(var[2]);



/* printing options */
format/m1/rd 8,4;

print " ";
print " sz_t   = ";; adfn1;;    print " asymp p-val  = ";; 1-cdfnc(adfn1); 
print " sz_rho = ";; adfn2;;    print " asymp p-val  = ";; 1-cdfnc(adfn2); 
print " pz_t   = ";; adfn3;;    print " asymp p-val  = ";; 1-cdfnc(adfn3); 
print " pz_rho = ";; adfn4;;    print " asymp p-val  = ";; 1-cdfnc(adfn4); 

format/m1/rd 6,0;
"";
"";
"Break dates";
"counrty  SZ_t  SZ_rho  PZ_t  PZ_rho";;
seqa(1,1,n)~brn1'~brn2'~brn3'~brn4'; 






/* procs */


proc (1) = adfs(y,x,br,p,mod);
local t,v,m,w,wl,dw,b,e,s,io,za,zt;

t = rows(y);

if br == 0;
if mod == 1; v = ones(t,1)~x;               endif;
if mod == 2; v = ones(t,1)~seqa(1,1,t)~x;   endif;
else;
if mod == 1; v = ones(t,1)~(zeros(br,1)|ones(t-br,1))~x;                endif;
if mod == 2; v = ones(t,1)~(zeros(br,1)|ones(t-br,1))~seqa(1,1,t)~x;    endif;
endif;

m   = eye(t)-v*inv(v'v)*v';
w   = m*y;                        /* eq error */
wl  = w[1:t-1];
dw  = w[2:t] - wl;

b   = inv(wl'wl)*(wl'dw);
e   = dw - wl*b;
s  	= (e'e)/rows(e);
if p == 0;  io = 0;     	    
else;       io = fejer(e,p);   endif;
s   = s + io + io';	
				
zt = inv(sqrt(wl'wl*s))*(wl'dw - t*io);
b  = inv(wl'wl)*(wl'dw - t*io);		
za = t*b;
 
retp(zt~za);
endp;




proc (1) = adfp(y,x,br,p,mod);
local t,v,m,w,wl,dw,b,e,s,za,zt;

t = rows(y);

if br == 0;
if mod == 1; v = ones(t,1)~x;               endif;
if mod == 2; v = ones(t,1)~seqa(1,1,t)~x;   endif;
else;
if mod == 1; v = ones(t,1)~(zeros(br,1)|ones(t-br,1))~x;                endif;
if mod == 2; v = ones(t,1)~(zeros(br,1)|ones(t-br,1))~seqa(1,1,t)~x;    endif;
endif;

m   = eye(t)-v*inv(v'v)*v';
w   = m*y;                        
wl  = w[p:t-1];
dw  = w[p+1:t] - wl;

if p == 1;
m   = eye(t-1);
else;
v   = lagp(w[1:t-1],p-1);
m   = eye(rows(v))-v*inv(v'v)*v';   
endif;

b   = inv(wl'm*wl)*(wl'm*dw);
e   = m*dw - m*wl*b;
s   = (e'e)/rows(e);
zt  = sqrt(wl'm*wl/s)*b;
za  = t*b;

retp(zt~za);
endp;



proc(2) = mom(mod,k);
local mu1,mu2,var1,var2;

mu1 = { -3.49890	-3.85190	-4.18972	-4.51141	-4.76502 ,
        -3.95088	-4.24887	-4.55330	-4.75105	-4.99327 };

mu2 = { -25.00168	-30.29004	-35.68431	-41.31450	-46.09029 ,
        -31.54225	-36.59367	-41.98604	-46.05167	-50.97382 };

var1 = { 0.39038    0.36498	  0.36336	0.35969	  0.32972 ,
         0.33317	0.34686	  0.32819	0.34866	  0.33897 };

var2 = { 72.38729	83.39962	99.14725	115.23283	117.25124 ,
         78.17861	96.40939	104.02695	125.02231	132.19716 };

retp((mu1[mod,k]|mu2[mod,k]),(var1[mod,k]|var2[mod,k]));
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


proc lrvar(u,p);
local lr,io,s;
io = fejer(u,p);
s  = (u'u)/rows(u);
lr = s + io + io';
retp(lr);
endp;



proc (1) = lagp(x,p);
local t,i,xl,dx,y,rio;

t = rows(x);

if p == 0;
dx = { };
retp(dx);
endif;

i = 0; 
xl = zeros(t-p,p+1);
do while i <= p;
xl[.,i+1] = x[p+1-i:t-i,.];
i = i + 1; 
endo;

i = 1; 
dx = zeros(t-p,p);
do while i <= p;
dx[.,i] = xl[.,i] - xl[.,i+1];
i = i + 1; 
endo;

retp(dx); 
endp;



proc lagn(x,n);
local y;
y = shiftr(x',n,(miss(0, 0))');
retp(y');
endp;


proc (1) = minpr(x);
local d,m,i,minx;

i = 1;
m = 1;
d = x[1];
do while i <= rows(x);
if x[i] < d;  
d = x[i];  
endif;
i = i + 1;
endo;

i = 2;
do while i <= rows(x);
if d == x[i]; 
m = i;  
goto stops; 
endif;
i = i + 1;
endo;

stops:

retp(m);
endp;
