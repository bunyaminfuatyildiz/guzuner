'Loop for: 
'Larsson, Lyhagen, Löthgren (2001). Likelihood-based cointegration tests in heterogeneous panels, Econometrics Journal (2001), volume 4, pp. 109–142.

'written by Saban Nazlioglu (snazlioglu@pau.edu.tr)
'05 August 2015

pageunstack(page=panelrank) crossid dateid @ lnfd lninc lnrem

scalar k=3		'number of endogenous
scalar n=9		'number of cross-sections
!mlag=2		'maksimum number of lags

matrix(6,6) mom

mom(1,1)=3.051
mom(2,1)=9.99
mom(3,1)=20.88
mom(4,1)=35.67
mom(5,1)=54.33
mom(6,1)=76.94

mom(1,2)=7.003
mom(2,2)=18.46
mom(3,2)=35.86
mom(4,2)=58.07
mom(5,2)=85.13
mom(6,2)=119.7

mom(1,3)=0.98
mom(2,3)=8.27
mom(3,3)=19.35
mom(4,3)=34.18
mom(5,3)=53.05
mom(6,3)=75.61

mom(1,4)=1.91
mom(2,4)=14.28
mom(3,4)=31.84
mom(4,4)=54.28
mom(5,4)=83.5
mom(6,4)=116.7

mom(1,5)=6.27
mom(2,5)=16.28
mom(3,5)=30.21
mom(4,5)=48.01
mom(5,5)=69.65
mom(6,5)=94.93

mom(1,6)=10.45
mom(2,6)=25.5
mom(3,6)=45.13
mom(4,6)=72.95
mom(5,6)=104.07
mom(6,6)=139.7

vector(k) LRm2	
vector(k) LRm3
vector(k) LRm4

for !i=1 to n
	vector(k+1) TRm2{!i}
	vector(k+1) TRm3{!i}
	vector(k+1) TRm4{!i}
	
		var var{!i}.ls 1 !mlag lnfd{!i} lninc{!i} lnrem{!i} 				'Unrestricted VAR(p) model
		freeze(tempa{!i}) var{!i}.laglen(!mlag, vname=lagv{!i})		'Creates optimal lag vector by criterion
		!optp=lagv{!i}(4,1) 											'Optimal lag by Schwarz
		var varc{!i}.ls 1 !optp lnfd{!i} lninc{!i} lnrem{!i}				'Unrestricted VAR(optp) model

		freeze(M2{!i}) varc{!i}.coint(b,!optp,save=test2{!i})			'Model 2
		freeze(M3{!i}) varc{!i}.coint(c,!optp,save=test3{!i})			'Model 3
		freeze(M4{!i}) varc{!i}.coint(d,!optp,save=test4{!i})			'Model 4
		
		'Store Trace statistics
		TRm2{!i}=@columnextract(test2{!i},3)							
		TRm3{!i}=@columnextract(test3{!i},3)
		TRm4{!i}=@columnextract(test4{!i},3)
next




for !r=1 to k

scalar ind=1
	for !i=1 to n
		vector(n) TRCm2{!r}
		vector(n) TRCm3{!r}
		vector(n) TRCm4{!r}

		TRCm2{!r}({!i})=TRm2{!i}(!r,1)
		TRCm3{!r}({!i})=TRm3{!i}(!r,1)
		TRCm4{!r}({!i})=TRm4{!i}(!r,1)

		LRm2({!r})=@mean(TRCm2{!r})
		LRm3({!r})=@mean(TRCm3{!r})
		LRm4({!r})=@mean(TRCm4{!r})

		'Set results into a table
			table(n*n,3*k+1) results
			results.setformat(@all) f.3
			results(1,2)        ="Model 2"
			results(1,2+k)    ="Model 3"
			results(1,2+2*k) ="Model 4"

			results(2,1)        	  ="Country"
			results(2,1+!r)       ="rank("+@str(!r-1)+")"
			results(2,1+!r+k)   ="rank("+@str(!r-1)+")"
			results(2,1+!r+2*k)="rank("+@str(!r-1)+")"

			setline(results,3)
			
			results.setformat(a) f.0
			results(ind+4,1)=!i

			results(ind+4,1+!r)       =TRCm2{!r}(ind)
			results(ind+4,1+!r+k)   =TRCm3{!r}(ind)
			results(ind+4,1+!r+2*k)=TRCm4{!r}(ind)

			setline(results,n+5)

			results(n+6,1)       	="LR_mean"
			results(n+6,1+!r)       =LRm2(!r)
			results(n+6,1+!r+k)   =LRm3(!r)
			results(n+6,1+!r+2*k)=LRm4(!r)

			results(n+7,1)="Mean"
			results(n+7,1+k-!r+1)       =mom(!r,1)
			results(n+7,1+k-!r+1+k)   =mom(!r,3)
			results(n+7,1+k-!r+1+2*k)=mom(!r,5)

			results(n+8,1)="Variance"
			results(n+8,1+k-!r+1)       =mom(!r,2)
			results(n+8,1+k-!r+1+k)   =mom(!r,4)
			results(n+8,1+k-!r+1+2*k)=mom(!r,6)
			
			results(n+9,1)="Panel-stat."
			results(n+9,1+!r)    =@sqrt(n)*(LRm2(!r)-mom(k+1-!r,1))/@sqrt(mom(k+1-!r,2))
			results(n+9,1+!r+k)  =@sqrt(n)*(LRm3(!r)-mom(k+1-!r,3))/@sqrt(mom(k+1-!r,4))
			results(n+9,1+!r+2*k)=@sqrt(n)*(LRm4(!r)-mom(k+1-!r,5))/@sqrt(mom(k+1-!r,6))
			
			setline(results,n+10)			

scalar ind=ind+1
next
next
show results

delete lagv* temp* test* var* m2* m3* m4*
delete ind k n trcm* trm* lrm* mom


