% Function for calculating multiscale entropy
% input: signal
% m: match point(s)
% r: matching tolerance
% factor: number of scale factor
% sampenc is available at http://people.ece.cornell.edu/land/PROJECTS/Complexity/sampenc.m

function e = fuzzymsentropy(input,m,mf,rn,local,tau,factor)

y=input;
y=y-mean(y);
y=y/std(y);

for i=1:factor
   s=coarsegraining(y,i);
   sampe=FuzEn_MFs(s, m, mf, rn, local,tau);
   e(i)=sampe;   
end
e=e';