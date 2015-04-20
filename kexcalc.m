function [results] = kexcalc(IntensityAA,IntensityBB,IntensityAB,IntensityBA,time,guess)
%KEXCALC calculates the exchange rate for a ZZ-exchange NMR experiment.

%Default guess values
default.kex = 2;
default.pA = 0.3;
default.R1A = 1;
default.R1B = 1;

%Use guesses from input structure; use default otherwise
if exist('guess','var') == 0
    guess = default;
end
if isfield(guess,'kex') == 0
    guess.kex = default.kex;
end
if isfield(guess,'pA') == 0
    guess.pA = default.pA;
end
if isfield(guess,'R1A') == 0
    guess.R1A = default.R1A;
end
if isfield(guess,'R1B') == 0
    guess.R1B = default.R1B;
end
guesses = [guess.kex;guess.pA;guess.R1A;guess.R1B];

%Subtract initial values from cross peaks and normalize intensities
ABsub = IntensityAB - IntensityAB(1);
BAsub = IntensityBA - IntensityBA(1);
total_intensity = IntensityAA(1) + IntensityBB(1);
normAA = IntensityAA/total_intensity;
normBB = IntensityBB/total_intensity;
normAB = ABsub/total_intensity;
normBA = BAsub/total_intensity;

%Make sure intensities are column vectors
[AAcol,AArow] = size(normAA);
if AArow > AAcol
    normAA = normAA';
end
[BBcol,BBrow] = size(normBB);
if BBrow > BBcol
    normBB = normBB';
end
[ABcol,ABrow] = size(normAB);
if ABrow > ABcol
    normAB = normAB';
end
[BAcol,BArow] = size(normBA);
if BArow > BAcol
    normBA = normBA';
end

%Call fminsearch with zzexchange to calculate the best fit values
[fitparameters] = fminsearch(@(x) zzexchange(x,normAA,normBB,normAB,normBA,time),guesses);
kex = fitparameters(1);
pA = fitparameters(2);
pB = 1 - pA;
R1A = fitparameters(3);
R1B = fitparameters(4);
[residual,IAA,IBB,IAB,IBA] = zzexchange(fitparameters,normAA,normBB,normAB,normBA,time);
results.kex = kex;
results.pA = pA;
results.pB = pB;
results.R1A = R1A;
results.R1B = R1B;
results.datafit = [normAA normBB normAB normBA IAA IBB IAB IBA];
results.fitIAA = IAA;
results.fitIBB = IBB;
results.fitIAB = IAB;
results.fitIBA = IBA;
results.dataIAA = normAA;
results.dataIBB = normBB;
results.dataIAB = normAB;
results.dataIBA = normBA;
results.residual = residual;

%Display results!
figure
hold on
plot(time,normAA,'d',time,normBB,'d',time,normAB,'d',time,normBA,'d')
plot(time,IAA,time,IBB,time,IAB,time,IBA)
ylabel('Intensity (au)')
xlabel('Time (s)')
title('k_{ex} ZZ-exchange fitting')
legend('I_{AA}','I_{BB}','I_{AB}','I_{BA}')
box on
fprintf('\n~~~~~~~~Fitting Results~~~~~~~~\nkex: %f\npA: %f\npB: %f\nR1A: %f\nR1B: %f\nresidual: %f\n',kex,pA,pB,R1A,R1B,residual)
end

function [residual,IAA,IBB,IAB,IBA] = zzexchange(guesses,IntensityAA,IntensityBB,IntensityAB,IntensityBA,time)
%Reconstruct input
kex = guesses(1);
pA = guesses(2);
pB = 1 - pA;
R1A = guesses(3);
R1B = guesses(4);
number_times = length(time);

%Calculate intensities
lambdaplus = (1/2)*(R1A+R1B+kex+((R1A-R1B+kex*(pB-pA))^2+4*pA*pB*kex^2));
lambdaminus = (1/2)*(R1A+R1B+kex-((R1A-R1B+kex*(pB-pA))^2+4*pA*pB*kex^2));

aAA = zeros(number_times,1);
aBB = zeros(number_times,1);
aAB = zeros(number_times,1);
aBA = zeros(number_times,1);

for n = 1:1:number_times
    aAA(n) = (1/2)*((1-(R1A-R1B+kex*(pB-pA))/(lambdaplus-lambdaminus))*exp(-lambdaminus*time(n))+(1+(R1A-R1B+kex*(pB-pA))/(lambdaplus-lambdaminus))*exp(-lambdaplus*time(n)));
    aBB(n) = (1/2)*((1+(R1A-R1B+kex*(pB-pA))/(lambdaplus-lambdaminus))*exp(-lambdaminus*time(n))+(1-(R1A-R1B+kex*(pB-pA))/(lambdaplus-lambdaminus))*exp(-lambdaplus*time(n)));
    aAB(n) = ((kex*pA)/(lambdaplus-lambdaminus))*((exp(-lambdaminus*time(n))-exp(-lambdaplus*time(n))));
    aBA(n) = ((kex*pB)/(lambdaplus-lambdaminus))*((exp(-lambdaminus*time(n))-exp(-lambdaplus*time(n))));
end

IAA = pA*aAA;
IBB = pB*aBB;
IAB = pB*aAB;
IBA = pA*aBA;

%Compare fit with data
deltaAA = IntensityAA - IAA;
deltaBB = IntensityBB - IBB;
deltaAB = IntensityAB - IAB;
deltaBA = IntensityBA - IBA;
residualAA = sum(deltaAA.^2);
residualBB = sum(deltaBB.^2);
residualAB = sum(deltaAB.^2);
residualBA = sum(deltaBA.^2);
residual = (residualAA + residualBB + residualAB +residualBA);
end
