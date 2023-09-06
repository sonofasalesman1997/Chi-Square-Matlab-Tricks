close all
clear
clc

totaldata = [1992,19; 1993,21; 1994,21; 1995,24; 1996,22; 1997,19; 1998,20; 1999,24; 2000,25; 2001,22; 2002,20; 2003,21; 2004,22; 2005,21; 2006,21; 2007,25; 2008,24; 2009,27; 2010,25; 2011,33; 2012,34; 2013,37; 2014,37; 2015,31; 2016,30; 2017,31; 2018,34; 2019,36; 2020,36; 2021,38; 2022,37; 2023,39];
data = [2002,20; 2003,21; 2004,22; 2005,21; 2006,21; 2007,25; 2008,24; 2009,27; 2010,25; 2011,33; 2012,34; 2013,37; 2014,37; 2015,31; 2016,30; 2017,31; 2018,34; 2019,36; 2020,36; 2021,38; 2022,37; 2023,39];
year = data(:,1);   % for graphical purposes
x = data(:,1)-1992; % for analysis
y = data(:,2);

dx = .5*ones(size(x));
dxneg = .5*ones(size(x));
dxpos = .5*ones(size(x));

dy = ones(size(y));
dyneg = ones(size(y));
dypos = ones(size(y));
errtot = sqrt(dy.^2+dx.^2);

figure(1)
errorbar(year,y,dyneg,dypos,dxneg,dxpos,'o')
title("State Trifectas")
xlabel("Year of our Lord")
ylabel("Total Trifectas")
grid on;grid minor
get(gcf)
set(gcf,'Position',[0 600 400 300])

% The model chosen for this: linear y = mx + b
mRange = -2:.01:3;  % range of your guesses for m
bRange = 8:.01:13; % range of your guesses for b
m = 1;  % first guess for m
b = 10; % first guess for b
dof = length(y)-2; %number of elements - number of parameters
ChiSqrGuess = sum(((y-(m.*x+b))./(errtot)).^2);
disp(ChiSqrGuess)

% Use the guess value to approximate for the closest minimum ChiSqr
for i = 1:length(bRange)
    ChiSqr_b(i) = sum(((y-(m.*x+bRange(i)))./(errtot)).^2);
end
for i = 1:length(mRange)
    ChiSqr_m(i) = sum(((y-(mRange(i).*x+b))./(errtot)).^2);
end

figure(2)
plot(bRange,ChiSqr_b,'.')
grid on;grid minor
title("Minimizing b")
xlabel("b Range")
ylabel("\chi^2_b")
get(gcf)
set(gcf,'Position',[400 600 400 300]);

figure(3)
plot(mRange,ChiSqr_m,'.')
grid on;grid minor
title("Minimizing m")
xlabel("m Range")
ylabel("\chi^2_m")
get(gcf)
set(gcf,'Position',[800 600 400 300])

[ChiSqrMin_b,index_b] = min(ChiSqr_b);
bguess = bRange(index_b);
disp(bguess)
[ChiSqrMin_m,index_m] = min(ChiSqr_m);
mguess = mRange(index_m);
disp(mguess)

% Always remember that the parameter value of the minimum chi square plus one must be the range of the uncertainty of each parameter
% You could also find uncertainty through the error propagation method and compare these
% Sometimes, must be calculated through uncertainty propagation!
% ex: f=f(x,y) => df = sqrt{ (δx * ∂f/∂x)^2 + (δy * ∂f/∂y)^2 }

figure(4)
plot(x,mguess.*x+bguess,'b.')
hold on
errorbar(x,y,dyneg,dypos,dxneg,dxpos,'o')
hold off
grid on;grid minor
title("Guess Values vs Data")
xlabel("Years since 1992")
ylabel("Total Trifectas")
get(gcf)
set(gcf,'Position',[1200 600 400 300])

figure(5)
plot(bRange,ChiSqr_b)
hold on
plot([bguess bguess 2*bguess],[ChiSqrMin_b ChiSqrMin_b+1 ChiSqrMin_b+1],'k')
xlim([8 11])
ylim([100 300])
grid on;grid minor
title("b Range zoomed in")
xlabel("b Range")
ylabel("\chi^2_b")
get(gcf)
set(gcf,'Position',[0 200 400 300])

figure(6)
plot(mRange,ChiSqr_m)
hold on
plot([mguess mguess 2*mguess],[ChiSqrMin_m ChiSqrMin_m+1 ChiSqrMin_m+1],'k')
xlim([0.9 1])
ylim([130 150])
grid on;grid minor
title("m Range zoomed in")
xlabel("m Range")
ylabel("\chi^2_m")
get(gcf)
set(gcf,'Position',[400 200 400 300])

% Create a matrix where each row/column correspond to differing m and b values respectively
for i = 1:length(ChiSqr_b)
    for j = 1:length(ChiSqr_m)
        ChiSqr_total(i,j) =  sum(((y-(mRange(j).*x+bRange(i))).^2)./(errtot.^2));
    end
end

% WARNING: this plot can take a long time to display if too many points
% NOTE: Only works if the range sizes are the same, figure only for display purposes
% NOTE: if more than 2 parameters, impossible to graph in 4+ dimensions
figure(7)
surf(bRange,mRange,ChiSqr_total)
grid on;grid minor
title("3D Topographic \chi^2 View")
xlabel("b Range")
ylabel("m Range")
zlabel("\chi^2 Values")
get(gcf)
set(gcf,'Position',[800 200 400 300])

% Now to find the TRUE overall minimum ChiSqr value from both parameters
ChiSqrMin_total = min(ChiSqr_total(:));
disp(ChiSqrMin_total)
[row, col] = find(ChiSqr_total == ChiSqrMin_total);
% The total reduced ChiSqr value
ChiSqrRed_total = ChiSqrMin_total/dof;
disp(ChiSqrRed_total)
% The true end goal: the parameters that match with lowest ChiSqrMin
truem = (mRange(col));
disp(truem)
trueb = (bRange(row));
disp(trueb)
