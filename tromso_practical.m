%% PART I
format long g; % display 15 digits after decimal point
% Load the daily co2 data:

fid1 = fopen('co2_mlo_surface-insitu_1_ccgg_DailyData.txt','r');
dum = textscan(fid1,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s\n','Headerlines',142);
clear fid1
t = datenum(dum{2},dum{3},dum{4});
x = dum{8};
N = length(t);
dt = t(2)-t(1);

%% 
% The first 136 values are missing values (-999.99) which we remove, and 
% we interpolate other missing data in the middle of the time series:

t(1:136) = [];
t = t - t(1);
x(1:136) = [];
q = x == -999.99;
x(q) = NaN;
x = fillbad(x);
N = length(t);
clear foo dum q

%% 
% We can now plot the data

f1 = figure;
hx = plot(t,x);
axis tight;
ylabel('ppm');
xlabel('time');

%% 
% 
% What are the basic statistics?
mx = mean(x);
vx = var(x);
disp([mx vx.^0.5]);

%% 
% The function hlines from jLab quickly adds horizontal lines to a plot:

hh = hlines([mx mx+vx.^0.5 mx-vx.^0.5]);

%%
% Let's plot the data histogram
% NOTE FOR OCTAVE USERS: the histogram function is not implemented in
% octave yet; use "hist" instead for approximate same results
f2 = figure;
h = histogram(x);
ylabel('count');
xlabel('CO_2 (ppm)');

% The function vlines from jLab simply adds horizontal lines to a plot:
vlines([mx mx+vx.^0.5 mx-vx.^0.5]);

%%
% Let's estimate the trend as a linear function of time as y1 = a*t + b
% For this we use Matlab polynomial estimator p_n t^n + p_{n-1} t^{n-1} +
% ... + p_1 t^1 + p_0
p1 = polyfit(t,x,1);
tr1 = polyval(p1,t);

% add linear trend to plot
figure(f1);
hold on
delete(hh);
h1 = plot(t,tr1,'linewidth',2);

%%
% I think the long-term trend may be modeled by a quadratic function of
% time: tr2 = a*t^2 + b*t + c
p2 = polyfit(t,x,2);
tr2 = polyval(p2,t);
% add to plot
figure(f1);
h2 = plot(t,tr2,'linewidth',2);
% add a legend to plot
legend([hx h1 h2],{'x','linear trend','quadratic trend'},'location','best');
title('CO2');

%%
% Let's plot the residual from the linear trend
f3 = figure;
hold on
h1 = plot(t,x-tr1);
h2 = plot(t,x-tr2);
legend([h1 h2],{'x - tr1','x - tr2'},'location','best');
title('CO2 minus trends');
ylabel('ppm');
xlabel('time');
axis tight

% new stats
x1 = x - tr1;
x2 = x - tr2;
m1 = mean(x1);
m2 = mean(x2);
v1 = var(x1);
v2 = var(x2);
s1 = sum(x1.^3)/N/v1^(3/2);
s2 = sum(x2.^3)/N/v2^(3/2);
k1 = sum(x1.^4)/N/v1^(4/2);
k2 = sum(x2.^4)/N/v2^(4/2);

disp([m2 v1 s1 k1 ; m2 v2 s2 k2]);
%note that the trend fitting has removed the mean value from x

%%
% what is the new histogram like?
figure
hold on
h1 = histogram(x-tr1);
h2 = histogram(x-tr2);
legend([h1 h2],{'x - tr1','x - tr2'},'location','best');
xlabel('CO_2 anomalies (ppm)');
ylabel('count');

%% 
% I think it is time to calculate spectrum estimates of all these series
% let's be naive and calculate the periodogram
dum = fft([x-mx x1-m1 x2-m2],[],1);
% see lecture on how to obtain periodogram estimate from Matlab fft command
spa = abs(dum).^2*dt/N;

% did we get this right? Check Parseval's theorem
disp(sum(spa,1)*(1/N*dt));
disp([vx v1 v2]);

% How does it look like?
figure
h = plot(spa);
legend(h,{'x','x1','x2'},'location','best');
ylog;
% this is not very useful here? What is the x-axis?

%%
% the jLab function fourier.m computes the Fourier frequencies for you, in radian per
% cycle
fp = fourier(dt,N); % that's N/2+1 frequencies
% check we obtain the expected frequencies
disp(fp([1 2 end])/(2*pi));
disp(1/(N*dt)); % this is the Rayleigh frequency

% As we saw in the lecture, for n > N/2, frequencies are negative, and the
% Fourier components contain redundant information, so we just need to
% consider the first N/2 components, and scale by 2

spa = 2*spa(1:N/2+1,:);
% did we get this right? Check Parseval's theorem again
disp(sum(spa,1)*fp(2)/(2*pi));
disp(var([x x-tr1 x-tr2]));

%%
% now we can plot our estimate
figure
hs = plot(fp/(2*pi),spa);
xlog;% this is a jLab function, a short cut for set(gca,'yscale','log')
ylog;
legend(hs,{'S_x','S_{x1}','S_{x2}'});
title('periodogram');
xlabel('Frequency (cpd)');ylabel('PSD (ppm^2 cpd^{-1})')
hv = vlines((365.25*[1 0.5 4/12 1/12]).^-1);

%%
% note that you can compute the periodogram quickly using jLab mspec
% function:
[f,spb] = mspec(dt,[x-mx x1-m1 x2-m2],[]);

figure
plot(f/(2*pi),spa(:,1))
ylog;
hold on
plot(f/(2*pi),spb(:,1),'--');
xlog
title('Periodogram');

%% let's now calculate a multitaper estimate
% define the multitaper parameters
NW = 3; % this is the time-frequency bandwidth; it sets the bias reduction
        % try to vary this from 2 to many in 0.5 increments
K = 2*NW-1; % this is a number of tapers, or windows, recommended based on NW
% calculate the "Slepian" tapers
psi = sleptap(N,NW,K);

% calculate the multitaper
% note that the mspec function automatically remove the mean of the data
[f,sm] = mspec(dt,[x-mx x1-m1 x2-m2],psi);

% plot the result
f4 = figure;
hold on
hs = plot(f/(2*pi),sm);
legend(hs,{'S_x','S_{x1}','S_{x2}'});
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');
ylog;
xlog;
axis tight
xtick(10.^(-5:0));
ytick(10.^(-4:2:8));
axis tight;
hv = vlines((365.25*[40 20 10 2 1 0.5 1/3 1/12]).^-1);

%%
% let's plot again x1 and x2 with confidence intervals

% Confidence intervals for spectral estimates; 
% see Bendat and Piersol book: as an example
alpha = 0.05; % significance level is 100*(1-alpha)
% need to provide values for inverse cumulative chi square function
% for K = 5 and alpha = 0.05, chi2inv(1-alpha/2,2*K) = 20.48
% and chi2inv(alpha/2,2*K) = 3.25
ciK = [2*K/chi2inv(1-alpha/2,2*K) 2*K/chi2inv(alpha/2,2*K)];

cc = lines(3); % this defines colors

f5 = figure;
hold on
hs2 = patch([f(2:end) ; flipud(f(2:end))]/(2*pi),[sm(2:end,2)*ciK(1) ; flipud(sm(2:end,2)*ciK(2))],'w');
set(hs2,'edgecolor','none','facecolor',whiten(cc(2,:)));
hs3 = patch([f(2:end) ; flipud(f(2:end))]/(2*pi),[sm(2:end,3)*ciK(1) ; flipud(sm(2:end,3)*ciK(2))],'w');
set(hs3,'edgecolor','none','facecolor',whiten(cc(3,:)));
xlog;ylog;
hs = plot(f/(2*pi),sm(:,2:3));
set(hs(1),'color',cc(2,:),'linewidth',2);
set(hs(2),'color',cc(3,:),'linewidth',2);
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');
xtick(10.^(-5:0));
hv = vlines((365.25*[40 20 10 2 1 0.5 1/3 1/12]).^-1);

% as expected, the two spectra differ at time scales longer than 
% please vary the number of tapers K and level alpha of confidence
% intervals


%%
% Can we estimate the annual cycle from the residuals of trend?

% let's try to fit sinusoidal functions
% I wrote a simple function to estimate by least squares a*cos(2 \pi f t + \phi)
[a,phi,r2a] = xy2cos(t,x2,1./[365.25]);
% r2a is the residual from the fit; subtract r2a from original time series to
% get estimated sinusoidal cycle
x2a = x2 - r2a;

figure
hold on
h = plot(t,[x2 x2a r2a]);
legend(h,{'x2','x2a','residual'});

[f,s] = mspec(dt,[x2 x2a r2a],psi);
figure;
hs = plot(f/(2*pi),s);
ylog;
xlog;
legend(hs,{'x2','x2a','residual'});

%%
[a,phi,r2s] = xy2cos(t,x2,1./[365.25./[1 2 3]]);
% get estimated seasonal cycle
x2s = x2 - r2s;

figure
hold on
plot(t,[x2 x2s r2s]);

[f,s] = mspec(dt,[x2 x2s r2s],psi);
figure;
hs = plot(f/(2*pi),s);
ylog;
xlog;
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');

% we may see later other methods to extract the seasonal signal

%% Filtering
% Let's estimate the low-frequency variability by applying a smoothing
% window. First, a theoretical interlude, see extra slides

% define several windows, length of 3 years
M = round(3*365);
w0 = ones(M,1);
w0 = w0/sum(w0);% boxcar window normalizing to one
w1 = jhanning(M);
w1 = w1/sum(w1);% Hanning window normalizing to one

figure;
h = plot(1:M,[w0 w1]);
legend(h,{'boxcar','Hanning'},'location','best');

% what are the Fourier transform of these windows?
dum = fft([w0 w1],N,1);
sw = 2*abs(dum(1:N/2+1,:)).^2*dt/N;

figure;
h = plot(f/(2*pi),sw);
xlog;ylog;
legend(h,{'boxcar','Hanning'},'location','best');
hv = vlines((M).^-1);
ylim(10.^[-15 -3]);

%%
% vfilt is a jLab function that conducts the convolution operation
x2f0 = vfilt(x2,w0,'mirror');
x2f1 = vfilt(x2,w1,'mirror');

figure;
hold on
h1 = plot(t,x2);
h2 = plot(t,[x2f0 x2f1],'linewidth',2);
legend([h1;h2],{'r2s','boxcar','hanning'},'location','best');

% let's examine what happens in the spectral domain
[f,s] = mspec(dt,[x2 x2f0 x2f1],psi);
figure;
hold on
hs = plot(f/(2*pi),s);
legend(hs,{'x2','boxcar','hanning'});
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');
ylog;
xlog;
axis tight
xtick(10.^(-5:0));
ytick(10.^(-4:2:8));
axis tight;
hv = vlines((M).^-1);

% please vary the length M of the window

%%
% let's reconduct this operation on the time series minus seasonal cycle
r2sf0 = vfilt(r2s,w0,'mirror');
r2sf1 = vfilt(r2s,w1,'mirror');

figure;
hold on
h1 = plot(t,r2s);
h2 = plot(t,[r2sf0 r2sf1],'linewidth',2);
legend([h1;h2],{'r2s','boxcar','hanning'},'location','best');

% let's examine what happens in the spectral domain
[f,s] = mspec(dt,[r2s r2sf0 r2sf1],psi);
figure;
hold on
hs = plot(f/(2*pi),s);
legend(hs,{'r2s','boxcar','hanning'});
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');
ylog;
xlog;
axis tight
xtick(10.^(-5:0));
ytick(10.^(-4:2:8));
axis tight;
hv = vlines((M).^-1);

%%
% what is the high-pass results?

r2h = r2s - r2sf1;

figure;
subplot(2,1,1)
hold on
h1 = plot(t,[r2s r2sf1]);
set(h1(2),'linewidth',2);
legend(h1,{'r2s','r2sf1'},'location','best');
xlim([0 N]);
subplot(2,1,2)
hold on
h2 = plot(t,r2h,'color',cc(3,:));
legend(h2,{'r2h'},'location','best');
xlim([0 N]);

% let's examine what happens in the spectral domain
[f,s] = mspec(dt,[r2s r2sf1 r2h],psi);
figure;
hold on
hs = plot(f/(2*pi),s);
legend(hs,{'r2s','r2sf1','rh2'});
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');
ylog;
xlog;
axis tight
xtick(10.^(-5:0));
ytick(10.^(-4:2:8));
hv = vlines((M).^-1);

%%
% band-pass filtering?
% make a Hanning? window times a complex exponential
% we need an even number of years (not sure why yet)
M = round(4*365);
w = hanning(M);
w = w/sum(w);
% annual filter
p1 = 365.25;
wea = w.*[exp(1i*2*pi*[0:length(w)-1]'./p1) + exp(-1i*2*pi*[0:length(w)-1]'./p1)];
% semi-annual
p2 = 365.25/2;
wesa = w.*[exp(1i*2*pi*[0:length(w)-1]'./p2) + exp(-1i*2*pi*[0:length(w)-1]'./p2)];
% tri-annual
p3 = 365.25/3;
weta = w.*[exp(1i*2*pi*[0:length(w)-1]'./p3) + exp(-1i*2*pi*[0:length(w)-1]'./p3)];

figure;
h = plot([wea wesa weta]);
legend(h,{'Annual','Semi-annual','Tri-annual'});
title('Filters');

% "seasonal" filter
wall = wea+wesa+weta;

x2sf = vfilt(x2,wall,'mirror');

figure;
hh = plot(t,[x2,x2sf,x2-x2sf]);
legend(hh,{'x2','x2sf','residual'});

% what is the "spectrum" of this filter?
dum = fft([wea wesa weta wall],N,1);
swe = 2*abs(dum(1:N/2+1,:)).^2*dt/N;

figure;
h = plot(f/(2*pi),swe,'linewidth',2);
set(h(end),'linestyle','--');
ylin;xlog;
legend(h,{'Annual','Semi-annual','Tri-annual','Sum = Seasonal'},'location','best');
xlim([0 0.5])
vlines([1 2 3]./365.26);

% but this filter is not perfect, use ylog to see this

% let's examine what happens in the spectral domain
[f,s] = mspec(dt,[x2 x2sf x2-x2sf],psi);
figure;
hold on
hs = plot(f/(2*pi),s);
legend(hs,{'x2','x2sf','residual'});
ylabel('PSD (ppm^2 cpd^{-1})');
xlabel('Frequency (1/day = cycle per day)');
ylog;
xlog;
axis tight
xtick(10.^(-5:0));
ytick(10.^(-4:2:8));
hv = vlines([p1 p2 p3].^-1);

% please redo the last section and change the length of the window to see
% what happens

%% Part II
close all;clear all;
% Let's now consider the time series of the North Atlantic meridional
% Overturning Circulation at 26N
% Data are from http://www.bodc.ac.uk/rapidmoc/

t = ncread('moc_transports.nc','time');
tDesc  = ncreadatt('moc_transports.nc','time','long_name');
tUnits = ncreadatt('moc_transports.nc','time','units');

x = ncread('moc_transports.nc','moc_mar_hc10');
xDesc = ncreadatt('moc_transports.nc','moc_mar_hc10','long_name');
xUnits = ncreadatt('moc_transports.nc','moc_mar_hc10','units');

% please conduct your own analyses of this time series using the tools
% explored earlier

