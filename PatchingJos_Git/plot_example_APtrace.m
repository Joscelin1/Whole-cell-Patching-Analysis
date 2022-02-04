 close all
% load WKY example
aAbfData = abfload('P:\Patching\20210715\2021_07_15_C02_0003_SingleAP_01.abf');

kk=1;
% this bit inseted by JS 22.07.21 to compensated for bridge ballance
% this is how much before and after the vertical step that is being
% removed
gapbefore = 2;
gapafter = 15;
%loop through the sweeps in the file


%the derivative of channel2 aka aDiffStim to find the edges of the square wave
aDiffStim = diff(aAbfData(:,2,kk));
[maxval, maxind] = max(aDiffStim);
[minval, minind] = min(aDiffStim);


channel1 = aAbfData(:,1,kk);
%to find difference in height of signal increase
offset = channel1(maxind+gapafter)- channel1(maxind-gapbefore);
%drop the middel chunk
channel1(maxind:minind) = channel1(maxind:minind) - offset;
%remove sections arround shift and replace with linear gradient
%section
channel1(maxind-gapbefore:maxind+gapafter-1) = linspace(channel1(maxind-gapbefore), channel1(maxind+gapafter), gapbefore+gapafter);
channel1(minind-gapbefore:minind+gapafter-1) = linspace(channel1(minind-gapbefore), channel1(minind+gapafter), gapbefore+gapafter);
%convert back to aDbf file so works with Jesses code
WKY = smooth(channel1);

WKY_abf = aAbfData(:,1,kk);


%% load SHR example
aAbfData = abfload('P:\Patching\20210617\2021_06_17_C02_0002_SingleAP_01.abf');

% this bit inseted by JS 22.07.21 to compensated for bridge ballance
% this is how much before and after the vertical step that is being
% removed
gapbefore = 2;
gapafter = 15;

%the derivative of channel2 aka aDiffStim to find the edges of the square wave
aDiffStim = diff(aAbfData(:,2,kk));
[maxval, maxind] = max(aDiffStim);
[minval, minind] = min(aDiffStim);

channel1 = aAbfData(:,1,kk);
%to find difference in height of signal increase
offset = channel1(maxind+gapafter)- channel1(maxind-gapbefore);
%drop the middel chunk
channel1(maxind:minind) = channel1(maxind:minind) - offset;
%remove sections arround shift and replace with linear gradient
%section
channel1(maxind-gapbefore:maxind+gapafter-1) = linspace(channel1(maxind-gapbefore), channel1(maxind+gapafter), gapbefore+gapafter);
channel1(minind-gapbefore:minind+gapafter-1) = linspace(channel1(minind-gapbefore), channel1(minind+gapafter), gapbefore+gapafter);
%convert back to aDbf file so works with Jesses code
SHR = channel1;


SHR_abf = aAbfData(:,1,kk);
%make baseline at begining same for both traces 
baseline_WKY = mean(WKY(1:1000));
baseline_SHR = mean(SHR(1:1000));
shift = diff([baseline_WKY, baseline_SHR]);
SHR_shift = smooth(SHR(:,1)- shift) ;

figure
hold on
plot(WKY(500:10000, 1), 'b', 'LineWidth', 2);
plot(SHR_shift(500:10000, 1), 'r', 'LineWidth', 2);
hold off

% figure
% hold on 
% plot(WKY_abf, 'b');
% plot(SHR_abf, 'r');
% hold off
