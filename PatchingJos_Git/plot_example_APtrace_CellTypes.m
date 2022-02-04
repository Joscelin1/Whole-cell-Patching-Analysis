 close all
% load WKY example
aAbfData = abfload('P:\Patching\20210715\2021_07_15_C01_0003_SingleAP_01.abf');

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


figure
hold on
plot(WKY(500:10000, 1), 'color', [0.7, 0.7, 0.7], 'LineWidth', 2);
hold off

