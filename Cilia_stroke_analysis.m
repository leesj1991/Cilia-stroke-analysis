clear
tic

%% Edit here
imageFullFileName = "third_ventricle_2.tif";
firstPoint = [ 389 210 ];
secondPoint = [ 525 833 ];

expansion = 5;
figureLimit = 100; % Hz
meanCalculationMin = 10; % Hz
meanCalculationMax = 40; % Hz

%% file read
dataName = split(imageFullFileName,".");
dataName = dataName(1);
info = imfinfo(imageFullFileName);
numberOfPages = length(info);
if numberOfPages <1000
    disp("Frame 부족");
    return
end
tempImage = imread(imageFullFileName, 1);
fullImage = zeros(size(tempImage,2), size(tempImage,1), numberOfPages);


for k = 1 : numberOfPages
    imageFullFileName + " loading " + round((k*100)/numberOfPages) + " %"
    fullImage(:,:,k) = imread(imageFullFileName, k)';
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all

%% ROI line

% 2
% firstPoint = [ 389 210 ];
% secondPoint = [ 525 833 ];

% 5
% firstPoint = [ 220 305 ];
% secondPoint = [ 451 821 ];

%% Parameters

subPlotNumber = 4;

%% make ROI list

x = [firstPoint(1) secondPoint(1)];
y = [firstPoint(2) secondPoint(2)];

if abs(x(1)-x(2)) < abs(y(1)-y(2))
    roiX = x(1):1*((x(2)-x(1))/abs((x(2)-x(1)))):x(2);
    roiY = round(y(1):(abs(y(1)-y(2))*((y(2)-y(1))/abs((y(2)-y(1)))))/(length(roiX)-1):y(2));
else
    roiY = y(1):1*((y(2)-y(1))/abs((y(2)-y(1)))):y(2);
    roiX = round(x(1):(abs(x(1)-x(2))*((x(2)-x(1))/abs((x(2)-x(1)))))/(length(roiY)-1):x(2));
end

roiSet = [roiX' roiY'];

%% alnaysis

fftOfAllRoi = [];

Fig = figure;
subplot(subPlotNumber,1,1);
hold on

for roiOrder = 1:size(roiSet,1)
    imageFullFileName + " analysis " + round((roiOrder*100)/size(roiSet,1)) + " %"
    roi = roiSet(roiOrder, :);
    fftOfAllExpansion = [];
    
    for x_cor = (roi(1)-expansion):(roi(1)+expansion)
        for y_cor = (roi(2)-expansion):(roi(2)+expansion)
            data = squeeze(fullImage(x_cor, y_cor, :));
            
            %plot(data)
            
            Fs = 1000;            % Sampling frequency
            T = 1/Fs;             % Sampling period
            L = 1000;             % Length of signal
            t = (0:L-1)*T;        % Time vector
            
            Y = fft(data);
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L/2))/L;
            f = f(2:end);
            
            fftOfAllExpansion = [fftOfAllExpansion P1(2:end)/sum(P1(2:end))];
        end
    end
    fftOfAllRoi = [fftOfAllRoi mean(fftOfAllExpansion, 2) / max(mean(fftOfAllExpansion, 2))];
    
    plot(f,smooth(mean(fftOfAllRoi,2)),'color',[0 0 0 5/size(roiSet,1)])
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    xlim([1 figureLimit]);
    
end

subplot(subPlotNumber,1,2);
plot(f,smooth(mean(fftOfAllRoi,2)), 'color',[0 0 0]);
xlabel('f (Hz)');
ylabel('|P1(f)|');
xlim([1 figureLimit]);

hmdata = fftOfAllRoi(1:figureLimit,:)';

subplot(subPlotNumber,1,3);
h = heatmap(hmdata);
h.Colormap = flipud(gray);
h.GridVisible = 'off';
cdx = h.XDisplayLabels;
cdy = h.YDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdx,1), size(cdx,2));
h.YDisplayLabels = repmat(' ',size(cdy,1), size(cdy,2));

subplot(subPlotNumber,1,4);
[M, I] = max(fftOfAllRoi);
[maxSorted, ISorted] = sort(I);
h = heatmap(hmdata(ISorted,:));
h.Colormap = flipud(gray);
h.GridVisible = 'off';
cdx = h.XDisplayLabels;
cdy = h.YDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdx,1), size(cdx,2));
h.YDisplayLabels = repmat(' ',size(cdy,1), size(cdy,2));

sgtitle({imageFullFileName 
    length(roiX) + " ROIs ( [" + firstPoint(1) + ", " + firstPoint(2) + "] ~ [" + secondPoint(1) + ", " + secondPoint(2) + "] )"
    "Mean frequency in" 
    "whole range = " + mean(I) + " Hz" 
    "between " + meanCalculationMin + " Hz ~ " + meanCalculationMax + " Hz = " + mean(I(I>meanCalculationMin & I<meanCalculationMax)) + " Hz"}, 'interpreter', 'none');

timeLabel = "_" + convertCharsToStrings(datestr(now,'yyyy-mm-dd_HH-MM-SS'));

print(Fig, dataName + timeLabel + ".pdf", '-dpdf', '-fillpage', '-painters');
close all;

imshow(fullImage(:,:,1)'/65520)
hold on
line(x, y,'Color','red','LineWidth',1)
print(dataName + timeLabel + "_roi.png", '-dpng');
close all;

"Analysis of " + dataName + timeLabel + ".tif is finished."

save(dataName + timeLabel + ".m");

toc

mean(I)
mean(I(I>meanCalculationMin & I<meanCalculationMax))
