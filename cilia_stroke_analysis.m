clear

%% Information _ data file

file_name = "sample.tif";   % Name of tif image sequence file
file_video_length = 1;      % Sec
file_video_fps = 1000;      % Hz

%% Information _ ROI line

roi_line = [                % point1-X  point1-Y  point2-X  point2-Y
    270 18  385 209
    381 242 461 492
    464 529 520 814
    529 866 537 1040
    ];

first_point = roi_line(:,1:2);
second_point = roi_line(:,3:4);

if sum(size(first_point) == size(second_point)) ~= 2
    error("The number of first points differs from the number of second points.");
end

%% Parameters _ ROI

distance_between_roi_center = 30;   % Pixels
roi_size = 5;                       % Pixels, ROI = ( Center of ROI - roi_size )  ~  ( Center of ROI + roi_size )

%% Parameters _ drawing

plot_x_limit = 100;                 % Hz
subplot_number = 4;

%% Parameters _ FFT

Fs = file_video_fps;                % Sampling frequency
T = 1/Fs;                           % Sampling period
L = file_video_length * Fs;         % Length of signal
t = (0:L-1)*T;                      % Time vector

%% Parameters _ etc

frequency_calculation_range_min = 10;   % Hz
frequency_calculation_range_max = 40;   % Hz

%% Load image file

data_name = split(file_name,".");
data_name = data_name(1);

data_info = imfinfo(file_name);
image_number_of_frame = length(data_info);
image_bit = data_info(1).BitsPerSample;
image_max_signal = power(2,image_bit)-1;

if image_number_of_frame ~= L
    error("Frame number of video is not matched with FFT parameter");
end

image_temp = imread(file_name, 1);
image = zeros(size(image_temp,2), size(image_temp,1), image_number_of_frame);

progress = 0;
for k = 1 : image_number_of_frame
    if round((k*100)/image_number_of_frame) > progress+1
        progress = progress+1;
        file_name + " loading " + round((k*100)/image_number_of_frame) + " %"
    end
    image(:,:,k) = imread(file_name, k)';
end

%% generate ROI list

[M,I] = sort(first_point(:,1));     % sorting ROI lines based on their first point
first_point = first_point(I,:);
second_point = second_point(I,:);

roi_list = [];                      % initialize roi_list

for i = 1 : size(first_point,1)
    roi_line_x = [first_point(i,1) second_point(i,1)];
    roi_line_y = [first_point(i,2) second_point(i,2)];
    
    if roi_line_x(1) > roi_line_x(2)    % make ROI generation starting point the point that has smaller x coordinate
        roi_line_x = flip(roi_line_x);
        roi_line_y = flip(roi_line_y);
    end
    
    slope = (roi_line_y(2)-roi_line_y(1)) / (roi_line_x(2)-roi_line_x(1));
    roi_x_interval = distance_between_roi_center / sqrt(slope^2 + 1);
    roi_y_interval = roi_x_interval * slope;
    roi_x = roi_line_x(1) : roi_x_interval : roi_line_x(2);
    roi_y = roi_line_y(1) : roi_y_interval : roi_line_y(2);
    
    roi_list = [roi_list; round([roi_x' roi_y'])];
end

%% alnaysis & plot & save

fft_of_all_roi = [];

Fig = figure('Position', [0 0 1600 800]);

subplot(subplot_number,1,1);
hold on;
title("FFT Result (black: raw, red: smoothed)")

progress = 0;
for roi_order = 1:size(roi_list,1)
    if round((roi_order*100)/size(roi_list,1)) > progress+1
        progress = progress+1;
        file_name + " analysis " + round((roi_order*100)/size(roi_list,1)) + " %"
    end
    roi = roi_list(roi_order, :);
    fft_of_all_expansion = [];
    
    for x_cor = (roi(1)-roi_size):(roi(1)+roi_size)
        for y_cor = (roi(2)-roi_size):(roi(2)+roi_size)
            data = squeeze(image(x_cor, y_cor, :));
            
            Y = fft(data);
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L/2))/L;
            f = f(2:end);
            
            fft_of_all_expansion = [fft_of_all_expansion P1(2:end)];
        end
    end
    
    fft_of_all_roi = [fft_of_all_roi mean(fft_of_all_expansion, 2) / max(mean(fft_of_all_expansion, 2))];
    
    plot(f,mean(fft_of_all_roi,2),'color',[0 0 0 5/size(roi_list,1)])
    plot(f,smooth(mean(fft_of_all_roi,2)),'color',[1 0 0 5/size(roi_list,1)])
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    xlim([1 plot_x_limit]);
    
end

subplot(subplot_number,1,2);
hold on;
title("Averaged FFT Result (black: raw, red: smoothed)")
plot(f,mean(fft_of_all_roi,2), 'color',[0 0 0]);
plot(f,smooth(mean(fft_of_all_roi,2)), 'color',[1 0 0]);
xlabel('f (Hz)');
ylabel('|P1(f)|');
xlim([1 plot_x_limit]);

data_heatmap = fft_of_all_roi(1:plot_x_limit,:)';

subplot(subplot_number,1,3);
h = heatmap(data_heatmap);
h.Colormap = flipud(gray);
h.GridVisible = 'off';
cdx = h.XDisplayLabels;
cdy = h.YDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdx,1), size(cdx,2));
h.YDisplayLabels = repmat(' ',size(cdy,1), size(cdy,2));
title("FFT Result (each row = each ROI)")

subplot(subplot_number,1,4);
[M, I] = max(fft_of_all_roi);
[maxSorted, ISorted] = sort(I);
h = heatmap(data_heatmap(ISorted,:));
h.Colormap = flipud(gray);
h.GridVisible = 'off';
cdx = h.XDisplayLabels;
cdy = h.YDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdx,1), size(cdx,2));
h.YDisplayLabels = repmat(' ',size(cdy,1), size(cdy,2));
title("FFT Result (each row = each ROI, reorder ROIs by maximum peak position)")

sgtitle({file_name
    size(roi_list,1) + " ROIs"
    "Mean frequency in"
    "whole range = " + mean(I) + " Hz"
    "between " + frequency_calculation_range_min + " Hz ~ " + frequency_calculation_range_max + " Hz = " + mean(I(I>frequency_calculation_range_min & I<frequency_calculation_range_max)) + " Hz"}, 'interpreter', 'none');

time_label = "_" + string(datetime,'yyyy-MM-dd_HH-mm-ss');

if(~exist("pdf_file", "dir"))
    mkdir("pdf_file");
end
disp("Saving pdf...")
print(Fig, pwd + "\pdf_file\" + data_name + time_label + "_result.pdf", '-dpdf', '-fillpage', '-vector');
disp("pdf saved.")

if(~exist("png_file", "dir"))
    mkdir("png_file");
end
disp("Saving png...")
print(Fig, pwd + "\png_file\" + data_name + time_label + "_result.png", '-dpng');
disp("png saved.")

close all;

Fig = figure;
imshow(image(:,:,1)'/image_max_signal)
hold on
for i = 1 : size(roi_list,1)
    line([roi_list(i,1)-(roi_size+1) roi_list(i,1)-(roi_size+1)], [roi_list(i,2)-(roi_size+1) roi_list(i,2)+(roi_size+1)],'Color','red','LineWidth',1)
    line([roi_list(i,1)-(roi_size+1) roi_list(i,1)+(roi_size+1)], [roi_list(i,2)-(roi_size+1) roi_list(i,2)-(roi_size+1)],'Color','red','LineWidth',1)
    line([roi_list(i,1)+(roi_size+1) roi_list(i,1)-(roi_size+1)], [roi_list(i,2)+(roi_size+1) roi_list(i,2)+(roi_size+1)],'Color','red','LineWidth',1)
    line([roi_list(i,1)+(roi_size+1) roi_list(i,1)+(roi_size+1)], [roi_list(i,2)+(roi_size+1) roi_list(i,2)-(roi_size+1)],'Color','red','LineWidth',1)
end

if(~exist("roi_file", "dir"))
    mkdir("roi_file");
end
disp("Saving roi...")
print(Fig, pwd + "\roi_file\" + data_name + time_label + "_roi.png", '-dpng');
disp("roi saved.")

close all;

if(~exist("mat_file", "dir"))
    mkdir("mat_file");
end
disp("Saving mat...")
save(pwd + "\mat_file\" + data_name + time_label + ".mat");
disp("mat saved.")

disp("Analysis of " + data_name + time_label + ".tif is finished.")
disp("whole range = " + mean(I) + " Hz");
disp("between " + frequency_calculation_range_min + " Hz ~ " + frequency_calculation_range_max + " Hz = " + mean(I(I>frequency_calculation_range_min & I<frequency_calculation_range_max)) + " Hz");