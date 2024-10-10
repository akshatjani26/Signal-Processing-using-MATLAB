%% Extract the frames from video 
clc
clear
v = VideoReader('ACM_2716_cylinder_water_tunneltrim_3.MP4');
n = 1;
while hasFrame(v)
    img = readFrame(v); 
    imwrite(rgb2gray(img),strcat('img_', num2str(n), '.png'));
    n=n+1; 
end 
%% Extracting the value of pixel for P1, P2, P3 for each frame
clear all 
clc
tstps = (0:(35/1059):35)'; 
tstps(end)=[]; 
%Enter the x & y coordinate of the pixel.
x1 = 2260;
x2 = 1930; 
x3 = 1580; 
y = 1397;
%Creation of an empty array
m1 = [ ]; 
m2 = [ ]; 
m3 = [ ]; 
i = 1; 
%While loop reads each and every frame (i.e. an image) of the video saved in the directory (images are saved in the form of img_1.png, img_2.png ...... upto img_1059.png)
while i<=1059
    I = imread(strcat('img_', num2str(i), '.png'));
    %impixel() reads the value of the pixel at the specified location in the image. 
    c1 = impixel(I, x1, y);
    c2 = impixel(I, x2, y);
    c3 = impixel(I, x3, y);
    %It adds up the pixel values as elements in the column in empty array 'm'
    m1(end+1, 1) = c1(1,1);
    m2(end+1, 1) = c2(1,1);
    m3(end+1, 1) = c3(1,1);
    i = i+1; 
end  
P1f = m1 - mean(m1); 
P2f = m2 - mean(m2); 
P3f = m3 - mean(m3); 
TS_f = [tstps P1f P2f P3f];
%% Normalized Correlation Coefficient
clc
tiledlayout(3,1)
%P1 & P1
[CC1, lags1] = xcorr(TS_f(:, 2),  'normalized');
%P1 and P2
[CC2, lags2] = xcorr(TS_f(:, 2), TS_f(:, 3),  'normalized');
%P1 and P3
[CC3, lags3] = xcorr(TS_f(:, 2), TS_f(:, 4), 'normalized');
xlim1 = find(lags1==0); 
xlim2 = find(lags2==0); 
xlim3 = find(lags3==0); 
CC1 = CC1(xlim1:end); 
CC2 = CC2(xlim2:end); 
CC3 = CC3(xlim3:end); 
lags1 = lags1(xlim1:end); 
lags2 = lags2(xlim2:end); 
lags3 = lags3(xlim3:end); 
table1 = [lags1' CC1 CC2 CC3]; 
nexttile
plot(lags1, CC1)
xlabel('lags')
ylabel('Correlation Coefficient')
title('Normalized Correlation on P1')
axis tight
nexttile
plot(lags2, CC2)
xlabel('lags')
ylabel('Correlation Coefficient')
title('Normalized Correlation on P1 and P2')
axis tight
nexttile
plot(lags3, CC3)
xlabel('lags')
ylabel('Correlation Coefficient')
title('Normalized Correlation on P1 and P3')
axis tight
xlabel('lags')
ylabel('Correlation Coefficient')
%% Finding the convection velocity
clc
delx = [0 330];
frame1 = find(CC1 == max(CC1));
frame2 = find(CC2 == max(CC2));
frame3 = find(CC3 == max(CC3));
t1 = lags1(frame1);
t2 = lags2(frame2);
t3 = lags3(frame3);
delt = [t1 t2]/30;
plot(delt, delx); 
title('del(x) vs del(tau)')
ylabel('delx')
xlabel('deltau')
%V = (delx(3) - delx(2))*30/(t3-t2); %Convection Velocity

%The convection velocity obtained is 228 pixels/second" 
