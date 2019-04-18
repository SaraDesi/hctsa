clc, clear all, close all

addpath('quaternion_library');

%% Orientamento testa

cartella = uigetdir();
full_path = fullfile(cartella,'TP','Dati_Grezzi.txt');
d = load(full_path);

FS = 100;
SamplePeriod = 1/256;
acc = d(:,[5,6,7]);
acc = acc./2^15*8; %g
gyr = d(:,[8,9,10]);
gyr = gyr./2^15*1000*pi/180; %rad/s
mag = d(:,[11,12,13]);
mag = mag./2^15*12; %Gauss
mag_raw = mag;
t = d(:,1)./FS;

load('calib_mag');
xc=((mag(:,1)-Os(1)))/Assis(1)/R;
yc=((mag(:,2)-Os(2)))/Assis(2)/R;
zc=((mag(:,3)-Os(3)))/Assis(3)/R;
mag = [xc yc zc];

% Sensor fusion
%AHRS = MadgwickAHRS('SamplePeriod', 1/100, 'Beta', 0.005);
AHRS = MahonyAHRS('SamplePeriod', SamplePeriod, 'Kp', 0.007);

quaternion = zeros(length(t), 4);
for j = 1:length(t)
    AHRS.Update(gyr(j,:), acc(j,:), mag(j,:));	% gyroscope units must be radians
    quaternion(j, :) = AHRS.Quaternion;
end

euler = quatern2euler(quaternConj(quaternion)) * (180/pi)/SamplePeriod/FS;	% use conjugate, 'quaternConj', for sensor frame relative to Earth and convert to degrees.
eulertemp(:,1) = euler(:,2);
eulertemp(:,2) = -euler(:,3);
eulertemp(:,3) = -euler(:,1);

plot(t,eulertemp(:,3))
pts_inizio_fine = ginput(2);
idx_if = knnsearch(t,pts_inizio_fine(:,1));
euler1 = eulertemp(idx_if(1):idx_if(2),3);
t_if = t(idx_if(1):idx_if(2));
figure
plot(t_if,euler1)

%% Traiettoria

full_path2 = fullfile(cartella,'RF','Position.txt');
R1 = readtable(full_path2);
R1.PosY = -R1.PosY;
figure
plot(R1.PosX,R1.PosY)
xlabel('metri')
ylabel('metri')
title('Traiettoria')
axis equal
hold on
pts = ginput(6); % tre ostacoli da aggirare
idx1 = knnsearch([R1.PosX R1.PosY],pts([1,2],:));
x1 = R1.PosX(idx1(1):idx1(2)) ;
y1 = R1.PosY(idx1(1):idx1(2)) ;
idx2 = knnsearch([R1.PosX R1.PosY],pts([3,4],:));
x2 = R1.PosX(idx2(1):idx2(2)) ;
y2 = R1.PosY(idx2(1):idx2(2)) ;
idx3 = knnsearch([R1.PosX R1.PosY],pts([5,6],:));
x3 = R1.PosX(idx3(1):idx3(2)) ;
y3 = R1.PosY(idx3(1):idx3(2)) ;
plot(x1,y1,'r')
plot(x2,y2,'r')
plot(x3,y3,'r')

m = (y1(end)-y1(1))/(x1(end)-x1(1)); % coeff angolare
n = y1(1)-m*x1(1);
pm = zeros(length(x1),2); % punti medi tra i punti e i loro simmetrici
pm(:,1) = (x1 + m*y1 - m*n)./(m^2 + 1);
pm(:,2) = m*pm(:,1) + n;
S = 2*pm - [x1 y1]; % punti simmetrici rispetto alla linea
ellipse_pts1 = [x1 y1; S];

m = (y2(end)-y2(1))/(x2(end)-x2(1));
n = y2(1)-m*x2(1);
pm = zeros(length(x2),2);
pm(:,1) = (x2 + m*y2 - m*n)./(m^2 + 1);
pm(:,2) = m*pm(:,1) + n;
S = 2*pm - [x2 y2];
ellipse_pts2 = [x2 y2; S];

m = (y3(end)-y3(1))/(x3(end)-x3(1));
n = y3(1)-m*x3(1);
pm = zeros(length(x3),2);
pm(:,1) = (x3 + m*y3 - m*n)./(m^2 + 1);
pm(:,2) = m*pm(:,1) + n;
S = 2*pm - [x3 y3];
ellipse_pts3 = [x3 y3; S];

ax = gca;
ellipse1 = fit_ellipse(ellipse_pts1(:,1),ellipse_pts1(:,2),ax);
ellipse2 = fit_ellipse(ellipse_pts2(:,1),ellipse_pts2(:,2),ax);
ellipse3 = fit_ellipse(ellipse_pts3(:,1),ellipse_pts3(:,2),ax);

%% Parametri spazio-temporali del passo

full_path = fullfile(cartella,'RF','Analisi_Passo.txt');
ap = readtable(full_path);
vel = mean(ap.Velocita__km_h_)*1000/60; %aggiungi la terza
lun = mean(ap.LunghezzaPasso_m_);

%% Parametri

[Tempo,Potenza,Var,DSJ,Max] = parametri_yaw(euler1,t_if,FS);
ra = mean([ellipse1.a,ellipse2.a]);
rb = mean([ellipse1.b,ellipse2.b]);