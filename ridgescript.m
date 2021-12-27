% script reads dump files with xyz position data
% saves b&w png images with ridges, valleys, and all ripples marked as points
%
% Author: Landon Cordova
% Rutgers University 
% 2021


set(0,'DefaultFigureVisible','off')

for qq=1:2
close all
clearvars -except qq
imgnames=["init_ridges.png", "init_peaks.png", "init_valleys.png"; "post_ridges.png", "post_peaks.png", "post_valleys.png"];
files = dir(fullfile('.'));
outputnames = ["outputinit", "outputpost"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input Paramaters
xbin=1000;                         %Number of linear interpolation points
ybin=900;                          % approx. multiply Nx|Ny by ??? for 1 atom per bin
sigma = 25;                        %Sigma value used in Guassian smoothing        
slope_thresh = 3e-3;              %First derivative test threshold
ratio_thresh = 1;                  %Ripple aspect ratio threshold
curvecut=0.01;
xborder=10;                        
yborder=125;					   %creates borders that are not analyzed to avoid boundary region

isCurved = 0;

%when looking at curved data gets cylinder fit from data
h= dlmread('fcurvature',' ',[1 0 1 0]);    %y cordinate of cylinder fit axis
k= dlmread('fcurvature',' ',[2 0 2 0]);    %z cordinate of cylinder fit axis
r = dlmread('fcurvature',' ',[0 0 0 0]);     %cylinder fit radius
if r>10000 
    r=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseName1 = 'dump.init.';
baseName2 = 'dump.postindent.';

%Grab atomic positions as x,y,z coordinates
if qq==1
    baseName = baseName1;
end
if qq==2
    baseName = baseName2;
    isCurved = 1;
end
finalstep=0;
largeststep=0;

%find final last dump file matching pattern
for w = 1:length(files)
    fileName = files(w).name;
    if  contains(fileName,baseName)
        step = str2num(erase(fileName,baseName));
        if step>largeststep
           finalstep=w;
           largeststep=step;
        end
    end
end
    
    %file = strcat(mydir,fileName);
    file = strcat(files(finalstep).name);
    
    %grab the step number
    step = finalstep;
    
    %grab the box info
    tmp = dlmread(file,' ',[5 0 7 1]);
    xlo = tmp(1,1);
    xhi = tmp(1,2);
    
    ylo = tmp(2,1);
    yhi = tmp(2,2);
    
    zlo = tmp(3,1);
    zhi = tmp(3,2);
    
    %grab the atomic positions
    tmp = dlmread(file,' ',9,0);
    x = (xhi-xlo)*tmp(:,3)+xlo;
    y = (yhi-ylo)*tmp(:,4)+ylo; 
    z = (zhi-zlo)*tmp(:,5)+zlo;
    if isCurved
        z =-(((y-h).^2+(z-k).^2).^(1/2)-r);
    end

%Create surface 
xv = linspace(min(x), max(x), xbin);
yv = linspace(min(y), max(y), ybin);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
Z(isnan(Z))=0;
Z = imgaussfilt(Z,sigma);
hx = (max(x)-min(x))/xbin;
hy = (max(y)-min(y))/ybin;

%Calculate 1st and 2nd derivatives for Hessian matrix
[Zx, Zy]=gradient(Z,hx,hy);
[Zxx, Zxy]=gradient(Zx,hx,hy);
[Zyx, Zyy]=gradient(Zy,hx,hy);

%Allocate memoory
[len1, len2]=size(Zxx);
e1=zeros(len1,len2); %major eigenvalue
e2=zeros(len1,len2); %minor eigenvalue
ea=zeros(len1,len2); %x component of eigenvector
eb=zeros(len1,len2); %y component of eigenvector
rp=false(len1,len2); %matrix of ridge points
vp=false(len1,len2); %matrix of valley points

%Calculate Hessian, eigenvalues and eigenvectors of Hessian
for i=1:len1
    for j=1:len2
        A=[Zxx(i,j) Zxy(i,j); Zyx(i,j) Zyy(i,j)];
        [V, D] = eig(A);
        if abs(D(1,1))>abs(D(2,2))
            e1(i,j)=D(1,1);
            e2(i,j)=D(2,2);
            C=V(:,1);
        else
            e1(i,j)=D(2,2);
            e2(i,j)=D(1,1);
            C=V(:,2);
        end
        ea(i,j)=C(1);
        eb(i,j)=C(2);
    end
end

%Threshold each point to determine if it is on a ridge or valley
for i=yborder:len1-yborder
    for j=xborder:len2-xborder
        lam1=abs(e1(i,j));
        lam2=abs(e2(i,j));
        if abs(Zx(i,j)*ea(i,j)+Zy(i,j)*eb(i,j))<slope_thresh && lam1>ratio_thresh*lam2
            if e1(i,j)<(-1*curvecut) || e2(i,j)<(-1*curvecut)
                rp(i,j)=1;
            elseif e1(i,j)>curvecut || e2(i,j)>curvecut
                vp(i,j)=1;
            end
        end
    end
end

%Plotting and saving
figure;
hold on
scatter(X(rp),Y(rp),'.','k');
scatter(X(vp),Y(vp),'.','k');
axis([min(x) max(x) min(y) max(y)])
axis off
hold off
daspect([1 1 2])
saveas(gcf,imgnames(qq,1))

figure;
hold on
scatter(X(rp),Y(rp),'.','k');
axis([min(x) max(x) min(y) max(y)])
axis off
hold off
daspect([1 1 2])
saveas(gcf,imgnames(qq,2))

figure;
hold on
scatter(X(vp),Y(vp),'.','k');
axis([min(x) max(x) min(y) max(y)])
axis off
hold off
daspect([1 1 2])
saveas(gcf,imgnames(qq,3))

%perform skeletonization and create connected graph of points for each ripple 
linedetect;
%sum component lengths of each skeletonized ripple
nodeanalyze;
%outputs information about ripple lengths
dlmwrite(outputnames(qq), myvalues,'delimiter',' ')
end




