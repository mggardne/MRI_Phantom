%#######################################################################
%
%                    * PHANTOM Star PLoT Program *
%
%          M-File which reads the T2* phantom images and plots a
%     selection of the images and segments the center (mid-scan) images.
%     Otsu's method is used to threshold the images with echo time
%     = 1 ms.
%
%          At the top of the vial, OD = 29.72 mm and ID = 27.81.  The
%     vial wall thickness is half the difference in diameters:
%     29.72-27.81 = 1.91 mm (difference in vial diameters) and
%     1.91 mm/2 = 0.955 mm (vial wall thickness).  The image pixel
%     size is 0.273 mm.  Converting the wall thickness to pixels =
%     0.955 mm/0.273 mm/pixel = 3.50 pixels.  Adding 10 pixels to
%     avoid partial volume effects = 3.5+10 = 13.5 pixels.  Reducing
%     the radii by 13.5 pixels will avoid the signal from the
%     polypropylene vial material.
%
%          The plots are saved to Postscript file:  phantoms?_plt1.ps
%     and PDF file:   phantoms?_plt.pdf.  A PDF Editor is needed to
%     convert the Postscript file to PDF and combine the two PDF files
%     into one PDF file:  phantoms?_plt.pdf.  ? is the initial series
%     number.
%
%          The centroids and radii of the regions of interest and the
%     masks for the vials are saved to Matlab MAT file:
%     phantoms?_plt.mat.  ? is the initial series number.
%
%     NOTES:  1.  Matlab MAT file dicom_lst2.mat must be in the current
%             directory or path.
%
%             2.  Matlab M-file circ_plt.m must be in the current
%             directory or path.
%
%     27-Jul-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Plot and Fitting Parameters
%
na = 144;               % Number of points in circles (2.5 degrees)
nvial = 6;              % Number of vials
% istrt = [16; 5];        % Start at slice 17 or 6
% nsl = [32; 10];         % Middle 32 or 10 slices
istrt = 12;             % Start at slice 13
nsl = 24;               % Middle 24 slices
%
% Load Image File Data from dicom_lst2.mat
%
load dicom_lst2.mat afiles ddirs etn isz nimages pspc sn stxt;
idv = contains(stxt,'T2*');            % Get T2* series
n = sum(idv);           % Number of series/echo times to analyze
%
etns = [etn{idv}];      % T2* echo times
[etns,ids] = sort(etns);               % Sorted echo times
idv = find(idv);
idv = idv(ids);         % Sorted index to T2* series
%
ddirr = ddirs(idv);     % Subdirectories for series
nfiles = nimages(idv);  % Numbers of files
d = diff(nfiles);
if any(d)
  error([' *** ERROR in phantoms_plt:  Series do not have the', ...
         ' same number of slices!']);
end
nfile = nfiles(1);      % Numbers of files/number of slices
%
afiless = afiles(idv);  % T2* files
isz = isz(idv,:);       % Image sizes in pixels
%
d = diff(isz);
if any(any(d))
  error([' *** ERROR in phantoms_plt:  Images do not have the', ...
         ' same size!']);
end
iszs = isz(1,:);
%
sn = sn(idv);           % Series numbers
sns = sn(1);            % Series number (first of the series)
snt = int2str(sns);     % Series number as text
pspc = pspc(idv,:);     % Pixel sizes
%
% Echo Times as Text
%
etxt = cell(1,n);
etxt{1} = sprintf('Echo Times = %4.2f,',etns(1));
etxt{2} = sprintf(' %4.2f,',etns(2));
for k = 3:n-1
   etxt{k} = sprintf(' %2i,',round(etns(k)));
end
etxt{n} = sprintf(' %2i ms',etns(n));
etxt = strcat(etxt{:});
etxt = strrep(etxt,'  ',' ');
%
% Vial In-Plane Dimensions
%
OD = 29.72;             % Vial outside diameter in mm
ID = 27.81;             % Vial inside diameter in mm
dr = (OD-ID)/2;         % Vial wall thickness in mm
npv = dr./min(pspc,[],2);    % Vial wall thickness in pixels
npv = ceil(100*npv)/100;     % Round up in 1/100ths of pixel
% npv = npv+2;            % Adding 2 pixels to avoid partial volume effects
npv = npv+10;           % Adding 10 pixels to avoid partial volume effects
%
IDp = ID./mean(pspc,2); % Vial inside diameter in pixels
llim = 0.85*IDp;        % Lower limit on vial diameters
ulim = 1.1*IDp;         % Upper limit on vial diameters
%
v = zeros([iszs nfile n]);
valmx = zeros(nfile,n);
%
% Loop through the Series/Echo Times
%
for j = 1:n
%
% Get Series Specific Directory and Image File Names
%
   ddir = ddirr{j};     % Subdirectory for series
   fnams = afiless{j};  % T2* image file names
%
% Get Images and Maximum Scaled Image Value
%
   for k = 1:nfile
      fnam = fullfile(ddir,fnams{k});
      info = dicominfo(fnam);
      sl = double(info.RescaleSlope);
      y0 = double(info.RescaleIntercept);
      img = dicomread(fnam);
      img = sl*double(img)+y0;
%
      v(:,:,k,j) = img;
      valmx(k,j) = max(img(:));
   end
%
end
%
% Postscript Plot File Name
%
pfile = ['phantoms' snt '_plt1.ps'];
%
% Plot Images with Echo Time = 1 ms
%
k1 = 3;                 % Third series - echo time = 1 ms (usually = 2)
scmx = 100*fix(max(valmx(:,k1))/100);       % Round maximum value down
%
figure;
orient tall;
montage(squeeze(v(:,:,:,k1)),'DisplayRange',[0 scmx],'Size',[8 6]);
pos = get(gca,'Position');
pos(2) = pos(2)/5;
set(gca,'Position',pos);
title({['Series ' snt ' Phantom Images']; ...
       ['Slices 1-' int2str(nfile)]; 'Echo Time = 1 ms'}, ...
       'FontSize',16,'FontWeight','bold');
print('-dpsc2','-r600','-fillpage',pfile);
%
% Plot All of the Echo Times in the Middle of the Scan
%
imid = nfile/2;         % Middle slice
scmx = 100*fix(max(valmx(imid,:))/100);     % Round maximum value down
%
figure;
orient landscape;
montage(squeeze(v(:,:,imid,:)),'DisplayRange',[0 scmx]);
pos = get(gca,'Position');
pos(2) = pos(2)/5;
set(gca,'Position',pos);
title({['Series ' snt ' Phantom Images']; 'MidScan Slice'; etxt}, ...
      'FontSize',16,'FontWeight','bold');
print('-dpsc2','-r600','-fillpage','-append',pfile);
%
% Plot Just the Mid Image with Echo Time = 1 ms
%
figure;
orient landscape;
img = squeeze(v(:,:,imid,k1));
imagesc(img);
axis image;
axis off;
colormap gray;
title({['Series ' snt ' Phantom Images']; ...
        'MidScan Slice'; 'Echo Time = 1 ms'},'FontSize',16, ...
        'FontWeight','bold');
print('-dpsc2','-r600','-fillpage','-append',pfile);
%
% Initialize Variables for the Middle Slices of Scan
%
lvl = zeros(nsl,1);     % Threshold levels for each slice
rm = zeros(nsl,1);      % Threshold metric
ctrs = zeros(nvial,2,nsl);   % Centroids of vials
rmin = zeros(nvial,nsl);     % Minor axis lengths (minimum radii)
vnum = int2str((1:nvial)');  % Vial number
vmsk = false(prod(iszs),nvial,nsl);    % Mask for vials
[xi,yi] = meshgrid(1:iszs(1),1:iszs(2));    % Pixel coordinates
xi = xi(:);
yi = yi(:);
npix = size(xi,1);      % Number of pixels in a slice
xyi = [xi yi];
%
for k = 1:nsl
%
% Threshold Slice
%
   isl = istrt+k;       % Index to slice
   img = squeeze(v(:,:,isl,k1));    % Echo time = 1 ms
   [lvl(k),rm(k)] = multithresh(img);
   bw = img>lvl(k);
%
% Get Properties of the Threshold Regions (Vials)
%
   rdat = regionprops('table',bw,'Centroid', ... % 'MajorAxisLength', ...
                      'MinorAxisLength');
   rsize = size(rdat,1);
   if rsize>6
     idv = rdat.MinorAxisLength>llim(j)&rdat.MinorAxisLength<ulim(j);
     rdat = rdat(idv,:);
     rsize = size(rdat,1);
     if rsize~=6
       error([' *** ERROR in phantoms_plt:  Incorrect number', ...
              ' of regions were identified!']);
     end
   elseif rsize<6
     error([' *** ERROR in phantoms_plt:  Incorrect number', ...
            ' of regions were identified!']);
   end
   ctrs(:,:,k) = rdat.Centroid;
   rmin(:,k) = rdat.MinorAxisLength/2-npv(j); % Reduce radius to avoid vial material
%
% Plot Threshold Regions
%
   [xp,yp] = circ_plt(rmin(:,k),ctrs(:,:,k),na);
   figure;
   orient landscape;
   imagesc(img);
   axis image;
   axis off;
   colormap gray;
   hold on;
   plot(xp,yp,'r-','LineWidth',1);
   plot(ctrs(:,1,k),ctrs(:,2,k),'r+','LineWidth',1);
   text(ctrs(:,1,k)+3,ctrs(:,2,k)+4,vnum,'Color','r','FontSize',12);
   title({['Series ' snt ' Phantom Images']; ...
          ['Slice ' int2str(isl)]; 'Echo Time = 1 ms'}, ...
          'FontSize',16,'FontWeight','bold');
   if k==round(nsl/2)
     print('-dpsc2','-r600','-fillpage','-append',pfile);
   end
%
% Create a Mask
%
   for l = 1:nvial
      rc = xyi-repmat(ctrs(l,:,k),npix,1);
      rc = sum(rc.*rc,2);              % Radius squared
      vmsk(:,l,k) = rc<rmin(l,k).*rmin(l,k);
   end
end
%
% Save Images and Masks
%
mf = ['phantoms' snt '_plt.mat'];   % MAT file name
save(mf,'ctrs','ddir','etns','etxt','fnams','idv','istrt','iszs', ...
     'lvl','n','nfile','nsl','nvial','rm','rmin','sns','snt','v', ...
     'vmsk');
%
% Plot the Six Vial Masks for the Last Middle Slice Image
%
figure;
orient landscape;
colormap gray;
%
rmx = double(max(img(:)));
%
for l = 1:6
   subplot(2,3,l);
   imagesc(img);
   axis image;
   axis off;
   hold on;
   h = imagesc(reshape(vmsk(:,l,k)*rmx,iszs));
   set(h,'AlphaData',0.5);
   title(['Vial ' int2str(l) ' Mask'],'FontSize',16, ...
         'FontWeight','bold');
end
%
sgtitle({['Series ' snt ' Phantom Images']; ['Slice ', ...
         int2str(isl)]},'FontSize',20,'FontWeight','bold');
%
print('-dpdf','-r600','-fillpage',['phantoms' snt '_plt.pdf']);
%
return