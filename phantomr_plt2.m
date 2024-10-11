%#######################################################################
%
%                  * PHANTOM Reduced PLoT 2 Program *
%
%          M-File which reads the T1rho phantom images and plots a
%     selection of the images and segments the center (mid-scan) images.
%     Otsu's method is used to threshold the images with spin lock time
%     = 0 ms.
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
%          The plots are saved to Postscript file:  phantomr?_plt1.ps
%     and PDF file:   phantomr?_plt2.pdf.  A PDF Editor is needed to
%     convert the Postscript file to PDF and combine the two PDF files
%     into one PDF file:  phantomr?_plt2.pdf.  ? is the series number.
%
%          The centroids and radii of the regions of interest and the
%     masks for the vials are saved to Matlab MAT file:
%     phantomr?_plt2.mat.  ? is the series number.
%
%     NOTES:  1.  Matlab MAT file dicom_lst.mat or dicom_lst2.mat must
%             be in the current directory or path.
%
%             2.  Matlab M-file circ_plt.m must be in the current
%             directory or path.
%
%     01-Jul-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Plot and Fitting Parameters
%
na = 144;               % Number of points in circles (2.5 degrees)
nvial = 6;              % Number of vials
istrt = [16; 5];        % Start at slice 17 or 6
nsl = [32; 10];         % Middle 32 or 10 slices
%
% Load Image File Data from dicom_lst.mat or dicom_lst2.mat
%
if exist('dicom_lst.mat','file')==2
  load dicom_lst.mat afiles idvr idx isz nimages pspc sn splcktc stxt;
  afiles = afiles(idx);
  nimages = nimages(idx);
  n = sum(idx);
  ddirs = repmat({'.'},n,1);
elseif exist('dicom_lst2.mat','file')==2
  load dicom_lst2.mat afiles ddirs idvr isz nimages pspc sn splcktc ...
                      stxt;
else
  error([' *** ERROR in phantomr_plt2:  dicom_lst.mat or ', ...
         'dicom_lst2.mat does not exist!']);
end
%
idv = ~startsWith(stxt,'WIP');         % Exclude Philips T1rho values
idvr = idvr&idv;
%
n = sum(idvr);          % Number of series to analyze
ddirr = ddirs(idvr);    % Subdirectories for series
nfiles = nimages(idvr); % Numbers of files
afiler = afiles(idvr);  % T1rho files
isz = isz(idvr,:);      % Image sizes in pixels
sn = sn(idvr);          % Series numbers
pspc = pspc(idvr,:);    % Pixel sizes
spltr = splcktc(idvr);  % T1rho spin lock times
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
% Loop through the Series
%
for j = 1:n
%
% Get Series Specific Indices
%
   ddir = ddirr{j};     % Subdirectory for series
   nfile = nfiles(j);   % Number of image files in this series
   splt = eval(['[' spltr{j} ']'])';   % T1rho spin lock times
   nsplt = size(splt,1);               % Number of spin lock times
   nsls = nfile/nsplt;  % Number of slices
   iszs = isz(j,:);     % Image size
   if nsls==64
     id = 1;            % Index to 64 slices series
   elseif nsls==20
     id = 2;            % Index to 64 slices series
   else
     error([' *** ERROR in phantomr_plt2:  Incorrect number', ...
                 ' of slices in series!']);
   end
   istrts = istrt(id);       % Start of middle slices
   nslm = nsl(id);           % Number of middle slices
   sns = sn(j);         % Series number
   snt = int2str(sns);  % Series number as text
%
% Spin Lock Times as Text
%
   sltxt = cell(1,nsplt);
   sltxt{1} = sprintf('Spin Lock Times =%2i,',splt(1));
   for k = 2:nsplt-1
      sltxt{k} = sprintf(' %2i,',splt(k));
   end
   sltxt{nsplt} = sprintf(' %2i ms',splt(nsplt));
   sltxt = strcat(sltxt{:});
%
% Postscript Plot File Name
%
   pfile = ['phantomr' snt '_plt1.ps'];
%
% Get Image File Names
%
   fnams = afiler{j};   % T1rho image file names
%
% Get Images and Maximum Scaled Image Value
%
   valmx = zeros(nfile,1);
   v = zeros([iszs nfile]);
   for k = 1:nfile
      info = dicominfo(fullfile(ddir,fnams{k}));
      sl = double(info.RescaleSlope);
      y0 = double(info.RescaleIntercept);
      img = dicomread(info);
      img = sl*double(img)+y0;
      v(:,:,k) = img;
      valmx(k) = max(img(:));          % Slice maximum
   end
%
   scmx = 10*fix(max(valmx)/10);       % Round maximum value down
%
% Plot All of the Images with Spin Lock Time = 0 ms
%
   v0 = v(:,:,1:nsplt:nfile);          % Spin lock time = 0 ms
%
   figure;
   orient landscape;
   montage(v0,'DisplayRange',[0 scmx]);
   pos = get(gca,'Position');
   pos(2) = pos(2)/5;
   set(gca,'Position',pos);
   title({['Series ' snt ' Phantom Images']; ...
          ['Slices 1-' int2str(nsls)]; 'Spin Lock Time = 0 ms'}, ...
          'FontSize',16,'FontWeight','bold');
   print('-dpsc2','-r600','-fillpage',pfile);
%
% Plot All of the Spin Lock Time in the Middle of the Scan
%
   imid = nfile/2;
   imid1 = imid-nsplt+1;
   figure;
   orient landscape;
   montage(v(:,:,imid1:imid),'DisplayRange',[0 scmx]);
   pos = get(gca,'Position');
   pos(2) = pos(2)/5;
   set(gca,'Position',pos);
   title({['Series ' snt ' Phantom Images']; 'MidScan Slice'; ...
          sltxt},'FontSize',16,'FontWeight','bold');
   print('-dpsc2','-r600','-fillpage','-append',pfile);
%
% Plot Just the Mid Image with Spin Lock Time = 0 ms
%
   figure;
   orient landscape;
   img = v(:,:,imid1);
   imagesc(img);
   axis image;
   axis off;
   colormap gray;
   title({['Series ' snt ' Phantom Images']; 'MidScan Slice'; ...
          'Spin Lock Time = 0 ms'},'FontSize',16,'FontWeight','bold');
   print('-dpsc2','-r600','-fillpage','-append',pfile);
%
% Middle Slices of Scan
%
   lvl = zeros(nslm,1);   % Threshold levels for each slice
   rm = zeros(nslm,1);    % Threshold metric
   ctrs = zeros(nvial,2,nslm); % Centroids of vials
   rmin = zeros(nvial,nslm);   % Minor axis lengths (minimum radii)
   vnum = int2str((1:nvial)');    % Vial number
   vmsk = false(prod(iszs),nvial,nslm);  % Mask for vials
   [xi,yi] = meshgrid(1:iszs(1),1:iszs(2)); % Pixel coordinates
   xi = xi(:);
   yi = yi(:);
   npix = size(xi,1);      % Number of pixels in a slice
   xyi = [xi yi];
%
   for k = 1:nslm
%
% Threshold Slice
%
      isl = istrts+k;      % Index to slice
      img = v0(:,:,isl);   % Spin lock time = 0 ms
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
          error([' *** ERROR in phantomr_plt2:  Incorrect number', ...
                 ' of regions were identified!']);
        end
      elseif rsize<6
        error([' *** ERROR in phantomr_plt2:  Incorrect number', ...
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
             ['Slice ' int2str(isl)]; 'Spin Lock Time = 0 ms'}, ...
             'FontSize',16,'FontWeight','bold');
      if k==round(nslm/2)
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
   mf = ['phantomr' snt '_plt2.mat'];  % MAT file name
   save(mf,'ctrs','ddir','fnams','idvr','istrts','iszs','lvl', ...
        'nfile','nslm','nsls','nsplt','nvial','rm','rmin','sns', ...
        'snt','splt','sltxt','v','v0','vmsk');
%
% Plot the Six Vial Masks for the Last Slice Image
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
   print('-dpdf','-r600','-fillpage',['phantomr' snt '_plt2.pdf']);
%
   pause;
%
   close all;
%
end
%
return