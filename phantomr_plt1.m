%#######################################################################
%
%                   * PHANTOM Star PLoT 1 Program *
%
%          M-File which reads the T1rho phantom images and plots a
%     selection of the images and segments the center (mid-scan) images.
%     Otsu's method is used to threshold the images with spin lock time
%     = 0 ms.
%
%          The radius of the one big vial is reduced by half to avoid
%     the signal from the edge of the vial.
%
%          The plots are saved to a Postscript file:  phantomr?_plt1.ps
%     where ? is the series number.
%
%          The centroids and radii of the region of interests and the
%     masks for the vial are saved to Matlab MAT file:
%     phantomr?_plt1.mat.  ? is the series number.
%
%     NOTES:  1.  Matlab MAT file dicom_lst2.mat must be in the current
%             directory or path.
%
%             2.  Matlab M-file circ_plt.m must be in the current
%             directory or path.
%
%     30-Nov-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Plot and Fitting Parameters
%
na = 144;               % Number of points in circles (2.5 degrees)
nvial = 1;              % Number of vials
% istrt = [16; 5];        % Start at slice 17 or 6
% nsl = [32; 10];         % Middle 32 or 10 slices
istrt = 24;             % Start at slice 24
nsl = 16;               % Middle 16 slices
%
% Load Image File Data from dicom_lst2.mat
%
load dicom_lst2.mat afiles ddirs idvr isz nimages pspc sn splcktc stxt;
%
idv = ~startsWith(stxt,'WIP');         % Exclude Philips T1rho values
idvr = idvr&idv;
%
n = sum(idvr);          % Number of possible series to analyze
ddirr = ddirs(idvr);    % Subdirectories for series
nfiles = nimages(idvr); % Numbers of files
afiler = afiles(idvr);  % T1rho files
isz = isz(idvr,:);      % Image sizes in pixels
sn = sn(idvr);          % Series numbers
pspc = pspc(idvr,:);    % Pixel sizes
spltr = splcktc(idvr);  % T1rho spin lock times
stxtr = stxt(idvr);     % T1rho series
stxtrc = cell(n,1);
for k = 1:n
   stxtrc{k} = ['Series ' int2str(sn(k)) ', ' stxtr{k}];
end
%
% Select Series to Analyze
%
[isel,ok] = listdlg('ListString',stxtrc,'SelectionMode','multiple', ...
                    'Name','T1rho Series','PromptString', ...
                    'Select T1rho Series to Analyze','ListSize', ...
                    [300 400]);
%
if ok<1
  warning(' *** WARNING in phantomr_plt1:  No series selected!');
  return
end
%
n = length(isel);       % Number of series to analyze
%
% Loop through the Series
%
for j = 1:n
%
% Get Series Specific Indices
%
   id = isel(j);
   ddir = ddirr{id};    % Subdirectory for series
%
   nfile = nfiles(id);  % Number of image files in this series
   splt = eval(['[' spltr{id} ']'])';  % T1rho spin lock times
   nsplt = size(splt,1);               % Number of spin lock times
   nsls = nfile/nsplt;  % Number of slices
   iszs = isz(id,:);    % Image size
   pspcs = pspc(id,:);  % Pixel size
   llim = 50/mean(pspcs);              % Lower limit on diameter
   ulim = 250/mean(pspcs);             % Upper limit on diameter
%
   sns = sn(id);        % Series number
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
   fnams = afiler{id};  % T1rho files
%
% Get Images and Maximum Scaled Image Value
%
   valmx = zeros(nfile,1);
   v = zeros([iszs nfile]);
%
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
% Initialize Variables for the Middle Slices of Scan
%
   lvl = zeros(nsl,1);  % Threshold levels for each slice
   rm = zeros(nsl,1);   % Threshold metric
   ctrs = zeros(nvial,2,nsl);     % Centroids of vials
   rmin = zeros(nvial,nsl);       % Minor axis lengths (minimum radii)
   vnum = int2str((1:nvial)');    % Vial number
   vmsk = false(prod(iszs),nvial,nsl);           % Mask for vials
   [xi,yi] = meshgrid(1:iszs(1),1:iszs(2));      % Pixel coordinates
   xi = xi(:);
   yi = yi(:);
   npix = size(xi,1);   % Number of pixels in a slice
   xyi = [xi yi];
%
   for k = 1:nsl
%
% Threshold Slice
%
      isl = istrt+k;    % Index to slice
      img = v0(:,:,isl);               % Spin lock time = 0 ms
      [lvl(k),rm(k)] = multithresh(img);
      bw = img>lvl(k);
%
% Get Properties of the Threshold Regions (Vials)
%
      rdat = regionprops('table',bw,'Centroid','MajorAxisLength', ...
                         'MinorAxisLength');
      rsize = size(rdat,1);
      if rsize>nvial
        idv = rdat.MajorAxisLength>llim&rdat.MinorAxisLength<ulim;
        rdat = rdat(idv,:);
        rsize = size(rdat,1);
        if rsize>nvial
          str = strel('square',2);
          bw = imerode(bw,str);        % Erode image - Use imopen?
          rdat = regionprops('table',bw,'Centroid','MajorAxisLength', ...
                             'MinorAxisLength');
          idv = rdat.MajorAxisLength>llim&rdat.MinorAxisLength<ulim;
          rdat = rdat(idv,:);
          rsize = size(rdat,1);
          if rsize~=nvial
            error([' *** ERROR in phantomr_pl1t:  Incorrect number', ...
                   ' of regions were identified!']);
          end
        elseif rsize<nvial
          error([' *** ERROR in phantomr_plt1:  Less than', ...
                           ' %i regions were identified!'],nvial);
        end
      elseif rsize<nvial
        error([' *** ERROR in phantomr_plt1:  Less than %i', ...
                         ' regions were identified!'],nvial);
      end
      ctrs(:,:,k) = rdat.Centroid;
      rmin(:,k) = rdat.MinorAxisLength/4;   % Reduce radius to avoid vial material
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
      if k==round(nsl/2)
        print('-dpsc2','-r600','-fillpage','-append',pfile);
      end
%
% Create a Mask
%
      for l = 1:nvial
         rc = xyi-repmat(ctrs(l,:,k),npix,1);
         rc = sum(rc.*rc,2);           % Radius squared
         vmsk(:,l,k) = rc<rmin(l,k).*rmin(l,k);
      end
   end
%
% Save Images and Masks
%
   mf = ['phantomr' snt '_plt1.mat'];  % MAT file name
   save(mf,'ctrs','ddir','fnams','idvr','istrt','iszs','lvl', ...
        'nfile','nsl','nsls','nsplt','nvial','rm','rmin','sns', ...
        'snt','splt','sltxt','v','v0','vmsk');
%
% Plot the Vial Masks for the Last Middle Slice Image
%
   figure;
   orient landscape;
   colormap gray;
%
   rmx = double(max(img(:)));
%
   for l = 1:nvial
%       subplot(2,3,l);
      imagesc(img);
      axis image;
      axis off;
      hold on;
      h = imagesc(reshape(vmsk(:,l,k)*rmx,iszs));
      set(h,'AlphaData',0.1);
%       title(['Vial ' int2str(l) ' Mask'],'FontSize',16, ...
%             'FontWeight','bold');
   end
%    
%    sgtitle({['Series ' snt ' Phantom Images']; ['Slice ', ...
%             int2str(isl)]},'FontSize',20,'FontWeight','bold');
%    
   title({['Series ' snt ' Phantom Images']; ['Slice ', ...
            int2str(isl)]; 'Vial 1 Mask'},'FontSize',16, ...
            'FontWeight','bold');
   %
   print('-dpsc2','-opengl','-r600','-fillpage','-append',pfile);
%
   if j<n
     pause;             % Review slice plots
     close all;
   end
%
end
%
return