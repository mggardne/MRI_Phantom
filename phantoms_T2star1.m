%#######################################################################
%
%                  * PHANTOM Star T2STAR 1 Program *
%
%          M-File which reads the MRI phantom data and vial mask to
%     calculate the T2* values in the one large vial.  T2* values are
%     calculated for all the pixels in the vial for each slice.  Also
%     all the pixels in the vial for each slice are used to calculate
%     single T2* values for the vials.
%
%          The T2* values and statistics are written to two MS-Excel
%     spreadsheets, phantoms?_T2star1.xlsx and phantomss?_T2star1.xlsx.
%     ? is the series number.  Vial results are written to
%     phantoms?_T2star1.xlsx and slice results are written to
%     phantomss?_T2star1.xlsx.
%
%          The plots are saved into Postscript files,
%     phantoms?_T2star1.ps.  ? is the series number.
%
%          The data and results are saved into MAT files,
%     phantoms?_T2star1.mat.  ? is the series number.
%
%     NOTES:  1.  Matlab MAT files phantoms?_plt1.mat must be in the
%             current directory or path.  ? is the series number.
%
%             2.  M-file exp_fun1.m must be in the current directory or
%             path.
%
%             3.  See phantoms_plt1.m for the segmentations of the
%             vials in the Series.
%
%     28-Oct-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
               'Algorithm','levenberg-marquardt','Jacobian','on', ...
               'UseParallel',true);
%
% T2star0 = [22; 11.5; 3];     % Initial T2* values (~average values)
% idx0 = (1:3)';
% idx0 = [idx0; flipud(idx0)]; % Index to initial vial T2*
% T2star0 = 11.5;         % Initial T2* values (~average value)
T2star0 = 25;           % Initial T2* values (~average value)
idx0 = 1;
fun = @exp_fun1;        % Exponential function
%
% Get T2* Phantom Segmentation MAT Files
%
fmats = dir('phantoms*_plt1.mat');
%
if size(fmats,1)<1
  warning([' *** WARNING in phantoms_T2star1:  No segmentation', ...
           ' MAT files found in current directory!']);
  return
end
%
fmats = {fmats.name}';
nf = size(fmats,1);     % Number of MAT files to process
%
% Spreadsheet Labels
%
lbls = {'' 'N' 'Minimum' 'Maximum' 'Fit' 'Mean' 'SD' 'COV'};
tl = cell2table(lbls);
%
% Loop through MAT Files
%
for o = 1:nf
%
% Load T2* Masks and Series Data
%
   matfile = fmats{o};
   load(matfile,'afiles','ddirr','etns','etxt','istrt','iszs','n', ...
                'nfile','nsl','nvial','sns','snt','v','vmsk');
%
   npix = prod(iszs);   % Number of pixels in an image
%
   vnams = cellstr([repmat('Vial ',nvial,1) int2str((1:nvial)')]);   % Vial names as text
%
% Postscript Plot File Name
%
   pfile = ['phantoms' snt '_T2star1.ps'];
%
% Get DICOM T2 Data for Slices
%
   rimgs = cell(nsl,nvial);  % # of slices, # of vials
%
% Loop through Slices
%
   for k = 1:nsl        % Slice
%
% Slice Information
%
      slk = k+istrt;         % Slice number
      sll = int2str(slk);    % Slice number as letters
%
% Loop through Echo Times
%
      rimgc = cell(n,nvial);           % # of echo times, # of vials
%
      for l = 1:n       % Loop through echo times
%
         img = squeeze(v(:,:,slk,l));  % T2 data for slice and echo time
%
         for m = 1:nvial   % Loop through vials
            rimgc{l,m} = img(vmsk(:,m,k))';
         end
%
      end
%
      for m = 1:nvial
         rimgs{k,m} = cat(1,rimgc{:,m});
      end
%
   end
%
   clear img rimgc;
%
% Set Up Arrays for Loop
%
   xdat = [etns ones(n,1)];       % Spin lock times in ms and exponential amplitude
%
   npxv = zeros(nvial,1);         % Number of pixels in curvefits
   T2starnls = zeros(nvial,1);    % Nonlinear least squares time constants
   nls_amp = zeros(nvial,1);      % Nonlinear least squares amplitudes
   nls_sse = zeros(nvial,1);      % Nonlinear least squares square root of sum of squared errors
   exit_flag = zeros(nvial,1);    % Nonlinear least squares exit flags
%
   npxvs = zeros(nsl,nvial);      % Number of pixels in curvefits for each slice
   T2starnlss = zeros(nsl,nvial); % Nonlinear least squares time constants for each slice
   nlss_amp = zeros(nsl,nvial);   % Nonlinear least squares amplitudes for each slice
   nlss_sse = zeros(nsl,nvial);   % Nonlinear least squares square root of sum of squared errors for each slice
   exit_flags = zeros(nsl,nvial); % Nonlinear least squares exit flags for each slice
%
   T2starnlsp = cell(nsl,nvial);  % Nonlinear least squares time constants for each pixel
   nlsp_amp = cell(nsl,nvial);    % Nonlinear least squares amplitudes for each pixel
   nlsp_sse = cell(nsl,nvial);    % Nonlinear least squares square root of sum of squared errors for each pixel
   exit_flagp = cell(nsl,nvial);  % Nonlinear least squares exit flags for each pixel
%
   sT2starvm = zeros(nsl,nvial);  % Mean for each slice
   sT2starvmn = zeros(nsl,nvial); % Minimum for each slice
   sT2starvmx = zeros(nsl,nvial); % Maximum for each slice
   sT2starvsd = zeros(nsl,nvial); % SD for each slice
%
   imgs = zeros(npix,nsl);        % T2* values by pixel by slice
%
% Calculate T2* for Each Vial
%
   for k = 1:nvial
      rimgv = cat(2,rimgs{:,k});
      npxv(k) = size(rimgv,2);
      amp0 = median(rimgv(1,:));
      rp0 = [amp0; T2star0(idx0(k))];
%
      etk = repmat(etns,npxv(k),1);
      rimgv = rimgv(:);
%
% Nonlinear Least Squares Exponential Fit to Get T2* Values
%
      [rp,~,~,eflag] = lsqcurvefit(fun,rp0,etk,rimgv,[],[],opt);
%
      T2starnls(k) = rp(2);
      nls_amp(k) = rp(1);
      d = exp_fun1(rp,etk)-rimgv;
      nls_sse(k) = sqrt(d'*d);
      exit_flag(k) = eflag;
%
% Calculate T2* for Each Vial in Each Slice
%
      for l = 1:nsl
         rimgv = rimgs{l,k};
         npxvs(l,k) = size(rimgv,2);
%
         etk = repmat(etns,npxvs(l,k),1);
         rimgvs = rimgv(:);
%
         [rp,~,~,eflag] = lsqcurvefit(fun,rp0,etk,rimgvs,[],[],opt);
%
         T2starnlss(l,k) = rp(2);
         nlss_amp(l,k) = rp(1);
         d = exp_fun1(rp,etk)-rimgvs;
         nlss_sse(l,k) = sqrt(d'*d);
         exit_flags(l,k) = eflag;
%
         T2rp = zeros(npxvs(l,k),1);
         ampp = zeros(npxvs(l,k),1);
         ssep = zeros(npxvs(l,k),1);
         eflagp = zeros(npxvs(l,k),1);
%
% Calculate T2* for Each Pixel in Each Vial in Each Slice
%
         parfor m = 1:npxvs(l,k)
%      
            rimgm = rimgv(:,m);
            [rp,~,~,eflag] = lsqcurvefit(fun,rp0,etns,rimgm,[],[],opt);
            T2rp(m) = rp(2);
            ampp(m) = rp(1);
            d = exp_fun1(rp,etns)-rimgm;
            ssep(m) = sqrt(d'*d);
            eflagp(m) = eflag;
%
         end
%
         T2starnlsp{l,k} = T2rp;
         nlsp_amp{l,k} = ampp;
         nlsp_sse{l,k} = ssep;
         exit_flagp{l,k} = eflagp;
%
% Calculate Statistics for Each Slice
%
         sT2starvm(l,k) = mean(T2rp);
         sT2starvmn(l,k) = min(T2rp);
         sT2starvmx(l,k) = max(T2rp);
         sT2starvsd(l,k) = std(T2rp);
         img = zeros(npix,1);
         img(vmsk(:,k,l)) = T2rp;
         imgs(:,l) = imgs(:,l)+img;
%
      end
   end
%
   sT2starvcov = 100*sT2starvsd./sT2starvm; % Calculate COV (%)
%
% Plot T2* Values for Each Slice
%
   for k = 1:nsl        % Slice
%
      figure;
      orient landscape;
      imagesc(reshape(imgs(:,k),iszs));
      axis image;
      colormap jet;
%      if strcmp(snt,'1601')            % 27 Aug 2021 series only
%        caxis([0 60]);
%      else
%        caxis([0 20]);
        caxis([15 35]);
%      end
      colorbar;
      title({['Series ' snt ' Phantom T2*']; ['Slice ' ...
            int2str(k+istrt)]; etxt},'FontSize',16, ...
            'FontWeight','bold');
      if k==1
        print('-dpsc2','-r600','-fillpage',pfile);
      else
        print('-dpsc2','-r600','-fillpage','-append',pfile);
      end
%
   end
%
% Calculate Vial Statistics
%
   T2starvm = zeros(nvial,1);
   T2starvmn = zeros(nvial,1);
   T2starvmx = zeros(nvial,1);
   T2starvsd = zeros(nvial,1);
%
   for k = 1:nvial
      T2starv = cell2mat(T2starnlsp(:,k));
      T2starvm(k) = mean(T2starv);
      T2starvmn(k) = min(T2starv);
      T2starvmx(k) = max(T2starv);
      T2starvsd(k) = std(T2starv);
   end
%
   T2starvcov = 100*T2starvsd./T2starvm;
%
% Write Vial T2* to MS-Excel Spreadsheet
%
   xlsfile = ['phantoms' snt '_T2star1.xlsx'];
%
   t = table(npxv,T2starvmn,T2starvmx,T2starnls,T2starvm,T2starvsd, ...
            T2starvcov,'RowNames',vnams,'VariableNames',lbls(2:8), ...
            'DimensionNames',{'Vial #','Variables'});
   writetable(t,xlsfile,'WriteMode','replacefile','WriteRowNames',true);
%
% Write Vial by Slice T2* to MS-Excel Spreadsheet
%
   xlsfile = ['phantomss' snt '_T2star1.xlsx'];
%
   lblss = repmat(lbls,1,nvial);
   nhdr = size(lblss,2);
   vnams_hdr = cell(1,nhdr);
   vnams_hdr{1} = 'Slice #';
   vnams_hdr(5:8:nhdr) = vnams;
   tlbls = cell2table([vnams_hdr; lblss]);
   writetable(tlbls,xlsfile,'WriteMode','replacefile', ...
              'WriteVariableNames',false);
%
   sl_lbls = [repmat('Slice ',nsl,1) int2str((istrt+1:istrt+nsl)')];
   sl_lbls = cellstr(sl_lbls);
   sl_lbls = strrep(sl_lbls,'  ',' ');
%
   dat0 = NaN(nsl,nvial-1);
   dat = [npxvs sT2starvmn sT2starvmx T2starnlss sT2starvm sT2starvsd, ...
          sT2starvcov dat0];
   nhdr = nhdr-1;       % With row names
   idv = [1:7:nhdr 2:7:nhdr 3:7:nhdr 4:7:nhdr 5:7:nhdr 6:7:nhdr ...
          7:7:nhdr];
   dat = dat(:,idv);
   ts = array2table(dat,'RowNames',sl_lbls);
   writetable(ts,xlsfile,'WriteMode','append','WriteVariableNames', ...
              false,'WriteRowNames',true);
%
% Save Data to MAT Files
%
   matfiles = ['phantoms' snt '_T2star1.mat'];
   save(matfiles,'T2starnls','T2starnlsp','T2starnlss','etns', ...
        'exit_flag','exit_flagp','exit_flags','idx0','imgs', ...
        'istrt','iszs','lbls','lblss','n','nfile','nls_amp');
   save(matfiles,'-append','nls_sse','nlsp_amp','nlsp_sse','nlss_amp', ...
        'nlss_sse','npix','npxv','npxvs','nsl','pfile','rimgs', ...
        'sT2starvm','sT2starvmn','sT2starvmx','sT2starvsd','sl_lbls');
   save(matfiles,'-append','sns','snt','t','tlbls','ts', ...
        'vnams','vnams_hdr');
%
% Close Plot Windows for this MAT File
%
   if o<nf
     pause;
     close all;
   end
%
end
%
return