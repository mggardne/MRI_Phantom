%#######################################################################
%
%                   * PHANTOM Rho T1rho 1 Program *
%
%          M-File which reads the MRI phantom data and vial mask to
%     calculate the T1rho values in the one large vial.  T1rho values
%     are calculated for all the pixels in the vial for each slice.
%     Also, all the pixels in the vial for each slice are used to
%     calculate single T1rho values for the vials.
%
%          The T1rho values and statistics are written to two MS-Excel
%     spreadsheets, phantomr?_T1rho1.xlsx and phantomrs?_T1rho1.xlsx.
%     ? is the series number.  Vial results are written to
%     phantomr?_T1rho1.xlsx and slice results are written to
%     phantomrs?_T1rho1.xlsx.
%
%          The plots are saved into Postscript files,
%     phantomr?_T1rho1.ps.  ? is the series number.
%
%          The data and results are saved into MAT files,
%     phantomr?_T1rho1.mat.  ? is the series number.
%
%     NOTES:  1.  Matlab MAT files phantomr?_plt1.mat must be in the
%             current directory or path.  ? is the series number.
%
%             2.  M-file exp_fun1.m must be in the current directory or
%             path.
%
%             3.  See phantomr_plt1.m for the segmentations of the
%             vial in the Series.
%
%     30-Nov-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
               'Algorithm','levenberg-marquardt','Jacobian','on', ...
               'UseParallel',true);
%
% T1rho0 = [34; 45; 66];  % Initial T1rho values (maximum expected)
% idx0 = [1; 3; 2];
% idx0 = [idx0; flipud(idx0)];           % Index to initial vial T1rhos
T1rho0 = 45;
idx0 = 1;
fun = @exp_fun1;        % Exponential function
%
% Get T1rho Phantom Segmentation MAT Files
%
fmats = dir('phantomr*_plt1.mat');
%
if size(fmats,1)<1
  warning([' *** WARNING in phantomr_T1rho1:  No segmentation', ...
           ' MAT files found in current directory!']);
  return
end
%
fmats = {fmats.name}';
% fmats = fmats([1 3]);   % Just do sagittal plane T1rho
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
% Load T1rho Masks and Series Data
%
   matfile = fmats{o};
   load(matfile,'istrt','iszs','nfile','nsl','nsls','nsplt', ...
                'nvial','sns','snt','splt','v','vmsk');
%
   npix = prod(iszs);   % Number of pixels in an image
%
   vnams = cellstr([repmat('Vial ',nvial,1) int2str((1:nvial)')]);   % Vial names as text
%
% Postscript Plot File Name
%
   pfile = ['phantomr' snt '_T1rho1.ps'];
%
% Read DICOM T1 Data for Slices and Spin Lock Times
%
   rimgs = cell(nsl,nvial);            % # of slices, # of vials
%
% Loop through Slices
%
   for k = 1:nsl        % Slice
%
% Slice Information
%
      slk = k+istrt;         % Slice number
      sll = int2str(slk);    % Slice number as letters
      nl = (slk-1)*nsplt;    % Index to last file in previous slice
%
% Loop through Spin Lock Times
%
      rimgc = cell(nsplt,nvial);  % # of spin lock times, # of vials
%
      for l = 1:nsplt   % Loop through spin lock times
%
% Read T1 Slice Image
%
         ni = nl+l;          % Index to file for this spin lock time
%
         img = v(:,:,ni);
%
         for m = 1:nvial
            rimgc{l,m} = img(vmsk(:,m,k))';
         end
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
   xdat = [splt ones(nsplt,1)];        % Spin lock times in ms and exponential amplitude
%
   npxv = zeros(nvial,1);              % Number of pixels in curvefits
   T1rhonls = zeros(nvial,1);          % Nonlinear least squares time constants
   nls_amp = zeros(nvial,1);           % Nonlinear least squares amplitudes
   nls_sse = zeros(nvial,1);           % Nonlinear least squares square root of sum of squared errors
   exit_flag = zeros(nvial,1);         % Nonlinear least squares exit flags
%
   npxvs = zeros(nsl,nvial);           % Number of pixels in curvefits
   T1rhonlss = zeros(nsl,nvial);       % Nonlinear least squares time constants for each slice
   nlss_amp = zeros(nsl,nvial);        % Nonlinear least squares amplitudes for each slice
   nlss_sse = zeros(nsl,nvial);        % Nonlinear least squares square root of sum of squared errors for each slice
   exit_flags = zeros(nsl,nvial);      % Nonlinear least squares exit flags for each slice
%
   T1rhonlsp = cell(nsl,nvial);        % Nonlinear least squares time constants for each pixel
   nlsp_amp = cell(nsl,nvial);         % Nonlinear least squares amplitudes for each pixel
   nlsp_sse = cell(nsl,nvial);         % Nonlinear least squares square root of sum of squared errors for each pixel
   exit_flagp = cell(nsl,nvial);       % Nonlinear least squares exit flags for each pixel
%
   sT1rhovm = zeros(nsl,nvial);        % Mean for each slice
   sT1rhovmn = zeros(nsl,nvial);       % Minimum for each slice
   sT1rhovmx = zeros(nsl,nvial);       % Maximum for each slice
   sT1rhovsd = zeros(nsl,nvial);       % SD for each slice
%
   imgs = zeros(npix,nsl);             % T1rho values by pixel by slice
%
% Calculate T1rho for Each Vial
%
   for k = 1:nvial
      rimgv = cat(2,rimgs{:,k});
      npxv(k) = size(rimgv,2);
      amp0 = median(rimgv(1,:));
      rp0 = [amp0; T1rho0(idx0(k))];
%
      sltk = repmat(splt,npxv(k),1);
      rimgv = rimgv(:);
%
% Nonlinear Least Squares Exponential Fit to Get T1rho Values
%
      [rp,~,~,eflag] = lsqcurvefit(fun,rp0,sltk,rimgv,[],[],opt);
%
      T1rhonls(k) = rp(2);
      nls_amp(k) = rp(1);
      d = exp_fun1(rp,sltk)-rimgv;
      nls_sse(k) = sqrt(d'*d);
      exit_flag(k) = eflag;
%
% Calculate T1rho for Each Vial in Each Slice
%
      for l = 1:nsl
         rimgv = rimgs{l,k};
         npxvs(l,k) = size(rimgv,2);
%
         sltk = repmat(splt,npxvs(l,k),1);
         rimgvs = rimgv(:);
%
         [rp,~,~,eflag] = lsqcurvefit(fun,rp0,sltk,rimgvs,[],[],opt);
%
         T1rhonlss(l,k) = rp(2);
         nlss_amp(l,k) = rp(1);
         d = exp_fun1(rp,sltk)-rimgvs;
         nlss_sse(l,k) = sqrt(d'*d);
         exit_flags(l,k) = eflag;
%
         T1rp = zeros(npxvs(l,k),1);
         ampp = zeros(npxvs(l,k),1);
         ssep = zeros(npxvs(l,k),1);
         eflagp = zeros(npxvs(l,k),1);
%
% Calculate T1rho for Each Pixel in Each Vial in Each Slice
%
         parfor m = 1:npxvs(l,k)
%
            rimgm = rimgv(:,m);
            [rp,~,~,eflag] = lsqcurvefit(fun,rp0,splt,rimgm,[],[],opt);
            T1rp(m) = rp(2);
            ampp(m) = rp(1);
            d = exp_fun1(rp,splt)-rimgm;
            ssep(m) = sqrt(d'*d);
            eflagp(m) = eflag;
%
         end
%
         T1rhonlsp{l,k} = T1rp;
         nlsp_amp{l,k} = ampp;
         nlsp_sse{l,k} = ssep;
         exit_flagp{l,k} = eflagp;
%
% Calculate Statistics for Each Slice
%
         sT1rhovm(l,k) = mean(T1rp);
         sT1rhovmn(l,k) = min(T1rp);
         sT1rhovmx(l,k) = max(T1rp);
         sT1rhovsd(l,k) = std(T1rp);
         img = zeros(npix,1);
         img(vmsk(:,k,l)) = T1rp;
         imgs(:,l) = imgs(:,l)+img;
%
      end
   end
%
   sT1rhovcov = 100*sT1rhovsd./sT1rhovm;    % Calculate COV (%)
%
% Plot T1rho Values for Each Slice
%
   splt_txt = sprintf('%i, ',splt);
   splt_txt = ['Spin Lock Times = ' splt_txt(1:end-2) ' ms'];
%
   for k = 1:nsl        % Slice
%
      figure;
      orient landscape;
      imagesc(reshape(imgs(:,k),iszs));
      axis image;
      colormap jet;
      caxis([30 60]);
      colorbar;
      title({['Series ' snt ' Phantom T1\rho']; ['Slice ' ...
            int2str(k+istrt)]; splt_txt},'FontSize',16, ...
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
   T1rhovm = zeros(nvial,1);
   T1rhovmn = zeros(nvial,1);
   T1rhovmx = zeros(nvial,1);
   T1rhovsd = zeros(nvial,1);
%
   for k = 1:nvial
      T1rhov = cell2mat(T1rhonlsp(:,k));
      T1rhovm(k) = mean(T1rhov);
      T1rhovmn(k) = min(T1rhov);
      T1rhovmx(k) = max(T1rhov);
      T1rhovsd(k) = std(T1rhov);
   end
%
   T1rhovcov = 100*T1rhovsd./T1rhovm;
%
% Write Vial T1rho to MS-Excel Spreadsheet
%
   xlsfile = ['phantomr' snt '_T1rho1.xlsx'];
%
   t = table(npxv,T1rhovmn,T1rhovmx,T1rhonls,T1rhovm,T1rhovsd, ...
            T1rhovcov,'RowNames',vnams,'VariableNames',lbls(2:8), ...
            'DimensionNames',{'Vial #','Variables'});
   writetable(t,xlsfile,'WriteMode','replacefile','WriteRowNames',true);
%
% Write Vial by Slice T1rho to MS-Excel Spreadsheet
%
   xlsfile = ['phantomrs' snt '_T1rho1.xlsx'];
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
   dat = [npxvs sT1rhovmn sT1rhovmx T1rhonlss sT1rhovm sT1rhovsd, ...
          sT1rhovcov dat0];
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
   matfiles = ['phantomr' snt '_T1rho1.mat'];
   save(matfiles,'T1rhonls','T1rhonlsp','T1rhonlss','exit_flag', ...
        'exit_flagp','exit_flags','idx0','imgs','istrt','iszs', ...
        'lbls','lblss','nfile','nls_amp','nls_sse','nlsp_amp');
   save(matfiles,'-append','nlsp_sse','nlss_amp','nlss_sse','npix', ...
        'npxv','npxvs','nsl','nsls','nsplt','pfile','rimgs', ...
        'sT1rhovm','sT1rhovmn','sT1rhovmx','sT1rhovsd','sl_lbls');
   save(matfiles,'-append','sns','snt','splt','t','tlbls','ts', ...
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