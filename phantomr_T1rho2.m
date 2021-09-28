%#######################################################################
%
%                 * PHANTOM Reduced T1rho 2 Program *
%
%          M-File which reads the MRI phantom data and vial masks
%     to calculate the T1rho values in the vials.  For the reduced
%     radii segmentations.  All the pixels in each slice are used to
%     calculate the T1rho in each vial.  Also all the pixels in each
%     slice are used to calculate the T1rho in the vials.
%
%          The T1rho values and statistics are written to two MS-Excel
%     spreadsheets, phantomr?_T1rho2.xlsx and phantomrs?_T1rho2.xlsx.
%     ? is the series number.  Vial results are written to
%     phantomr?_T1rho2.xlsx and slice results are written to
%     phantomrs?_T1rho2.xlsx.
%
%          The data and results are saved into MAT files,
%     phantomr?_T1rho2.mat.
%
%     NOTES:  1.  Matlab MAT files dicom_lst.mat or dicom_lst2.mat and
%             phantomr?_plt2.mat must be in the current directory or
%             path.  ? is the series number.
%
%             2.  M-file exp_fun1.m must be in the current directory or
%             path.
%
%             3.  See phantomr_plt2.m for the segmentations of the
%             vials in the Series.
%
%     02-Jul-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
               'Algorithm','levenberg-marquardt','Jacobian','on', ...
               'UseParallel',true);
%
T1rho0 = [34; 45; 66];  % Initial T1rho values (maximum expected)
idx0 = [1; 3; 2];
idx0 = [idx0; flipud(idx0)];           % Index to initial vial T1rhos
fun = @exp_fun1;        % Exponential function
%
% Load Image File Data from dicom_lst.mat or dicom_lst2.mat
%
if exist('dicom_lst.mat','file')==2
  load dicom_lst.mat idvr sn stxt;
elseif exist('dicom_lst2.mat','file')==2
  load dicom_lst2.mat idvr sn stxt;
else
  error([' *** ERROR in phantomr_plt2:  dicom_lst.mat or ', ...
         'dicom_lst2.mat does not exist!']);
end
%
idv = ~startsWith(stxt,'WIP');         % Exclude Philips T1rho values
idvr = idvr&idv;
%
n = sum(idvr);          % Number of series to analyze
sn = sn(idvr);          % Series numbers
%
% Spreadsheet Labels
%
lbls = {'' 'N' 'Minimum' 'Maximum' 'Fit' 'Mean' 'SD' 'COV'};
tl = cell2table(lbls);
%
% Loop through the Series
%
for j = 1:n
%
% Get Series Number
%
   sns = sn(j);         % Series number
   snt = int2str(sns);  % Series number as text
%
% Load T1rho Masks
%
   matfile = ['phantomr' snt '_plt2.mat'];
   load(matfile,'istrts','iszs','nfile','nslm','nsls','nsplt', ...
                'nvial','splt','v','vmsk');
%
   npix = prod(iszs);   % Number of pixels in an image
%
   vnams = cellstr([repmat('Vial ',nvial,1) int2str((1:nvial)')]);   % Vial names as text
%
% Postscript Plot File Name
%
   pfile = ['phantomr' snt '_T1rho2.ps'];
%
% Read DICOM T1 Data for Slices and Spin Lock Times
%
   rimgs = cell(nslm,nvial);      % # of slices, # of vials
%
% Loop through Slices
%
   for k = 1:nslm       % Slice
%
% Slice Information
%
      slk = k+istrts;        % Slice number
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
         nf = nl+l;          % Index to file for this spin lock time
%
         img = v(:,:,nf);
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
   npxvs = zeros(nslm,nvial);          % Number of pixels in curvefits
   T1rhonlss = zeros(nslm,nvial);      % Nonlinear least squares time constants for each slice
   nlss_amp = zeros(nslm,nvial);       % Nonlinear least squares amplitudes for each slice
   nlss_sse = zeros(nslm,nvial);       % Nonlinear least squares square root of sum of squared errors for each slice
   exit_flags = zeros(nslm,nvial);     % Nonlinear least squares exit flags for each slice
%
   T1rhonlsp = cell(nslm,nvial);       % Nonlinear least squares time constants for each pixel
   nlsp_amp = cell(nslm,nvial);        % Nonlinear least squares amplitudes for each pixel
   nlsp_sse = cell(nslm,nvial);        % Nonlinear least squares square root of sum of squared errors for each pixel
   exit_flagp = cell(nslm,nvial);      % Nonlinear least squares exit flags for each pixel
%
   sT1rhovm = zeros(nslm,nvial);       % Mean for each slice
   sT1rhovmn = zeros(nslm,nvial);      % Minimum for each slice
   sT1rhovmx = zeros(nslm,nvial);      % Maximum for each slice
   sT1rhovsd = zeros(nslm,nvial);      % SD for each slice
%
   imgs = zeros(npix,nslm);            % T1rho values by pixel by slice
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
      for l = 1:nslm
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
   for k = 1:nslm       % Slice
%
      figure;
      orient landscape;
      imagesc(reshape(imgs(:,k),iszs));
      axis image;
      colormap jet;
      caxis([26 90]);
      colorbar;
      title({['Series ' snt ' Phantom T1\rho']; ['Slice ' ...
            int2str(k+istrts)]; splt_txt},'FontSize',16, ...
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
   xlsfile = ['phantomr' snt '_T1rho2.xlsx'];
%
   t = table(npxv,T1rhovmn,T1rhovmx,T1rhonls,T1rhovm,T1rhovsd, ...
            T1rhovcov,'RowNames',vnams,'VariableNames',lbls(2:8), ...
            'DimensionNames',{'Vial #','Variables'});
   writetable(t,xlsfile,'WriteMode','replacefile','WriteRowNames',true);
%
% Write Vial by Slice T1rho to MS-Excel Spreadsheet
%
   xlsfile = ['phantomrs' snt '_T1rho2.xlsx'];
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
   sl_lbls = [repmat('Slice ',nslm,1) int2str((istrts+1:istrts+nslm)')];
   sl_lbls = cellstr(sl_lbls);
   sl_lbls = strrep(sl_lbls,'  ',' ');
%
   dat0 = NaN(nslm,nvial-1);
   dat = [npxvs sT1rhovmn sT1rhovmx T1rhonlss sT1rhovm sT1rhovsd, ...
         sT1rhovcov dat0];
   nhdr = nhdr-1;       % With row names
   idv = [1:6:nhdr 2:6:nhdr 3:6:nhdr 4:6:nhdr 5:6:nhdr 6:6:nhdr];
   dat = dat(:,idv);
   ts = array2table(dat,'RowNames',sl_lbls);
   writetable(ts,xlsfile,'WriteMode','append','WriteVariableNames', ...
              false,'WriteRowNames',true);
%
% Save Data to MAT Files
%
   matfiles = ['phantomr' snt '_T1rho2.mat'];
   save(matfiles,'T1rhonls','T1rhonlsp','T1rhonlss','exit_flag', ...
        'exit_flagp','exit_flags','idx0','imgs','istrts','iszs', ...
        'lbls','lblss','nfile','nls_amp','nls_sse','nlsp_amp');
   save(matfiles,'-append','nlsp_sse','nlss_amp','nlss_sse','npix', ...
        'npxv','npxvs','nslm','nsls','nsplt','pfile','rimgs', ...
        'sT1rhovm','sT1rhovmn','sT1rhovmx','sT1rhovsd','sl_lbls');
   save(matfiles,'-append','sns','snt','splt','t','tlbls','ts', ...
        'vnams','vnams_hdr');
%
% Close Plot Windows for this Series
%
   close all;
%
end
%
return