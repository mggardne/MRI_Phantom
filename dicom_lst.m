%#######################################################################
%
%                        * DICOM LiST Program *
%
%          M-File which reads a DICOMDIR file and collects information
%      about the MRI series and DICOM files.
%
%     NOTES:  1.  Reads Philips MRI DICOMDIR files.  May not work with
%             other types of images.
%
%             2.  Program also reads and collects information from the
%             first DICOM file in each series.
%
%             3.  A table of information is displayed to the screen,
%             written to a MS-Excel spreadsheet and to a Matlab MAT
%             file.  The spreadsheet and MAT file are written in the
%             same directory as the DICOMDIR file.  The MAT file
%             contains additional variables and all the file names in
%             each series and patient position in MRI coordinates.
%
%             4.  The program traps for some, but not all parameters in
%             the DICOM file header.
%
%     17-Aug-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Directory with DICOM Images
%
[ddirfile,ddir] = uigetfile({'DICOMDIR', 'DICOMDIR File'; '*.*', ...
                             'All Files (*.*)'},'Pick a DICOMDIR file');
ddirfile = fullfile(ddir,ddirfile);
%
% Check for DICOMDIR File
%
if exist(ddirfile,'file')
%
% Get File Names and Acquired Pixel Spacing
%
  d0 = images.dicom.parseDICOMDIR(ddirfile);
  ns = size(d0.Patients.Studies.Series,2);
%
  nimages = zeros(ns,1);               % Number of DICOM files
  arow = zeros(ns,1);                  % Image row size
  acol = zeros(ns,1);                  % Image column size
  apixel = zeros(ns,2);                % Pixel spacing
  afile1 = cell(ns,1);                 % First DICOM file in series
  afile2 = cell(ns,1);                 % Last DICOM file in series
  afiles = cell(ns,1);                 % All DICOM file in series
  im_type = cell(ns,1);                % Image type in series
%  
  for k = 1:ns
     nimages(k) = size(d0.Patients.Studies.Series(k).Images,2);
     if nimages(k)>0
       n = nimages(k);
       arow(k) = d0.Patients.Studies.Series(k).Images(1).Payload.Rows;
       acol(k) = d0.Patients.Studies.Series(k).Images(1).Payload. ...
                 Columns;
       apixel(k,:) = d0.Patients.Studies.Series(k).Images(1). ...
                     Payload.PixelSpacing';
       im_type{k} = d0.Patients.Studies.Series(k).Images(1). ...
                    Payload.ImageType;
       afile1{k} = d0.Patients.Studies.Series(k).Images(1). ...
                     Payload.ReferencedFileID;
       afile2{k} = d0.Patients.Studies.Series(k).Images(n). ...
                     Payload.ReferencedFileID;
       afiles{k} = cell(n,1);
       for l = 1:n
          afiles{k}{l} = d0.Patients.Studies.Series(k).Images(l). ...
                         Payload.ReferencedFileID;
       end
     end
  end
else
  uiwait(msgbox('DICOMDIR file not found!','Information', ...
                'warn','modal'));
end
%
% Get DICOM Information from Files
%
idx = nimages>0;        % Actual series and not just space between series
nseries = sum(idx);
id = find(idx);
%
adur = zeros(nseries,2);               % Acquisition duration (s)
adurs = cell(nseries,1);               % Acquisition duration string (mm:ss)
etn = cell(nseries,1);                 % Echo times as numbers for UTE T2 star sequences
ets = cell(nseries,1);                 % Echo times as a string for UTE T2 star sequences
psz = zeros(nseries,2);                % Rows and Columns
isz = zeros(nseries,2);                % Width and Height
sthk = zeros(nseries,1);               % Slice Thickness
sspc = zeros(nseries,1);               % Spacing Between Slices
pos = zeros(nseries,3);
pspc = zeros(nseries,2);               % Pixel Spacing
ptxt = cell(nseries,1);                % Protocol Name
stxt = cell(nseries,1);                % Series Description
sl = zeros(nseries,1);                 % Rescale Slope
rinterc = zeros(nseries,1);            % Rescale Intercept
splt = zeros(nseries,1);               % Spin lock time (Trigger Time)
sn = zeros(nseries,1);                 % Series number
sdat = datetime(zeros(nseries,1),0,0,'Format','dd-MMM-yyyy HH:mm:ss');
%
for k = 1:nseries
%
   l = id(k);
%
   fnam = afiles{l}{1};
   fnam = fullfile(ddir,fnam);
   if exist(fnam,'file')
     info = dicominfo(fnam);
%
     if isfield(info,'AcquisitionDuration')
       adur(k) = info.AcquisitionDuration;
       adurs{k} = char(duration(seconds(adur(k)),'Format','mm:ss'));
     else
       adurs{k} = '00:00';
     end
     if isfield(info,'Rows')
       psz(k,:) = [info.Rows info.Columns];
       isz(k,:) = [info.Width info.Height];
     end
     if isfield(info,'SliceThickness')
       sthk(k) = info.SliceThickness;
       pspc(k,:) = info.PixelSpacing';
       sspc(k) = info.SpacingBetweenSlices;
     end
     if isfield(info,'ImagePositionPatient')
       pos(k,:) = info.ImagePositionPatient';
     end
     if isfield(info,'RescaleSlope')
       sl(k) = info.RescaleSlope;
       rinterc(k) = info.RescaleIntercept;
     end
     ptxt{k} = info.ProtocolName;
     stxt{k} = info.SeriesDescription;
     sn(k) = info.SeriesNumber;
     if isfield(info,'TriggerTime')
       splt(k) = info.TriggerTime;
     end
     if isfield(info,'EchoTrainLength')
       n = info.EchoTrainLength;
       if n>nimages(l)
         n = nimages(l);
       end
       ets{k} = info.Private_2001_1025;
       et = NaN(n,1); % Echo times for this sequence
       et(1) = info.EchoTime;
       for m = 2:n
          fnam = afiles{l}{m};
          fnam = fullfile(ddir,fnam);
          if exist(fnam,'file')
            info1 = dicominfo(fnam);
            et(m) = info1.EchoTime;
          else
            break;
          end
       end
       etu = unique(et);
       idnnan = ~isnan(etu);
       etn{k} = etu(idnnan);
       if size(etn{k},1)==1&&n>1
         idslash = strfind(ets{k},'/');
         if ~isempty(idslash)
           et1 = eval(ets{k}(1:idslash-1));
           det1 = eval(ets{k}(idslash+1:end));
           etn{k} = (et1:det1:(n-1)*det1+et1)';
         else
           etn{k} = eval(ets{k});
         end
       end
       ets{k} = [ets{k} '(' int2str(n) ')'];
     end
     dat = info.SeriesDate;
     tim = info.SeriesTime;
     sdat(k) = datetime([dat tim],'InputFormat','yyyyMMddHHmmss.SSSSS');
   end
%
end
%
% Get Spin Lock Times from Series Description
%
idvr = contains(stxt,' SL ');
%
splcktc = cellstr(repmat('None',nseries,1));     % Cell array for all series
%
if any(idvr)
  splckt = extractBetween(stxt(idvr),' SL ','ms');    % Spin lock times as text
  splckt = strrep(splckt,'_',',');
  splcktc(idvr) = splckt;
end
%
% Include UTE T2star Echo Times
%
idvu = contains(splcktc,'None')&~cellfun(@isempty,ets);
%
if any(idvu)
  splcktc(idvu) = ets(idvu);
end
%
% Put Series Data into a Table and Get Column Names
%
colnams = {'Series#','SeriesDateTime','SeriesDescription', ...
           'ProtocolName','ImageType','AcquisitionDuration', ...
           'SpinLockTimes','#ofFiles','FirstFile','LastFile','Rows', ...
           'Columns','PixelX','PixelY','SliceThickness', ...
           'SpacingBetweenSlices','ImageRows','ImageColumns', ...
           'Width','Height','ImagePixelX','ImagePixelY', ...
           'RescaleSlope','RescaleIntercept'};
t0 = table(sn,sdat,string(stxt),string(ptxt),string(im_type(idx)), ...
           char(adurs),string(splcktc),nimages(idx), ...
           string(afile1(idx)),string(afile2(idx)),arow(idx), ...
           acol(idx),round(apixel(idx,1),3),round(apixel(idx,2),3), ...
           sthk,sspc,psz(:,1),psz(:,2),isz(:,1),isz(:,2), ...
           round(pspc(:,1),3),round(pspc(:,2),3),sl,rinterc, ...
           'VariableNames',colnams)
%
% Write Table to Spreadsheet
%
xlsnam = fullfile(ddir,'dicom_lst.xlsx');
writetable(t0,xlsnam);
%
% Save MAT File
%
clear d0 dat det1 et1 etu idnnan idslash info info1 k l m n tim;
matnam = fullfile(ddir,'dicom_lst.mat');
save(matnam);
%
return