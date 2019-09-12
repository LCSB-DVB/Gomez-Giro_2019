%% Golgi Morphometrics Opera 3D 
clear
clc

%% User inputs

Barcodes = {'S:\OperaQEHS\OperaDB\Gemma Giro\GG_20180404_60X_Endo_GM130_2\GG_20180404_60X_Endo_GM130_2', ...
    'S:\OperaQEHS\OperaDB\Gemma Giro\GG_20180404_60X_Endo_GM130_3\GG_20180404_60X_Endo_GM130_3'};


%% Run documentation

FolderThisAnalysis = ['S:\HCS_Platform\Data\JavierJarazo\For_Gemma\Epithelial_Golgi_', datestr(clock,'yyyymmdd_HH_MM_SS')];
mkdir(FolderThisAnalysis)
FileNameShort = mfilename;
newbackup = sprintf('%s_log.m',[FolderThisAnalysis, '\', FileNameShort]);
FileNameAndLocation=[mfilename('fullpath')];

currentfile = strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);
PreviewPath = [FolderThisAnalysis, '\Previews'];
mkdir(PreviewPath)
VersionMatlab = version;
save([FolderThisAnalysis, '\', 'MatlabVersion.mat'], 'VersionMatlab');

f_LogDependencies(FileNameShort, FolderThisAnalysis); 
ObjectsAllFields = {};

% delete(gcp('nocreate'))
% parpool()

for n = 1:size(Barcodes,2)


        Barcode = Barcodes{n};
        InfoTable = f_InfoTable(Barcodes{n}, 'S:\OperaQEHS\OperaDB\Gemma Giro\GG_20180404_60X_Endo_GM130_2\Settings\GM130.lay');
        
        parfor l = 1:height(InfoTable)
        %for l = 1:2
            cube = readflexcube(InfoTable.files{l}, 'PlaneCount', 6); % Read 4-D image cube
            ch1 = cube.data(:, :, :, 1); %Hoechst vol(ch1)
            ch2 = cube.data(:, :, :, 2); % GM130 vol(ch2)
            
               
            InfoTableThis = InfoTable(l,:);
            ObjectsAllFields{n,l} = ImageAnalysisGolgi3D_PerField(ch1,ch2,InfoTableThis, PreviewPath);
            
        end
        
ObjectsAll = vertcat(ObjectsAllFields{:});
ObjectsName = InfoTable.Barcode{1};
save([FolderThisAnalysis, filesep, ObjectsName], 'ObjectsAll')

AreaNamesToKeep = ObjectsAll.Properties.VariableNames(~ismember(ObjectsAll.Properties.VariableNames, {'AreaVectors', 'NodeDegreeVector'}));
Summary = ObjectsAll(:, AreaNamesToKeep); % remove vectors
writetable(Summary, [FolderThisAnalysis, filesep, ObjectsName, '.csv'], 'WriteVariableNames', true); % Saving as comma separated file     
        
end

