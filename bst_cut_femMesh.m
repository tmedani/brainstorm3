function bst_cut_femMesh(iSubject)
%% SECTION 1 : Get the data
disp('cut fem mesh')
bst_progress('start', 'Cutting Mesh','Load mesh...')

%% Get the input data
ProtocolInfo     = bst_get('ProtocolInfo');
ProtocolSubjects = bst_get('ProtocolSubjects');
% Default subject
if (iSubject == 0)
    sSubject = ProtocolSubjects.DefaultSubject;
    % Normal subject
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end

%% Get the mesh file
% Get the conductivity values
FemFiles = file_fullpath(sSubject.Surface(sSubject.iFEM).FileName);
[filepath,name,ext]  = fileparts(FemFiles);
% Load the mesh
headFEM=  load(FemFiles);

% set a threshold point :
maxHeadPointZ = max(headFEM.Vertices(:,3));
minHeadPointZ = min(headFEM.Vertices(:,3));
delta = (maxHeadPointZ - minHeadPointZ);


ratio = 1/3; % par to remove from the mesh
% Ask user for the node shifting
[res, isCancel]  = java_dialog('input', 'Ratio of the bottom head to remove:', ...
    'Botton ratio to remove', [], num2str(ratio));
if isCancel || isempty(str2double(res))
    return
end
ratio = str2double(res);

threshold = minHeadPointZ+ratio*delta;
z0 = threshold;

bst_progress('start', 'Cutting Mesh','Cutting mesh...')

tic,
newFace = headFEM.Elements;
newTissue = headFEM.Tissue;

zCoordinate =  headFEM.Vertices(:,3);
temp = zCoordinate(newFace) <= z0;
toRemove = find(sum(temp,2));

newFace(toRemove,:) = [];
newTissue(toRemove,:)  = [];


surfIN.vertices = headFEM.Vertices;
surfIN.faces = newFace;
surfIN.tissue = newTissue;
[NewSurf,UsedV]=delete_unused_vertices(surfIN);
timeOfCuting = toc

bst_progress('start', 'Cutting Mesh','save mesh to database ...')

% Create output structure
FemMat = db_template('femmat');
FemMat = headFEM;
FemMat.Vertices = NewSurf.vertices;
FemMat.Elements =NewSurf.faces;
FemMat.Tissue =NewSurf.tissue;
FemMat.Comment = sprintf('FEM cut %dV (%d, %d layers)', length(FemMat.Vertices), ratio, length(FemMat.TissueLabels));
% Add history
FemMat = bst_history('add', FemMat, 'process_cut_fem');
% Save to database
FemFile = file_unique(bst_fullfile((filepath), sprintf('tess_fem_%s_%dV.mat', 'cut_mesh', length(FemMat.Vertices))));
bst_save(FemFile, FemMat, 'v7');
db_add_surface(iSubject, FemFile, FemMat.Comment);

            bst_progress('stop', 'Cutting mesh...');

end
