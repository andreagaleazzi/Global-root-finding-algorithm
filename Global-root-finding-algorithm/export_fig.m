function export_fig(fig, folder, filename)

% If folder does not exist create it
folder_PNG = strcat(folder,'/PNG/');
folder_PDF = strcat(folder,'/PDF/');

if ~exist(folder, 'dir')
    mkdir(folder);
end
if ~exist(folder_PNG, 'dir')
    mkdir(folder_PNG);
end
if ~exist(folder_PDF, 'dir')
    mkdir(folder_PDF);
end

path_EMF = strcat(folder,'/',filename,'.emf');
path_PNG = strcat(folder_PNG,'/',filename,'.png');
path_PDF = strcat(folder_PDF,'/',filename,'.pdf');

if exist(path_EMF, 'file')
    delete (path_EMF);
end
if exist(path_PNG, 'file')
    delete (path_PNG);
end
if exist(path_PDF, 'file')
    delete (path_PDF);
end

exportgraphics(fig, path_EMF)
exportgraphics(fig, path_PDF)
exportgraphics(fig, path_PNG, 'Resolution', 200)

end

