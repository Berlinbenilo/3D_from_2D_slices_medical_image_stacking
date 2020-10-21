% No of input images
n = 41; 
X = zeros(512,512,n);
folder = 'C:\Users\ADHI\Downloads\surface2volume';
filePattern = fullfile(folder, '*.jpg');
srcFiles = dir(filePattern);

for i = 2:n
    fileName = ['C:\Users\ADHI\Downloads\surface2volume\', srcFiles(i).name];
    I=imread(fileName);
    I = rgb2gray(I);
    I = imresize(I,[512,512]);
    I = imclearborder(I);
    X(:,:,i) = I;
end


D = uint8(X);
D = squeeze(D);
D = padarray(D,[5 5 5],'both');
      
% Create an isosurface
Ds = smooth3(D);
surface = isosurface(Ds,5);
      
% Display the surface
figure;
hiso = patch('Vertices',surface.vertices,...
                   'Faces',surface.faces,...
                   'FaceColor',[1,.75,.65],...
                   'EdgeColor','none');

view(45,30) 
axis tight 
daspect([1,1,.4])
lightangle(45,30); 
set(gcf,'Renderer','zbuffer'); lighting phong
isonormals(Ds,hiso)
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)

% Reconstruct the volume and display it as montage
OV = surface2volume(surface,[],1);
nDims = size(OV);