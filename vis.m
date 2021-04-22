%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% visualization of the instabile mode %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
addpath('freefem_matlab_octave_plot-public/release-v2.0/demos/ffmatlib')
addpath('Data')
addpath('Eigenvalues')
addpath('Eigenvectors')
addpath('Flows')

%loading data from ff++
%load mesh
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('Cyl.msh');
%load connectivity
VPh=ffreaddata('cyl_VPh.txt');
%load baseflow
[Ub] = ffreaddata('Vectored_baseflow.txt');
[VhP2,ub]=ffvectorget({'P2','P2','P1'}, VPh, Ub, 1);
[VhP2,vb]=ffvectorget({'P2','P2','P1'}, VPh, Ub, 2);
[VhP1,pb]=ffvectorget({'P2','P2','P1'}, VPh, Ub, 3);
%load eigevalues
id = fopen('eig.txt');
d = textscan(id,'%f %f');
e = cell2mat(d);
fclose('all');

%number of eigenvalues calculated
nEV = length(e(:,1));

%load eigenvectors and build eigenvalues in cartesian form
la = zeros(1,nEV);
up = zeros(length(ub),nEV);
vp = zeros(length(vb),nEV);
pp = zeros(length(pb),nEV);
for i =1:nEV
    la(i) = e(i,1)+1i*e(i,2);
    name = sprintf('vectorized_eigenvectors_%i.txt',i);
    %     namep=sprintf('p_eigenvector_%i.txt',i);
    %     nameu=sprintf('u_eigenvector_%i.txt',i);
    %     namev=sprintf('v_eigenvector_%i.txt',i);
    [Up] = ffreaddata(name);
    [~, up(:,i)]=ffvectorget({'P2','P2','P1'}, VPh, Up, 1);
    [~, vp(:,i)]=ffvectorget({'P2','P2','P1'}, VPh, Up, 2);
    [~, pp(:,i)]=ffvectorget({'P2','P2','P1'}, VPh, Up, 3);
end

%selecting the eigenvalue withe the greatest real part and the relative
%modes

[~,i] = max(real(la));

lai = la(i);
upi = up(:,i);
uvi = vp(:,i);
ppi = pp(:,i);

%linearized flow
A = 25; %a sufficently big amplitude to visualize the perturbation

% temp = 0:0.1:5;
% figure
% for i = 1: length(temp)
%     u = ub+A*real(upi*exp(lai*temp(i)));
%     ffpdeplot(p,b,t, ...
%         'VhSeq',VhP2, ...
%         'XYData',u, ...
%         'Mesh','off', ...
%         'CBTitle','u', ...
%         'XLim',[-5 20],'YLim',[-5 5], ...
%         'Title','Instabile mode',...
%         'MColor','b','Boundary','off','ColorMap',jet);
%     drawnow
%     hold on
% end

phi=linspace(0,2*pi,50);
phase=exp(1i*phi);
clear M
for i=1:length(phase)
    utemp=real(up*phase(i));
    vtemp=real(vp*phase(i));
    ffpdeplot(p,b,t,'VhSeq',VhP2,'XYData',utemp,'Mesh','off','MColor','b','Boundary','off','ColorMap',jet);
    drawnow
end
