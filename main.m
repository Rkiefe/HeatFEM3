%{
	Solves the heat equation with linear thermal conductivity
	
	The boundary condition: k grad T . n = h(Text-T)
	h is the heat transfer coefficient (reciprocal to the thermal insulance)
%}

close all
clear
clc

% ------------ User defined variables and constants ------------
L = [10,10,10]; % Object dimensions

Time = 5e-3;	% time of experiment, s

rho_in = 7.9;			% density, g/cm^3
rho_out = 1.225e-3; 	% g/cm^3 wikipedia

k_in = 8.8; 			% conductivity, W/m/K (https://www.matweb.com/search/datasheet_print.aspx?matguid=750a3dd8d69b44b79468fbaf72a2beef)
k_out = 2.262e-2;		% wikipedia

Cp_in = 0.3; 			% J/g/K
Cp_out = 1.012;			% J/g/K wikipedia

Ti_f = @(x,y,z) (293);	% Initial temperature function
Text = 280;				% Outside temperature, K
hCoef = 1e9; 			% Heat transfer coefficient, W/m^2/K

heatSource = 0;	% volumetric heat source (W/cm^3)

% Solver settings
dt = 1e-3; 	% time step, s

% Mesh parameters
Hmax = 0; % 0 -> let matlab choose
Hmin = 0; % ...

% Model view settings
viewSettings = containers.Map;
viewSettings('FaceAlpha') = 0;
viewSettings('FaceLabels') = 'on';
% viewSettings('EdgeLabels') = 'off';

% ---------------------------------------------

% >> Make an empty model
model = createpde;

% Handle that model through my own class - geometry
sketch = geometry();

% >> Add Container
container.dimensions = 5*L; % cm
container.centered = true;
sketch.addCuboid(model,container);

% >> Add Refrigerant
cube.centered = true;
cube.dimensions = L; % cm
sketch.addCuboid(model,cube)

% >> Import stl files
stl_folder = "..\..\..\..\FEMCE\STL_files\";
stl_file = "simple_SnowFlake.stl";
ref.scale = 0.001;
ref.centered = true;
% ref.axis = [1,0,0];
% ref.angle = 90;
% sketch.import(model,stl_folder+stl_file,ref)

% >> Plot model
sketch.plotModel(model,viewSettings); fig = gcf;
pause()
close(fig);

% >> Mesh models
meshSettings = containers.Map; % Matlab's dictionary equivalent
meshSettings('Hmax') = Hmax; % 10/100
meshSettings('Hmin') = Hmin; % 0.2/100
% meshSettings('Hgrad') = 2; % 0.1/100
% meshSettings('Hedge') = 13:(12*8);
% meshSettings('refinedOrder') = 0.1/100;

disp("Generating mesh...")
sketch = sketch.mesh(model,meshSettings); clear meshSettings

% Process mesh
disp("Processing mesh...")
mesh = manageMesh();
mesh = mesh.getNumberOfElementsInside(model);
mesh = mesh.processMesh(model);

disp("Number of elements: "+mesh.nt)
disp("Number of elements inside: "+mesh.nInside)
pause()

% >> Boundary conditions
g = zeros(mesh.numFaces,1);
g(6+[1:6]) = hCoef; % [3,4,5,8]

% >> Initial temperature
T = zeros(mesh.nv,1);
for i = 1:mesh.nv
	r = mesh.p(1:3,i);

	T(i) = Ti_f(r(1),r(2),r(3));
end
Ti = T;

% >> Stiffness matrix
kTh = zeros(mesh.nt,1) + k_out;
kTh(mesh.InsideElements) = k_in;

A = stiffnessMatrix(mesh,kTh);

% >> Mass matrix
c = zeros(mesh.nt,1) + rho_out*Cp_out;
c(mesh.InsideElements) = rho_in*Cp_in;
M = massMatrix(mesh,c);

% Load Vector | skiped if heatSource == 0
if heatSource ~= 0
	q_V = zeros(mesh.nt,1); q_V(mesh.InsideElements) = heatSource;
	b = loadVector(mesh,q_V);
else
	b = zeros(mesh.nv,1);
end

% Boundary matrix
R = boundaryMatrix(mesh,g);

% Boundary vector
r = boundaryVector(mesh,Text*g);

% >> Run
t = 0;
while t < Time
	% Solve matrix equation
	T_new = (M+dt*(A+R))\(dt*(b+r) + M*T);

	% Update
	t = t + dt;
	T = T_new;

	disp("At "+ 100*t/Time +" %")
end



InsideNodes = zeros(mesh.nv,1);
for i = 1:mesh.nInside
	k = mesh.InsideElements(i);
	InsideNodes(mesh.t(:,k)) = 1;;
end

% >> Plots
scatter3(mesh.p(1,InsideNodes>0),mesh.p(2,InsideNodes>0),mesh.p(3,InsideNodes>0),40,Ti(InsideNodes>0),'filled')
cbar = colorbar;
cbar.Label.String = "T (K)";
axis equal
title("Initial Temperature")

figure
scatter3(mesh.p(1,InsideNodes>0),mesh.p(2,InsideNodes>0),mesh.p(3,InsideNodes>0),40,T(InsideNodes>0),'filled')
cbar = colorbar;
cbar.Label.String = "T (K)";
axis equal
title("Final Temperature")

function plotFaces(mesh,options)
	% Plot the mesh triangles with colors

	arguments
		mesh;
		options.fig = figure;
		options.container = 1:6;
		options.FaceAlpha = 0;
	end
	hold on

    % Add the refrigerant faces in red
    container = 1:6;
    for it = 1:length(mesh.surfaceT)
        tr = mesh.surfaceT(1:3,it);
        if ~ismember(mesh.surfaceT(4,it),container)
            patch = fill3(mesh.p(1,tr),mesh.p(2,tr),mesh.p(3,tr),'r');
        else
        	patch = fill3(mesh.p(1,tr),mesh.p(2,tr),mesh.p(3,tr),'b');
        end
    	patch.FaceAlpha = options.FaceAlpha;
    	patch.EdgeAlpha = 0;
    end

    view(3)
end

function data = loadData(folder)
	data = struct;
	files = dir(folder);
	for i = 1:numel(files)
		[~,name,ext] = fileparts(files(i).name);
		if ~isequal(ext,'.') && ~isempty(ext)
			data.(name) = importdata(folder+files(i).name);
			% data.(name) = readmatrix(folder+files(i).name);

			if isstruct(data.(name))
				data.(name) = data.(name).data;
			end
		end
	end
end

function M = massMatrix(mesh,c)
	M_k = 1/20 * [2,1,1,1;...
				  1,2,1,1;...
				  1,1,2,1;...
				  1,1,1,2]; % 3D mass matrix

	M = zeros(mesh.nv);
	for k = 1:mesh.nt
		nds = mesh.t(:,k);

		M(nds,nds) = M(nds,nds) + mesh.VE(k)*M_k*c(k);
	end
end

function A = stiffnessMatrix(mesh,mu)
	A = zeros(mesh.nv);
	for k = 1:mesh.nt
		
		% Nodes of the element
		nds = mesh.t(:,k);

		% For each node
		for i = 1:length(nds)
			[~,bi,ci,di] = abcd(mesh.p,nds,nds(i)); % basis function parameters

			for j = i:length(nds) % matrix is symetric
				[~,bj,cj,dj] = abcd(mesh.p,nds,nds(j));

				A(nds(i),nds(j)) = A(nds(i),nds(j)) + mu(k)*mesh.VE(k)*(bi*bj + ci*cj + di*dj);
				A(nds(j),nds(i)) = A(nds(i),nds(j)); % A is symmetric
			end
		end
	end
end

function b = loadVector(mesh,F)
	b = zeros(mesh.nv,1);
	for s = 1:length(mesh.surfaceT)
		nds = mesh.surfaceT(1:3,s);
		
		k = mesh.surface2element(s);
		k_nds = mesh.t(:,k);

		center = mean(mesh.p(1:3,nds),2);
		for i = 1:length(nds)
			nd = nds(i);
			[ai,bi,ci,di] = abcd(mesh.p,k_nds,nd);

			b(nd) = b(nd) + mesh.VE(k)*F(k)*(ai + bi*center(1) + ci*center(2) + di*center(3));
		end
	end
end

function R = boundaryMatrix(mesh,h)
	R_k = 1/12 *[2,1,1;...
				 1,2,1;...
				 1,1,2]; % 2D mass matrix

	R = zeros(mesh.nv);
	for s = 1:length(mesh.surfaceT)
		nds = mesh.surfaceT(1:3,s);
		
		% coordinates of the nodes of the surface triangle
		r = mesh.p(1:3,nds); 
		
		% Area of the triangle
		areaT = areaTriangle(r(1,:),r(2,:),r(3,:));
		
		R(nds,nds) = R(nds,nds) + areaT*R_k*h(mesh.surfaceT(end,s));
	end
end

function b = boundaryVector(mesh,F)
	b = zeros(mesh.nv,1);
	for s = 1:length(mesh.surfaceT)
		nds = mesh.surfaceT(1:3,s);

		k = mesh.surface2element(s);
		k_nds = mesh.t(:,k);

		% coordinates of the nodes of the surface triangle
		r = mesh.p(1:3,nds); 
		centroid = mean(r,2);

		% Area of the triangle
		areaT = areaTriangle(r(1,:),r(2,:),r(3,:));

		for i = 1:length(nds)
			[ai,bi,ci,di] = abcd(mesh.p,k_nds,nds(i));
			phi_i = ai + bi*centroid(1) + ci*centroid(2) + di*centroid(3);

			b(nds(i)) = b(nds(i)) + F(mesh.surfaceT(end,s))*phi_i*areaT;
		end
	end
end

% Area of the 3D triangle
function A = areaTriangle(xt,yt,zt)
    ons = [1, 1, 1];
    A = 0.5*sqrt(det([xt;yt;ons])^2 + det([yt;zt;ons])^2 + det([zt;xt;ons])^2);
end % Area of the 3D triangle

% Basis function coef.
function [a,b,c,d] = abcd(p,nodes,nd)

    nodes(nodes==nd) = [];

    x = [p(1,nd),p(1,nodes)]';
    y = [p(2,nd),p(2,nodes)]';
    z = [p(3,nd),p(3,nodes)]';

    psi = [1;0;0;0];

    M = [ones(4,1),x,y,z];

    aux = M\psi;

    a = aux(1);
    b = aux(2);
    c = aux(3);
    d = aux(4);

end % Basis function coef.
