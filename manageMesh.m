classdef manageMesh
	properties
		nInside
		numFaces
		VE
		p
		t
		nv
		nt
		InsideElements
		surfaceT
		surface2element
	end

	methods

		% Get the number of elements inside a model mesh
		function obj = getNumberOfElementsInside(obj,model)

			% Number of elements inside the geometry
			obj.nInside = 0;
			for f = 2:model.Geometry.NumCells
			    obj.nInside = obj.nInside + numel(findElements(model.Mesh,"region","Cell",f));
			end

		end

		function obj = processMesh(obj,model)
			% Volume of each element
			[~,obj.VE] = volume(model.Mesh); obj.VE = obj.VE';

			% Index of elements that are not the container  
			obj.InsideElements = [];
			for f = 2:model.Geometry.NumCells
			    obj.InsideElements = [obj.InsideElements,findElements(model.Mesh,"region","Cell",f)];
			end

			% Nodes and elements
			obj.p = model.Mesh.Nodes;
			obj.t = model.Mesh.Elements;

			nv = length(obj.p); % number of nodes
			nt = length(obj.t); % number of elements

			obj.nv = nv;
			obj.nt = nt;

			% Number of faces in the geometry
			obj.numFaces = model.Geometry.NumFaces;

			% add a row for node border index
			obj.p = [obj.p;zeros(1,nv)];

			% Surface triangles and their respective elements
			obj.surfaceT = [];
			obj.surface2element = [];
			
			% Create the surface triangles
			for f = 1:obj.numFaces
			    % Find nodes attached to a surface
			    N_ID = findNodes(model.Mesh,"region","Face",f);

			    % Find the elements attached to this node
			    EF = findElements(model.Mesh,"attached",N_ID);

			    % Find what elements have 3 nodes in the current target surface
			    for i = 1:length(EF)
			        k = EF(i);
			        nds = obj.t(:,k);

			        nds_s = intersect(N_ID,nds);

			        if numel(nds_s) == 3
			            obj.surfaceT = [obj.surfaceT,[nds_s;f]];
			            obj.surface2element(size(obj.surfaceT,2)) = k;
			        end
			    end
			end

			[~,idx]=unique(sort(obj.surfaceT(1:3,:)',2),'rows','stable');
			obj.surfaceT = obj.surfaceT(:,idx);
			obj.surface2element = obj.surface2element(idx);
		end

	end

	methods(Access = private, Static)
		
		% Additional functions
	    function counts = bordCount(borderArray,numFaces)
		    bins = 1:numFaces+1;
		    occur = [];
		    for ind = 1:numel(borderArray)
		        occur = [occur,borderArray{ind}];
		    end

		    % Count how many nodes belong to each surface
		    counts = histcounts(occur,bins);

	    end

	    function [ndBord,localIndex] = belong2Surface(borderArray,elNodes)
		    ndBord = [];
		    localIndex = [];
		    for ind = 1:numel(borderArray)
		        if ~isempty(borderArray{ind})
		            ndBord = [ndBord,elNodes(ind)];
		            localIndex = [localIndex,ind];
		        end
		    end
	    end
	end
end