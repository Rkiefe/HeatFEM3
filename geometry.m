classdef geometry
	properties
		Hmax
		Hmin
		MeshGradation
	end

	methods
		% Makes the mesh of the model according to a dynamic set of settings
        function obj = mesh(obj,model,settings)
			if nargin < 2
				generateMesh(model);
				return
			end

			command = '';
			for key = keys(settings)
			    % Convert the value to string explicitly
			    value = settings(key{1});
				
				if isequal('refinedOrder',key{1})
					continue
				end

				if isequal('Hedge',key{1})
					command = [command, '"', key{1}, '",{['];

					for i = 1:length(value)-1
						command = [command, num2str(value(i), '%d') ,','];
					end
					command = [command,num2str(value(end), '%d'),'],',num2str(settings('refinedOrder')),'},'];

				elseif isnumeric(value)
			        value =  num2str(value, '%d');
			        command = [command, '"', key{1}, '",', value, ','];
			    
			    else
			   		command = [command, '"', key{1}, '","', value, '",'];
				end
			end

			% Remove the trailing comma
			command = command(1:end-1);

			% Construct the final command variable
			commandCell = eval(['{', command, '}']);

			generateMesh(model,"GeometricOrder","linear", commandCell{:});

			obj.Hmax = model.Mesh.MaxElementSize;
			obj.Hmin = model.Mesh.MinElementSize;
			obj.MeshGradation = model.Mesh.MeshGradation;
						
        end % Make the mesh
	end

	methods(Static)

		% Adds cuboids to a model. If no geometry is defined, the geo added is the container
        function message = addCuboid(model,cuboid)
        	message = '';

        	% Cuboid must be a struct
        	if ~isequal(class(cuboid),'struct')
        		disp("Cuboid must be a struct with dimensions and position parameters")
        		return
        	end

        	if ~isempty(find(cuboid.dimensions==0,1))
        		message = 'Invalid dimensions for the geometry. No singularities please.';
    			return
    		end

    		g1 = multicuboid(cuboid.dimensions(1),cuboid.dimensions(2),cuboid.dimensions(3));
			
			% First scale the object if necessary
			if isfield(cuboid, 'scale') 
				g1 = scale(g1,cuboid.scale);
			end

    		% Then rotate if necessary
    		if isfield(cuboid, 'angle') && isfield(cuboid, 'axis') 
    			rotate(g1,cuboid.angle,[0,0,0],cuboid.axis);
			end

			% Move to desired position
			if ~isfield(cuboid,'centered') 
				g1 = translate(g1,cuboid.position);
			elseif ~cuboid.centered
				g1 = translate(g1,cuboid.position);
			elseif cuboid.centered
				center = mean(g1.Vertices);
				g1 = translate(g1,-center);
			end

			% addcell(g1)
			if isempty(model.Geometry)
				message = 'Model was empty. Container now defined. Now add the refrigerant.';
				model.Geometry = g1;
				container = g1;
			else
				message = 'Model already had a container. Adding the refrigerant.';
				try
					addCell(model.Geometry,g1);
				catch exception
					message = 'Unable to add cell. ';
					message = [message,exception.message];
				end
			end

        end % Adds cuboids \ makes the container

        function message = import(model,file,settings)
        	
        	message = '';

        	if ~isfield(settings,'position') && ~isfield(settings,'centered')
        		disp("Missing the desired position")
        		return
        	end


        	% Import model
        	g2 = importGeometry(file);

        	% Scale if necessary
        	if isfield(settings,'scale')
        		g2 = scale(g2,settings.scale);
        	end

        	% Then rotate if necessary
    		if isfield(settings, 'angle') && isfield(settings, 'axis') 
    			rotate(g2,settings.angle,[0,0,0],settings.axis);
			end

        	% Move to desired position
        	if ~isfield(settings,'centered') 
				g2 = translate(g2,settings.position);
			elseif ~settings.centered
				g2 = translate(g2,settings.position);
			elseif settings.centered
				if isempty(g2.Vertices)
					message = [message,' Cannot center spheres'];
					return
				end

				center = mean(g2.Vertices)
				g2 = translate(g2,-center);
			end

			% Add a container
			if isempty(model.Geometry)
			
				% Make a container for the stl
				X = 0.5*10*(max(g2.Vertices(:,1))-min(g2.Vertices(:,1)));
				Y = 0.5*10*(max(g2.Vertices(:,2))-min(g2.Vertices(:,2)));
				Z = 0.5*10*(max(g2.Vertices(:,3))-min(g2.Vertices(:,3)));

				gContainer = multicuboid(X,Y,Z);
				Container_center = mean(gContainer.Vertices);
				gContainer = translate(gContainer,-Container_center);

				% Add the container
				model.Geometry = gContainer;

				% Add the stl
				try
					addCell(model.Geometry,g2);
					message = [message,'Created a container and imported the model'];
				catch exception
					message = 'Unable to add cell. ';
					message = [message,exception.message];
				end
				
				return
			end

			try
				if isempty(model.Geometry)
					model.Geometry = g2;
				else
					addCell(model.Geometry,g2);
				end

				message = 'File imported with success.';
			catch exception
				message = 'Unable to add cell. ';
				message = [message,exception.message];
			end
			
        end

        function plotMesh(model,settings)
        	figure

        	if nargin < 2
				pdeplot3D(model.Mesh);
				return
			end

			command = '';
			for key = keys(settings)
			    % Convert the value to string explicitly
			    value = settings(key{1});
			    if isnumeric(value)
			        value =  num2str(value, '%d');
			        command = [command, '"', key{1}, '",', value, ','];
			    else
			   		command = [command, '"', key{1}, '","', value, '",'];
				end
			end

			% Remove the trailing comma
			command = command(1:end-1);

			% Construct the final command variable
			commandCell = eval(['{', command, '}']);

			
			pdeplot3D(model.Mesh,commandCell{:}); fig = gca;
        end

        function plotModel(model,settings)
        	figure
        	if nargin < 2
				pdegplot(model)
				return
			end

			command = '';
			for key = keys(settings)
			    % Convert the value to string explicitly
			    value = settings(key{1});
			    if isnumeric(value)
			        value =  num2str(value, '%d');
			        command = [command, '"', key{1}, '",', value, ','];
			    else
			   		command = [command, '"', key{1}, '","', value, '",'];
				end
			end

			% Remove the trailing comma
			command = command(1:end-1);

			% Construct the final command variable
			commandCell = eval(['{', command, '}']);

			pdegplot(model,commandCell{:})
        end
    end
end
