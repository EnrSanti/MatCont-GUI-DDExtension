classdef GUIPlotOutputter < handle
    properties
        lastindex = 0;
        dimension;
        
        fx
        fy
        fz
        
        axeshandle
        group 
        getData
        
        pointplot
        
        plotconfig;
        pointloader;
        curvelabel;
    end
    
    methods
        function obj = GUIPlotOutputter(axeshandle, fx, fy, fz, pointplotvector, plotconfig, pointloader, curvelabel)
            
            obj.fx = fx;
            obj.fy = fy;
            obj.fz = fz;
            obj.axeshandle = axeshandle;
            
            if isempty(fz)
                obj.dimension = 2;
                datagetter = @(data, s, ind) {fx(data{:}, s, ind, ind(end)),  fy(data{:}, s, ind, ind(end))};

            else
                obj.dimension = 3;
                datagetter = @(data, s, ind) {fx(data{:}, s, ind, ind(end)),  fy(data{:}, s, ind, ind(end)), fz(data{:}, s, ind, ind(end))};
            end  
            obj.pointplot = any(pointplotvector);
            
            if obj.pointplot 
                obj.getData = @(data, s, ind) dimensionEqualizer(datagetter(data, s, ind));
                
            else
               obj.getData = datagetter; 
            end
            
            
            obj.reset();
             
            if nargin < 7
                pointloader = [];
                curvelabel = '';
            end
            obj.pointloader = pointloader;
            obj.plotconfig = plotconfig; 
            obj.curvelabel = curvelabel;
            
            
        end

        function output(obj, contdata, s, ind)
            
            
            
            i = ind(end);
            if i > obj.lastindex
                j = max(1, obj.lastindex + obj.pointplot); %explanation needed -> ok
                
                if ~isempty(obj.plotconfig.modification)
                    for k = j:(i - 1 +obj.pointplot)
                        data = obj.getData(contdata, s, k:k+1-obj.pointplot);
                        modification = obj.plotconfig.modification(contdata, k+1-obj.pointplot);
                        line(data{:}, 'Parent', obj.group,  obj.plotconfig.curve{:}, modification{:});
                    end
                else
                    
                    data = obj.getData(contdata, s, j:i);
                    line(data{:}, 'Parent', obj.group,  obj.plotconfig.curve{:});
                end
                obj.lastindex = i;

                if ~isempty(s) && i == s.index && s.index > 1
                   obj.handleSpecial(contdata, s, s.index)
                end
            end
            
            
        end
                
        function outputPoint(~, varargin)  %not implemented, only in Simulation.
        end    
        
        function handleSpecial(obj, contdata, s, index)
            data = obj.getData(contdata, s, index);
            labeldata = cell(1, length(data));
            for k = 1:length(data)
                labeldata{k} = data{k}(1,:);
            end
            
            if obj.pointplot
                line(data{:},'Parent', obj.group, obj.plotconfig.specialcurve{:});
                %line(labeldata{:},'Parent', obj.group, 'Color', 'red', 'marker', 's', 'linestyle', 'none');
            else
                line(labeldata{:},'Parent', obj.group, obj.plotconfig.specialpoint{:},'ButtonDownFcn', @(o,e) selectPoint(o, obj.pointloader, obj.curvelabel), 'UserData', s);
            end
            
            
            skewdirection = obj.getSkewDirection(contdata, s, index, labeldata);
            obj.placeLabel(s, labeldata, skewdirection);
        end
        
        function vector = getSkewDirection(obj, contdata, s, index, labeldata)
            d = axis(obj.axeshandle);
            
            vector = [0 0 0];
            for dim = 1:obj.dimension
                i = dim - 1;
                vector(dim) = d(2 + 2*i) - d(1 + 2*i);
            end
            vector = 0.012*vector;
            
            if index > 2 && obj.dimension == 2 && ~obj.pointplot
                prevpoint = obj.getData(contdata, s, index-1);
                direction = [labeldata{1} - prevpoint{1}; labeldata{2} - prevpoint{2}];
                direction(1) = direction(1) / (d(2) - d(1));
                direction(2) = direction(2) / (d(4) - d(3));
                direction = [-direction(2); direction(1)];
                direction(1) = direction(1) * (d(2) - d(1));
                direction(2) = direction(2) * (d(4) - d(3));        
                direction =   sqrt( (0.012*(d(2) - d(1)))^2 +  (0.012*(d(4) - d(3)))^2 ) * direction/norm(direction);
                %line(obj.axeshandle, [labeldata{1} - direction(1), labeldata{1} + direction(1)],  [labeldata{2} - direction(2), labeldata{2} + direction(2)], 'color', 'green')
                vector = direction;
            end

        end
        
        function placeLabel(obj, s, labeldata, skewdirection)
               
                for dim = 1:obj.dimension
                   labeldata{dim} = labeldata{dim} + skewdirection(dim);
                end
                
                
                th = text(labeldata{:}, strip(s.label) , 'Parent' ,obj.axeshandle, 'ButtonDownFcn' , @(o,e) selectPoint(o, obj.pointloader, obj.curvelabel), 'UserData', s, obj.plotconfig.label{:}); 
        end
        
        
        function plotSolution(obj, solution)
           global session;
           obj.reset();
           
           %-_-_-_-_-_-_%
           sessionGDS=session.settings.fields.system; 
           if(sessionGDS.sys_type=="DDE")
               %recalculate first RE row(S)
               [ntst,ncol]=solution.getDiscretization();
               d1=sessionGDS.no_RE;

               if(d1>0) %se non ha re non si fa
                   %in globals.campi prendi nomi, togli cds, recupera
                   %struct corretta
                   contentGlobals=fieldnames(solution.globals);
                   structDs=contentGlobals{1};
                   if(strcmp(structDs,"cds"))
                       structDs=contentGlobals{2};
                   end
                   
                   %take the columns
                   [r,c]=size(solution.x);
                   d2=sessionGDS.dim-sessionGDS.no_RE;
                   %handle to get the functions in the m file
                   hdl=session.settings.fields.system.handle;
                   [~,rhsFun]=hdl();
                   for indexCol=1:c
                       for indexBlock=1:ntst*ncol+1
                            for indexRowRE=1:d1
                                %getting the right parameters, converting
                                %to a row vector & num2cell
                                paramsRow=(solution.globals.(structDs).P0)';
                                parametersLine=num2cell(paramsRow);      
                                activeParams=solution.globals.(structDs).ActiveParams;
                                
                                parametersLine(activeParams)=num2cell(solution.x(end-length(activeParams)+1:end));
                                %parametersLine(activeParams)=num2cell(solution.x(paramsRow(end-length(activeParams)+1:end)));
                                %negli appunti che ho segnato ieri, era
                                %come nella linea commentata, ma non mi
                                %sembra sia corretto quel paramRow, nel
                                %caso mi sia sbagliato adesso, decommentare
                                %la linea e commentare quella sopra (182)
                                solution.x((d2*(sessionGDS.no_discretizationPoints+1)+d1*sessionGDS.no_discretizationPoints)*(indexBlock-1)+d2*(sessionGDS.no_discretizationPoints+1)+indexRowRE,indexCol)=rhsFun{indexRowRE}(0,solution.x((d2*(sessionGDS.no_discretizationPoints+1)+d1*sessionGDS.no_discretizationPoints)*(indexBlock-1)+(1:(d2*(sessionGDS.no_discretizationPoints+1)+d1*sessionGDS.no_discretizationPoints)),indexCol),parametersLine{:});
                            end
                        end
                   end
               end
           end
           %-_-_-_-_-_-_%
           
           
           obj.output({solution.x, solution.h, solution.f}, solution.s(1), 1:size(solution.x, 2));
           
           for i = 2:length(solution.s)-1
              obj.handleSpecial({solution.x, solution.h, solution.f},solution.s(i), solution.s(i).index);
               
           end
        end
        
        function reset(obj)
           obj.lastindex = 0; 
           if ~obj.pointplot
               obj.group = hggroup('Parent' , obj.axeshandle );
           else
              obj.group = obj.axeshandle; 
           end
        end
        
    end
end

function data = dimensionEqualizer(data)
    dims = cellfun(@(x) size(x, 1), data);
    dim = size(data{find(dims ~= 1, 1)}, 1);
    
    for i = find(dims == 1)
        data{i} = repmat(data{i}, dim, 1);
    end
end

function selectPoint(guiobject, pointloader, label)
    if isempty(pointloader)
        disp(repmat('-', 1, 80))
        disp(guiobject.UserData);
        if isfield(guiobject.UserData, 'data')
           disp(guiobject.UserData.data); 
        end
        disp(repmat('-', 1, 80))
    else
        if isvalid(pointloader)
            pointloader.loadPoint(label, guiobject.UserData);
        else
            fprintf(2, 'Connection with current session has been lost, unable to load in point\n');
            
        end
    end

end


