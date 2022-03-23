classdef CLSystem
    
    properties
        %-_-_-_-_-_-_%
        sys_type = '';
        no_discretizationPoints = 0;
        no_delays = 0;
        %-_-_-_-_-_-_%
        name = ''
        coordinates = {}
        parameters = {};
        dim = 0;
        time = 't';
        handle = [];
        userfunctions = []
        ufdata = [];
        diagramlocation = '';
        derstr
        equations = {};
    end
    
    
    methods
        %-_-_-_-_-_-_%
        function sys_type = getSys_type(obj)
           sys_type = obj.sys_type; 
        end
        function Mpoints = getNo_discretizationPoints(obj)
           Mpoints = obj.no_discretizationPoints; 
        end
        function nDelays = getNo_delays(obj)
           nDelays = obj.no_delays; 
        end
        %-_-_-_-_-_-_%
        function name = getName(obj)
           name = obj.name; 
        end
        function r = getParameters(obj)
            r = obj.parameters;
        end
        function r = getCoordinates(obj)
           r = obj.coordinates; 
        end
        
        function obj = CLSystem(datfile)
            if nargin == 0
               return; 
            end
                
                
            location = which(datfile);
            assert(~isempty(location), ['system does not exist: ', datfile]);
            location(end-3:end) = [];
            datfile = load(datfile); % carica il .mat
            gds = datfile.gds;
            
            %-_-_-_-_-_-_%
            %if the loaded struct contains the field sys_type then
            %initialize the field of the returned object with that value,
            %otherwise, catch exception, and initialize the sys type to
            %"ODE"
            try
                obj.sys_type = gds.sys_type;
                % se è ode non carico neanche
                if(strcmp(obj.sys_type,"DDE"))
                    obj.no_discretizationPoints=gds.no_discretizationPoints;
                    obj.no_delays=gds.no_delays;
                    disp("DDE "+obj.no_discretizationPoints);
                end
            catch 
                obj.sys_type = systemType.getType(1);  %"ODE"
                disp("legacy system");
            end 
            %-_-_-_-_-_-_%
            
            
            obj.name = gds.system;
            obj.coordinates = datfile.gds.coordinates(:,1)';
            obj.parameters = datfile.gds.parameters(:, 1)';
            obj.dim = gds.dim;
            obj.time = gds.time{1};
            obj.equations = gds.equations;
            if isfield(gds, 'userfunction')
                obj.userfunctions = gds.userfunction;
                obj.ufdata = gds.options.UserfunctionsInfo;
                for i = 1:length(obj.userfunctions)
                   obj.ufdata(i).valid = ~isempty(obj.userfunctions{i}); 
                   %if ~obj.ufdata(i).valid
                   obj.ufdata(i).state = 0;
                   obj.ufdata(i).label = strip(obj.ufdata(i).label);
                   %end
                end
            else
                obj.userfunctions = {};
                obj.ufdata = [];
            end
            
            assert (length(obj.userfunctions) == length(obj.ufdata), 'userfunction mismatch')

            obj.diagramlocation = location;
            obj.handle = str2func(obj.name);
            %compute der-string: ex: SSSNN, SNNNN, SSSSS
            derstr = repmat(' ', 1, 5);
            for i=1:5
                for j=1:4
                    if (gds.der(j,i)==1)
                        switch j
                            case 1
                                derstr(i)='N';
                            case 3
                                derstr(i)='R';
                            case 4
                                if (exist('sym')==2)
                                    derstr(i)='S';
                                else
                                    derstr(i)='F';
                                end
                        end
                    end
                end
            end
            obj.derstr = derstr;
            obj.diagramInit();
            %qui ho obj :)
        end
        
        function s = getUFString(obj)
            s = '(none)';
            
            if ~isempty(obj.ufdata)
               valids = [obj.ufdata.valid];
               labels = {obj.ufdata(valids).label};
                
                if ~isempty(labels)
                   s = strjoin(labels, ', '); 
                end
            end
            
            
        end
        
        function diagramInit(obj)
            if isempty(obj.diagramlocation); return; end
             %create diagram directory if not exists
            if (exist( obj.diagramlocation, 'dir') ~= 7)
                mkdir(obj.diagramlocation);
            end
            diagrams = CLDiagram.getDiagramList(obj.diagramlocation);
            if (isempty(diagrams))
                mkdir(fullfile(obj.diagramlocation, 'diagram'));
            end           
        end
        
        function s = getDerInfo(obj)
           s = obj.derstr; 
        end
        
        function t = getTimeName(obj)
           t = obj.time; 
        end
        function eq = getEquations(obj)
           eq = obj.equations; 
        end
        function d = getUFData(obj)
            d = obj.ufdata;
        end
        function labels = getUFLabels(obj)
            if isempty(obj.ufdata)
                labels = {};
            else
                labels = {obj.ufdata.label};
            end
        end
        
        %Settings Interface:
        
        function value = getValue(obj)
            if isempty(obj.name)
                value = [];
            else
                value = obj;
            end
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            valid = 0;
            msg = 'not allowed';
        end
        
        function id = getGroupID(obj)
            id = 0;
        end
        function id = getSubGroupID(obj)
            id = 0;
        end
        function id = getItemID(obj)
            id = 1;
        end
        function dim = getDim(obj)
           dim = obj.dim; 
        end
        function s = toString(obj)
            s = obj.name;
        end
        
        function newobj = copy(obj, ~)
            newobj = obj;
        end
        function h = getHelpStr(obj)
            h = 'no help';
        end
        function b = isVisible(obj)
            b = 1;
        end
        function setVisible(obj, bool)
        end
       function row = getIDs(obj)
          row = [obj.getGroupID(), obj.getSubGroupID(), obj.getItemID()];
       end
       function box = renderGUI(varargin)
          box = []; 
       end
       function p = getDiagramPath(obj)
           p = obj.diagramlocation;
           
       end
       function t = getValueType(~)
            t = 'NONE';
       end
       function b = sanityCheck(obj, settings)
           %fprintf(2, 'remember, after every loading of a state, system chould be updated\n');
           b = 1;
       end
    end
    
end
