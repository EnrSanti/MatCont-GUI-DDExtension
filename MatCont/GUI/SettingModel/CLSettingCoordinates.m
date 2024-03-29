classdef CLSettingCoordinates < CLSetting

    properties
       coordinates 
    end
    
    methods
        function obj = CLSettingCoordinates(settings, coordinates)
            
            global session;
            dim = length(coordinates);
            %-_-_-_-_-_-_%
            %try if the system is dde one, if so, then set the dimension to
            %the actual dimension of the discretizing system
            try
               if(session.settings.fields.system.sys_type=="DDE")
                    no_RE=session.settings.fields.system.no_RE;
                    m=session.settings.fields.system.no_discretizationPoints;
                   
                    dim = (((dim-no_RE)*(m+1))+m*no_RE);  %modifica d2*(m+1)+d1*m
               end
            catch exception
                dim=dim;
            end
            %-_-_-_-_-_-_%
            obj = obj@CLSetting('coord', zeros(1, dim), InputRestrictions.vector(dim), 2, 1, length(coordinates) + 1, '~~~');
            obj.coordinates = coordinates;      
            
            if ~isempty(settings)
                obj.installGhosts(settings);
            end
        end
        function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
            %GUIrender = obj.validitycheck.getGUIRenderer();
            %box = LayoutBox({GUIrender(panelhandle, settings, label, options{:},'BackgroundColor', [0.90 0.90 0.90] )});
            box = [];
            
        end
        
        function installGhosts(obj, settings)
            for i = 1:length(obj.coordinates)
                coordinate = obj.coordinates{i};
                settings.addSetting(['co_', coordinate], CLSettingCoordinate(obj, i));
            end
            
            
        end
        
       function newobj = copy(obj, newsettings)
           newobj = CLSettingCoordinates(newsettings, obj.coordinates);
           newobj.value = obj.value;
           obj.copyOver(newobj);
       end        
        
       function b = sanityCheck(obj, settings)
            b = obj.internalSanityCheck(settings);
           
            if ~b
                for i = 1:length(obj.coordinates)
                    coordinate = obj.coordinates{i};
                    settings.removeSetting(['co_', coordinate]);
                end
            end
            
       end
       
       function b = internalSanityCheck(obj, settings)
          system = settings.system;
          if isempty(system); b = 0; return; end
          coordinates = system.getCoordinates();
          if length(coordinates) ~= length(obj.coordinates); b = 0; return; end
          
          for k = 1:length(obj.coordinates)
                if ~strcmp(obj.coordinates{k}, coordinates{k})
                    b = 0;
                    return
                end
          end
          b = 1;
       end          
    end
end
