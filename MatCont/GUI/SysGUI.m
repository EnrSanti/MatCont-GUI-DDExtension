classdef SysGUI
    
    methods(Static)
        
        function outhandle = new()
            
            
            global path_sys
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys filesep()];
            %salva poi (dalla gui)
            systems_standalone('init');
            %ha creato il gds vuoto
            fhandle = systems_standalone();
            if nargout == 1
               outhandle = fhandle; 
            end
        end
        
        
        function outhandle = edit(systemname)
            
            if (~ischar(systemname))
                systemname = func2str(systemname);
            end
            
            global path_sys;
            global gds;
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys filesep()];
            load( [path_sys  systemname '.mat' ]);  %overwrites gds.
            
            %-_-_-_-_-_-_%
            %if the struct loaded didn't containt the sys_type field, that
            %is added with the value of "ODE" (legacy system)
            if(~isfield(gds,'sys_type'))
                gds.sys_type = "ODE";
            end   
            %-_-_-_-_-_-_%
            
            %crea finestra gui
            fhandle = systems_standalone();
            if nargout == 1
                outhandle = fhandle;
            end
            
        end
        
        function outhandle = userfunctions(systemname)
            
             disp("uf called");
           
            
            if (~ischar(systemname))
                systemname = func2str(systemname);
            end
            global path_sys
            global gds;
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys '/'];
            systemmatfile = [path_sys  systemname '.mat' ];
            load(systemmatfile); %overwrites 'gds'.
            %caricato il sistema..
            
            fhandle = userfun_standalone();
            if nargout == 1
                outhandle = fhandle;
            end
            
        end
        
        
        function sys = gui_loader(systemcmd)
            global gds;
            gds = [];
            disp("creo gui");
            fhandle = systemcmd(); 
            if isvalid(fhandle); uiwait(fhandle); end
            %ha già scritto il file a questo punto
            if ~isempty(gds) && isfield(gds, 'ok') && gds.ok
                sysname = gds.system;
                [path_sys, ~, ~] = fileparts(which('standard.m'));
                syspath = fullfile(path_sys, [sysname '.mat']);
                sys = CLSystem(syspath);
                
            else
                sys = [];
            end
        end
        
    end

end
