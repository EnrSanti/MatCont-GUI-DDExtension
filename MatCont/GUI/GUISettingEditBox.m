classdef GUISettingEditBox < handle
    
    properties
        handle
        eventlistener = [];
        settings
        settingname
        OKColor;
    end
    
    methods
        
         function obj = GUISettingEditBox(parent, settings, settingname, varargin )
             

             obj.settingname = settingname;
             obj.settings = settings;
            

             
            obj.handle = uicontrol(parent , 'Style' , 'edit'  , 'Unit' , 'Pixels'  , 'String' , '' , 'Callback' , @(src,ev) obj.newValue()  , varargin{:} );  
            obj.OKColor = get(obj.handle, 'BackgroundColor');
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.setBackgroundOK();
            obj.eventlistener = settings.addlistener('settingChanged' , @(srv,ev) obj.settingChanged()); 
            obj.settingChanged();
           
        end
        
        function newValue(obj)
            string =  get(obj.handle ,'String');

            try
                if ~isempty(string)
                    x = evalin('base' , string);
                else
                    x = []; 
                end
                setting = obj.settings.getSetting(obj.settingname);
                
                [valid, errormsg] = setting.setValue(x); % è QUI!.... .-.
%               
                %-_-_-_-_-_-_%
                %extracting all the data
                global session;
                m=session.settings.fields.system.no_discretizationPoints;
                dim=session.settings.fields.system.dim;
                no_RE=session.settings.fields.system.no_RE;

                RE_coords=[""];
                actualIndex=[0];
                %if the value modified is related to a coordinate
                if(isequal(class(setting),'CLSettingCoordinate'))
                    %for each RE, update index & save the old one
                    for indiceRE=1:no_RE
                        %getting the name of the RE coordinate
                        %(co_nameCoordiante)
                        RE_coords(indiceRE)="co_"+setting.coordmodel.coordinates{dim-no_RE+indiceRE};
                        actualIndex(indiceRE)=obj.settings.fields.(RE_coords(indiceRE)).index;
                        obj.settings.fields.(RE_coords(indiceRE)).index=(dim-no_RE)*(m+1)+actualIndex(indiceRE)-(dim-no_RE);
                    end
                end
                %-_-_-_-_-_-_%
                 if (~valid)
                    obj.performErrorDisplay();
                    fprintf(2, sprintf('[%s] ERROR(%s): %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), obj.settingname, errormsg, string));
                    obj.settingChanged();
                else
                    obj.settings.refresh();
                end
                %-_-_-_-_-_-_%
                %if a coordinate has been modified, shift back all the RE
                %coords indexes
                if(isequal(class(setting),'CLSettingCoordinate'))
                    for indiceRE=1:no_RE
                        obj.settings.fields.(RE_coords(indiceRE)).index=actualIndex(indiceRE);
                    end
                end
                %-_-_-_-_-_-_%
            catch error
          
                fprintf(2, sprintf('[%s] ERROR(%s): %s: %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), obj.settingname, error.message, string));
                obj.performErrorDisplay();
                obj.settingChanged(); %restore original value
            end
       
  
        end
        function performErrorDisplay(obj)
            obj.setBackgroundERROR();
            pause(0.3)
            obj.setBackgroundOK();
            
        end
        function setBackgroundOK(obj)
           set(obj.handle, 'BackgroundColor' , [1 1 1]);
           set(obj.handle, 'BackgroundColor' , obj.OKColor);
        end

        function setBackgroundERROR(obj)
           set(obj.handle, 'BackgroundColor' ,  [1    0.3    0.3]);
        end
                       
        
        function settingChanged(obj)
          %-_-_-_-_-_-_%
            %extracting all the data
            global session;
            setting=session.settings;
            if(session.settings.fields.system.sys_type=="DDE")
                m=session.settings.fields.system.no_discretizationPoints;
                dim=session.settings.fields.system.dim;
                no_RE=session.settings.fields.system.no_RE;

                RE_coords=[""];
                actualIndex=[0];
                %if the value modified is related to a coordinate
                
                    %for each RE, update index & save the old one
                    for indiceRE=1:no_RE
                        %getting the name of the RE coordinate
                        %(co_nameCoordiante)
                        RE_coords(indiceRE)="co_"+setting.fields.coord.coordinates{dim-no_RE+indiceRE};
                        actualIndex(indiceRE)=dim-no_RE+indiceRE;
                        if(RE_coords(indiceRE)==obj.settingname)
                            obj.settings.fields.(RE_coords(indiceRE)).index=(dim-no_RE)*(m+1)+actualIndex(indiceRE)-(dim-no_RE);
                            set(obj.handle , 'String' , obj.settings.getSetting(obj.settingname).toString()); 
                            obj.settings.fields.(RE_coords(indiceRE)).index=actualIndex(indiceRE);
                            return;
                        end
                    end
               %if it's not a RE-coord
               set(obj.handle , 'String' , obj.settings.getSetting(obj.settingname).toString()); 
               return;
            end
            %-_-_-_-_-_-_%

            set(obj.handle , 'String' , obj.settings.getSetting(obj.settingname).toString()); 
        end
        
        function destructor(obj)
            delete(obj.eventlistener);
            delete(obj);
        end
        
        function e = Extent(obj)
           e = obj.handle.Extent; 
        end
        function set(obj, varargin)
           %disp(varargin{2});
           set(obj.handle, varargin{:}); 
        end
    end
    
end
