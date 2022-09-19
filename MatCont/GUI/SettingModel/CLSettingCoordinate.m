classdef CLSettingCoordinate < CLSetting

    properties
       index
       coordmodel
    end
    
    methods
        function obj = CLSettingCoordinate(coordmodel, index, subcat)
            if nargin < 3
                subcat = 1;
            end
            coordinate = coordmodel.coordinates{index};
            obj = obj@CLSetting(coordinate, coordmodel.value(index), InputRestrictions.NUM, 2, subcat, index, '~~~~');
            obj.index = index;
            obj.coordmodel = coordmodel;
            
        end
        
        function [valid, msg] = setValue(obj, newvalue) %è qui parte 2...
            global settings;
            [valid, msg] = obj.validitycheck.validate(newvalue);
            if valid
                if(settings.fields.system.sys_type=="DDE")
                    m=settings.fields.system.no_discretizationPoints;
                    DDE_no=settings.fields.system.dim-settings.fields.system.no_RE;
                    RE_no=settings.fields.system.no_RE;
                    if(obj.index>DDE_no) %allora la coordinata modificata è di una RE
                        offsetCoord=DDE_no*(m+1);
                                                                              %se non c'è scritto end c'è un motivo non cambiarlo.
                        obj.coordmodel.value((obj.index-DDE_no+offsetCoord:RE_no:((m+1)*DDE_no+m*RE_no))) = newvalue;
                    else
                        obj.coordmodel.value((obj.index:DDE_no:DDE_no*(m+1))) = newvalue;
                    end
                else
                    obj.coordmodel.value(obj.index) = newvalue;
                end
            end
        end
        
       function value = getValue(obj)
           value = obj.coordmodel.value(obj.index);
       end        
      function newobj = copy(~, ~)
           newobj = [];
      end                
        
      function b = isVisible(obj)
         b= obj.coordmodel.isVisible(); 
      end
          
      
    end
end
