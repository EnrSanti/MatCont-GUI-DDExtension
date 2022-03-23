%-_-_-_-_-_-_%
%the class contains a static array which represent the various types of
%differential equation systems the user can use
classdef systemType
   properties (Constant)
       type = ["ODE","DDE"];
   end
   methods (Static)
      %given the an index (1..len(type)) returns the content of the
      %specified array cell
      function out = getType(index) 
         out=systemType.type(index);
      end
      %returns the whole array of types
      function out = getSysTypes()
          out=systemType.type;
      end
      %returns given a string in the array, it's position
      function out = getPositionOfType(sys_type)
         out = find(strcmp(systemType.type,sys_type));
      end
   end
end
%-_-_-_-_-_-_%
  