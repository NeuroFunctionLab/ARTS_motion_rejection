% Delete the file if exsit

function f=delete_if_exist(varargin)

 for ii=1:nargin
       if exist(varargin{ii}) ~= 0 
           delete(varargin{ii});
       end;
end;
       
