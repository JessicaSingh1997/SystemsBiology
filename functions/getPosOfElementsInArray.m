function pos = getPosOfElementsInArray(elements, array)

pos = cell2mat(arrayfun(@(x)find(strcmp(x,array)),elements,'UniformOutput',false))';

end


