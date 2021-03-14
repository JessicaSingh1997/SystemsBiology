function [ecuaciones,ecuaciones2]=getRxn_cobraFormat(model,pos,names)
if nargin<2 || isempty(pos); pos = 1:length(model.rxns); end;
if nargin<3 || isempty(names); names = 0; end;
if names==1; 
    nombreCompuestos = model.metNames;
elseif names==2
    nombreCompuestos = model.metFormulas;
else
    nombreCompuestos=model.mets;
end
if isfield(model,'metNames')
    nombreCompuestos2=model.metNames;
else
    nombreCompuestos2=model.mets;
end

if iscell(pos)
    pos = cell2mat(arrayfun(@(x)find(strcmp(x,model.rxns)),pos,'UniformOutput',false))';
elseif ischar(pos)
    pos = find(strcmp(model.rxns, pos));
end

ecuaciones=[];
ecuaciones2=[];
for i=1:length(pos)
    reversible=false;
    if (model.lb(pos(i))<0 && model.ub(pos(i))>0) || (model.lb(pos(i))>0 && model.ub(pos(i))<0)
        reversible=true;
    end
    
    %Se encuentran las posiciones de los compuestos de la reaccion i
    PosicionCompuestos=find(model.S(:,pos(i)));
    
    %Se crea la ecuaciÛn que describe la reacciÛn
    reactantes='';
    productos='';
    reactantes2='';
    productos2='';
    YaTieneUnReactante=false;
    YaTieneUnProducto=false;
    
    %Se va agregando cada compuesto a la secciÛn de la ecuaciÛn que
    %corresponda (reactantes o productos)
    for j=1:length(PosicionCompuestos)
        if model.S(PosicionCompuestos(j),pos(i))<0
            if YaTieneUnReactante
                reactantes=[reactantes ' + ' num2str((model.S(PosicionCompuestos(j),pos(i)))*-1) ' ' nombreCompuestos{PosicionCompuestos(j)} ];
                reactantes2=[reactantes2 ' + ' num2str((model.S(PosicionCompuestos(j),pos(i)))*-1) ' ' nombreCompuestos2{PosicionCompuestos(j)} ];
            else
                reactantes=[reactantes num2str((model.S(PosicionCompuestos(j),pos(i)))*-1) ' ' nombreCompuestos{PosicionCompuestos(j)} ];
                reactantes2=[reactantes2 num2str((model.S(PosicionCompuestos(j),pos(i)))*-1) ' ' nombreCompuestos2{PosicionCompuestos(j)} ];
                YaTieneUnReactante=true;
            end
        else
            if YaTieneUnProducto
                productos=[productos ' + ' num2str(model.S(PosicionCompuestos(j),pos(i))) ' ' nombreCompuestos{PosicionCompuestos(j)} ];
                productos2=[productos2 ' + ' num2str(model.S(PosicionCompuestos(j),pos(i))) ' ' nombreCompuestos2{PosicionCompuestos(j)} ];
            else
                productos= [productos num2str(model.S(PosicionCompuestos(j),pos(i))) ' ' nombreCompuestos{PosicionCompuestos(j)} ];
                productos2= [productos2 num2str(model.S(PosicionCompuestos(j),pos(i))) ' ' nombreCompuestos2{PosicionCompuestos(j)} ];
                YaTieneUnProducto=true;
            end
        end
    end
    
    if reversible
        ecuacion=[reactantes ' <=> ' productos];
        ecuacion2=[reactantes2 ' <=> ' productos2];
    else
        ecuacion=[reactantes ' -> ' productos];
        ecuacion2=[reactantes2 ' -> ' productos2];
    end
    ecuaciones=[ecuaciones;{ecuacion}];
    ecuaciones2=[ecuaciones2;{ecuacion2}];
end

end

