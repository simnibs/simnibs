function [TH2,T]=cat_io_struct2table(S,F)
  if ~exist('F','var') || isempty(F) || (ischar(F) && F=='*') || ...
    (iscell(F) && ~all(cellfun('isempty',strfind(F,'*'))))
    F = get_fieldnames(S);
  end

  TH=cell(1,numel(F)); 
  for fi=1:numel(F)
    TH{fi} = strrep(F{fi},'.','_');
  end

  T=cell(numel(S),numel(F)); TH2=TH;
  for si=1:numel(S)
    fields=''; ef=0;
    for fi=1:numel(F)
      try
        eval(sprintf('iF=numel(S(si).%s);',F{fi}));
        eval(sprintf('cF=isnumeric(S(si).%s);',F{fi}));
        if cF && iF>1 && iF<10
          for iFi=1:iF
            TH2{fi+ef}=sprintf('%s_%d',TH{fi},iFi);
            eval(sprintf('T{si,fi+ef} = S(si).%s(iFi);',F{fi}));
            if iFi<iF, ef=ef+1; end   
          end
          if ef>100; return; end
        else
          TH2{fi+ef}=TH{fi};
          eval(sprintf('T{si,fi+ef} = S(si).%s;',F{fi}));
        end
      catch
        fields=sprintf('%s,%s',fields,F{fi});  
 %       fprintf('Miss field %d - ''%s'' in %d!\n',fi,F{fi},si);
        T{si,fi+ef} = [];
        TH2{fi+ef}  = TH{fi};
      end
    end
    if ~isempty(fields)
      fprintf('Miss field [%s] in %d!\n',fields(2:end),si);
    end
  end
end
function fnS = get_fieldnames(S)
  fnS = fieldnames(SN);
  for fnSi=1:numel(fnS)
    if isfield(S,fnS{fnSi}) 
      if RepByEmpty || ~isempty(SN.(fnS{fnSi}))
        if isstruct(SN.(fnS{fnSi})) 
          S.(fnS{fnSi}) = get_fieldnames(S.(fnS{fnSi}),SN.(fnS{fnSi}),RepByEmpty);
        else
          S.(fnS{fnSi}) = SN.(fnS{fnSi});
        end
      end
    else
      S.(fnS{fnSi}) = SN.(fnS{fnSi});
    end
  end
end