function subsystems = getSubSystem_str(model)
subsystems = cell(length(model.rxns),1);
for i=1:length(model.rxns)
    if iscell(model.subSystems(i)) & ~isempty(model.subSystems{i})
        subsystems{i} = join(model.subSystems{i},'//');
    else
        subsystems{i} = '';
    end
end
end