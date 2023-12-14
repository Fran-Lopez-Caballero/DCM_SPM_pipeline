function parsave(fullfilename,variable_to_save, name_var_save) %#ok<*INUSL>
    eval([name_var_save ' = variable_to_save;']);
    eval(['save(''' fullfilename ''',''' name_var_save ''');'])
end