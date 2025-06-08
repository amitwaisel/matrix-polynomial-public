function [result_name] = get_result_name(str)
    result_name = str;
    result_name = strrep(result_name, '.', '_');
    result_name = strrep(result_name, ' ', '_');
end

