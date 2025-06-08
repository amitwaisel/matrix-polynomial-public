function [result_title] = get_result_title(result_name)
    result_title = strrep(result_name, '_', ' ');
    % Capitalize first character and lower case remaining characters.
    if upper(result_title(1)) == 'X'
        result_title = lower(result_title(2:end));
    else
        result_title = [upper(result_title(1)), lower(result_title(2:end))];
    end
    
    if result_title(1) == '1' && result_title(2) == ' ' && lower(result_title(end)) == 's'
        result_title = result_title(1:end-1);
    end
end

