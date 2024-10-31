function replaced_matrix = replace_zeros_with_mean(matrix)
    % Initialize the output matrix
    replaced_matrix = matrix;
    
    % Iterate over each column of the matrix
    for col = 1:size(matrix, 2)
        % Find indices of zero values in the current column
        zero_indices = find(matrix(:, col) < 1e-3);
        
        % Iterate over each zero value found in the current column
        for i = 1:length(zero_indices)
            index = zero_indices(i);
            
            % Find the nearest non-zero elements before and after the zero value
            prev_index = find(matrix(1:index, col) ~= 0, 1, 'last');
            next_index = find(matrix(index:end, col) ~= 0, 1, 'first') + index - 1;
            
            % Calculate the mean of the nearest non-zero elements
            if ~isempty(prev_index) && ~isempty(next_index)
                mean_value = (matrix(prev_index, col) + matrix(next_index, col)) / 2;
            end
            if  isempty(prev_index) && ~isempty(next_index)
                mean_value = matrix(next_index, col);  
            end
            if ~isempty(prev_index) && isempty(next_index)
                mean_value = matrix(prev_index, col);
            end

                % Replace the zero value with the mean
                replaced_matrix(index, col) = mean_value;
            
        end
    end
end



