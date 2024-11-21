% Example vectors (as doubles)
large_time_vector = 0:0.1:10;  % 0 to 10 in steps of 0.1
small_time_vector = 0:0.5:10;  % 0 to 10 in steps of 0.5

% Find indices of small_time_vector in large_time_vector
indices = find(ismember(large_time_vector, small_time_vector));

% Display the indices
disp(indices);