function parforstatus(N, res, init)
    %%% print a status message in a parfor loop at a given resolution.
    
    if nargin < 3
        init = 0;
    end
    if init
        f = fopen('parforprogress.tmp', 'w');
        fclose(f);
        return
    else
        % add count
        f  = fopen('parforprogress.tmp', 'a');
        fprintf(f, '1\n');
        fclose(f);
    end
    
    % chunck sizes corresponding to the desired resolution and the given
    % data length
    sz = max(round(N*res), 1);
    
    % number of precesses completed
    f  = fopen('parforprogress.tmp', 'r');
    [~, count] = fscanf(f, '%d');
    fclose(f);
    if mod(count,sz) == 0
        fprintf('%4.1f%% complete.\n', (count/N)*100)
    end
    
end