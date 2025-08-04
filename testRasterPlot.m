% Generating raster plot visualization of calcium imaging data for Sora

%% load data
dataFile = 'Sample Data-';   % data fiile
dataPath = 'data';                                                  % data path
[~, sheets] = xlsfinfo(fullfile(dataPath, [dataFile, '.xlsx']));    % get sheet info
nSheet = length(sheets);                                            % number of sheets
disp('Loading data...')
dat = {};
header = {};
for i = 1 : nSheet                                                  % load every sheet
    [dat{i}, header{i}] = xlsread(fullfile(dataPath, [dataFile, '.xlsx']), sheets{i});
    header{i} = header{i}(1, :);
    disp(['Sheet [', sheets{i}, '] loaded.']);
end

%% event detection parameters
par.tWin    = 10/60;        % sliding window in minutes
par.thr     = 5.5;            % threshold, number of MAD
par.ref     = 2/60;         % refractory of calcium events in minutes

%% event detection
[t, evt, sig, nCell] = deal({});            % container for results
dat_ = dat;
for i = 1 : nSheet
    disp(['Processing time series from sheet [', sheets{i}, ']...']);
    evt{i} = {}; sig{i} = {};
    t{i} = dat{i}(:, 1);                    % time
    nCell{i} = size(dat{i}, 2)-1;           % number of cells
    for j = 1 : nCell{i}
        x = dat{i}(:, j+1);                 % raw signal
        [evt{i}{j}, sig{i}{j}] = detectExtremeMAD(t{i}, x, par);    % detection
        dat_{i}(:, j+1) = sig{i}{j};
        disp(['Cell #', int2str(j), ' of ', int2str(nCell{i}), ' done.']);
    end
end
 
%% visualization parameters
nVert       = 6;
nScale      = 1/10;
tRange      = [0, 60];
tSigma      = 15/60;        % in minutes
tRes        = 1e3;
tSess       = 10;           % in minutes
rMax        = 3;            % in Hertz
fracRaster  = 0.8;

%% visualization
hf = figure;
tPSTH = linspace(tRange(1), tRange(2), tRes);
for i = 1 : nSheet
    haTrace{i} = subplot(nVert, nSheet, (0:nVert-4)*nSheet+i);
    haRaster{i} = subplot(nVert, nSheet, (nVert-3)*nSheet+i);
    haPSTH{i} = subplot(nVert, nSheet, (nVert-2)*nSheet+i);
    haHist{i} = subplot(nVert, nSheet, (nVert-1)*nSheet+i);
    % compute event rates
    nEvt{i} = zeros(nCell{i}, 2);
    for j = 1 : nCell{i}
        nEvt{i}(j, 1) = sum(evt{i}{j}<tSess);
        nEvt{i}(j, 2) = sum(evt{i}{j}>=tSess);
    end
    nEvt{i} = nEvt{i} / tSess;
    [~, cellOrder] = sort(nEvt{i}(:, 1));
    % traces
    axes(haTrace{i}); hold on;
    for j = 1 : nCell{i}
        plot(t{i}, sig{i}{cellOrder(j)}*nScale+j);
    end
    % rasters
    axes(haRaster{i}); hold on;
%     line([10, 10], [0.5, nCell{i}+0.5], 'Color', 'r');
    for j = 1 : nCell{i}
        % plot(tRange, [j, j], 'k-');
        % if ~isempty(evt{i}{j}), plot(evt{i}{j}, zeros(size(evt{i}{j}))+j, 'r+'); end
        for e = evt{i}{cellOrder(j)}
            line([e, e], [j-fracRaster/2, j+fracRaster/2], 'Color', 'b');
        end
    end
    % psths
    axes(haPSTH{i}); hold on;
    t_ = cell2mat(evt{i});
    PSTH = zeros(size(tPSTH));
    for j = 1 : length(t_), PSTH = PSTH + normpdf(tPSTH, t_(j), tSigma)/nCell{i}; end
    area(tPSTH, PSTH, 'EdgeColor', 'none', 'FaceColor', 'k');
    % hist
    axes(haHist{i}); hold on;
    bar(1:2, mean(nEvt{i}), 'k');
    errorbar(1:2, mean(nEvt{i}), std(nEvt{i}), 'k+');
    % format
    title(haTrace{i}, sheets{i});
    ylabel(haTrace{i}, 'Cell #');
    ylabel(haRaster{i}, 'Cell #');
    ylabel(haPSTH{i}, 'Calcium event rate [event/cell/min]');
    xlabel(haPSTH{i}, 'Time [min]');
    xlabel(haHist{i}, 'Time [min]');
    ylabel(haHist{i}, 'Average event rate [event/cell/min]');
    set(haTrace{i}, 'XLim', tRange, 'yLim', [0, nCell{i}+1], 'Box', 'on');
    set(haRaster{i}, 'XLim', tRange, 'yLim', [0, nCell{i}+1], 'Box', 'on');
    set(haPSTH{i}, 'XLim', tRange, 'yLim', [0, rMax], 'Box', 'on')
    set(haHist{i}, 'XLim', [0, 3], 'XTick', 1:2, 'XTickLabel', {'0 - 10', '10 - 20'}, 'Box', 'on')
end

%% save analysis results
outputFile = [dataFile, ' - processed time series'];
disp(['Saving processed time series data to [', outputFile, '.xlsx]...']);
if exist(fullfile(dataPath, [outputFile, '.xlsx']), 'file')
    delete(fullfile(dataPath, [outputFile, '.xlsx']));
    disp('Output file already exists, now deleted.');
end
% copyfile(fullfile(dataPath, [dataFile, '.xlsx']), fullfile(dataPath, [outputFile, '.xlsx']));
% disp('New output file created.');
for i = 1 : length(sheets)
    % [~, ~, Raw] = xlsread(fullfile(dataPath, [outputFile, '.xlsx']), sheets{i});
    % [Raw{:, :}] = deal('');
    % xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), Raw, sheets{i});
    xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), header{i}, sheets{i}, 'A1');
    xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), dat_{i}, sheets{i}, 'A2');
    disp(['Sheet [', sheets{i}, '] saved.']);
end

outputFile = [dataFile, ' - event rates'];
disp(['Saving event rate data to [', outputFile, '.xlsx]...']);
if exist(fullfile(dataPath, [outputFile, '.xlsx']), 'file')
    delete(fullfile(dataPath, [outputFile, '.xlsx']));
    disp('Output file already exists, now deleted.');
end
% copyfile(fullfile(dataPath, [dataFile, '.xlsx']), fullfile(dataPath, [outputFile, '.xlsx']));
% disp('Copied output file.');
for i = 1 : length(sheets)
    % [~, ~, Raw] = xlsread(fullfile(dataPath, [outputFile, '.xlsx']), sheets{i});
    % [Raw{:, :}] = deal('');
    % xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), Raw, sheets{i});
    xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), header{i}, sheets{i}, 'A1');
    xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), {'from 0 to 10'; 'from 10 to 20'}, sheets{i}, 'A2');
    xlswrite(fullfile(dataPath, [outputFile, '.xlsx']), nEvt{i}', sheets{i}, 'B2');
    disp(['Sheet [', sheets{i}, '] saved.']);
end


