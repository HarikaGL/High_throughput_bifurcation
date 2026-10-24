% ======================= TR_02_scan.m ==================================================
% Which parameter change moves a parameter set from one class to another?
%
% SETS   up to nPer sets of every class in every sampling iteration, drawn at random from
%        the v3 ensemble (Use_For_ML = 1; portable random choice, same in MATLAB/Octave).
% SCAN   every single parameter (moves 1-14), and in Iterations 3-4 the symmetric pairs
%        (k4, k8) (move 15) and (k1, k11) (move 16), is multiplied by 10^d,
%        d in {-1, -0.5, 0.5, 1}; every result is labelled exactly.
% EXACT  the quantity of an exact condition is set to a margin by changing one parameter:
%        x = -0.3 via k8 (Classes 1, 10 of Iteration 1; detached branch appears),
%        y = +0.3 via k6 and via k4 (Classes 3, 9 of Iteration 2; detached branch appears),
%        k4 = k8 = 2 k2 (Class 10 of Iterations 3, 4; branch point appears),
%        log10(KM1/KM5) = +0.3 / -0.3 via k5 (Class 6; -> Class 7 / 8).
% SYMMETRY BREAKING  single-parameter moves of Class 6-8 sets (moves 1-14): is the branch
%        point lost, and does the (x, y) quadrant (criterion 10) predict whether the new
%        set has a detached branch?
%
% Set itList (e.g. [1 3]) before running to split the work; the summary is written when
% the results of all four iterations exist.
% Output (Class_transitions/results/): TR_scan_it<i>.mat, TR_scan.csv, TR_exact_moves.csv,
%        TR_scan_summary.txt (best single move for every observed transition, exact moves,
%        symmetry breaking)
% ======================================================================================
HERE = fileparts(mfilename('fullpath')); if isempty(HERE), HERE = pwd; end
addpath(fullfile(fileparts(HERE), 'support'), fullfile(HERE, 'tr_lib'));
resDir = fullfile(HERE, 'results'); if ~exist(resDir, 'dir'), mkdir(resDir); end
if ~exist('nPer', 'var'), nPer = 120; end
if ~exist('itList', 'var'), itList = 1:4; end
ds = [-1 -0.5 0.5 1];
pn = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8','k4&k8','k1&k11'};

E = tr_load_ensemble(HERE);
for it = itList
    t0 = tic;
    cls = unique(E.class(E.use & E.iteration == it))';
    SC = []; EX = []; SETS = [];
    for c = cls
        pool = find(E.use & E.iteration == it & E.class == c);
        if numel(pool) < 50, continue; end
        pick = pool(sort(tr_pick(numel(pool), nPer, 1300 + 10*it + c)));
        moves = 1:14; if it >= 3, moves = 1:16; end
        rows = cell(numel(pick), 1); exr = cell(numel(pick), 1);
        parfor s = 1:numel(pick)
            p0 = E.p(pick(s), :); R = []; X = [];
            for mv = moves
                for d = ds
                    switch mv
                        case 15, p = tr_move(p0, 'pair', [8 14], d);
                        case 16, p = tr_move(p0, 'pair', [3 4], d);
                        otherwise, p = tr_move(p0, 'mul', mv, d);
                    end
                    [cl, S] = tr_label(p);
                    df = Inf; if ~isempty(S), df = S.detached_fold_ATOT; end
                    Q = tr_conditions(p);
                    R = [R; it, c, pick(s), mv, d, cl, Q.x, Q.y, Q.det, df]; %#ok<AGROW>
                end
            end
            % exact-condition moves: id, type, value
            ex = {};
            if it == 1 && any(c == [1 10]), ex = {1, 'setx', -0.3}; end
            if it == 2 && any(c == [3 9]),  ex = {2, 'sety', 0.3; 3, 'setyk4', 0.3}; end
            if it >= 3 && c == 10,          ex = {4, 'bp', 0.3}; end
            if it == 3 && c == 6,           ex = {5, 'km', 0.3; 6, 'km', -0.3}; end
            for q = 1:size(ex, 1)
                p = tr_move(p0, ex{q, 2}, ex{q, 3});
                [cl, S] = tr_label(p); df = Inf; if ~isempty(S), df = S.detached_fold_ATOT; end
                X = [X; it, c, pick(s), ex{q, 1}, ex{q, 3}, cl, df]; %#ok<AGROW>
            end
            rows{s} = R; exr{s} = X;
        end
        SC = [SC; vertcat(rows{:})]; EX = [EX; vertcat(exr{:})]; %#ok<AGROW>
        SETS = [SETS; repmat([it c], numel(pick), 1), pick, E.row(pick), E.p(pick, :)]; %#ok<AGROW>
        fprintf('iteration %d, class %d: %d sets scanned\n', it, c, numel(pick));
    end
    save(fullfile(resDir, sprintf('TR_scan_it%d.mat', it)), 'SC', 'EX', 'SETS', '-v7');
    fprintf('iteration %d done in %.0f s\n', it, toc(t0));
end

%% ---------------- summary (when all iterations are done) ----------------------------------
if all(arrayfun(@(i) exist(fullfile(resDir, sprintf('TR_scan_it%d.mat', i)), 'file') == 2, 1:4))
    SC = []; EX = []; SETS = [];
    for i = 1:4
        L = load(fullfile(resDir, sprintf('TR_scan_it%d.mat', i)));
        SC = [SC; L.SC]; EX = [EX; L.EX]; SETS = [SETS; L.SETS]; %#ok<AGROW>
    end
    tr_write_csv(fullfile(resDir, 'TR_scan.csv'), 'iter,src,set,move,d,dst,x,y,crit10,detached_fold_ATOT', SC);
    tr_write_csv(fullfile(resDir, 'TR_exact_moves.csv'), 'iter,src,set,move_id,value,dst,detached_fold_ATOT', EX);
    tr_write_csv(fullfile(resDir, 'TR_sets.csv'), ['iter,class,set,row,' strjoin(E.names, ',')], SETS);
    fid = fopen(fullfile(resDir, 'TR_scan_summary.txt'), 'w');
    for f = [1 fid]
        fprintf(f, 'Single-parameter scan: %d labels, %d source sets\n\n', size(SC,1), size(SETS,1));
        % best single move for every (iteration, source, destination) seen in the scan
        fprintf(f, 'Best single move for every transition (fraction of the source sets reaching the target)\n');
        fprintf(f, 'iter src -> dst   move     d      reached   other classes reached by that move\n');
        key = unique(SC(SC(:,6) ~= SC(:,2), [1 2 6]), 'rows');
        BEST = [];
        for q = 1:size(key, 1)
            it = key(q,1); s = key(q,2); t = key(q,3);
            bestF = -1;
            for mv = unique(SC(SC(:,1) == it & SC(:,2) == s, 4))'
                for d = ds
                    k = SC(:,1) == it & SC(:,2) == s & SC(:,4) == mv & SC(:,5) == d;
                    if ~any(k), continue; end
                    fr = mean(SC(k,6) == t);
                    if fr > bestF, bestF = fr; bm = mv; bd = d; kk = k; end
                end
            end
            if bestF < 0.05, continue; end
            dist = SC(kk, 6); u = unique(dist)';
            oth = strjoin(arrayfun(@(v) sprintf('%d: %.0f%%', v, 100*mean(dist == v)), u, 'UniformOutput', false), ', ');
            fprintf(f, '%3d %3d -> %-3d  %-7s %+5.1f   %5.1f%%    %s\n', it, s, t, pn{bm}, bd, 100*bestF, oth);
            BEST = [BEST; it, s, t, bm, bd, bestF]; %#ok<AGROW>
        end
        tr_write_csv(fullfile(resDir, 'TR_best_moves.csv'), 'iter,src,dst,move,d,fraction', BEST);
        % exact moves
        fprintf(f, '\nExact-condition moves (fraction of the source sets in each class after the move)\n');
        nm = {'x = -0.3 via k8', 'y = +0.3 via k6', 'y = +0.3 via k4', 'k4 = k8 = 2 k2', ...
              'log10(KM1/KM5) = +0.3 via k5', 'log10(KM1/KM5) = -0.3 via k5'};
        for id = 1:6
            for it = 1:4
                for s = unique(EX(EX(:,4) == id & EX(:,1) == it, 2))'
                    k = EX(:,4) == id & EX(:,1) == it & EX(:,2) == s; dst = EX(k, 6);
                    u = unique(dst)';
                    fprintf(f, '  %-30s iteration %d, Class %2d (n = %3d): %s\n', nm{id}, it, s, sum(k), ...
                        strjoin(arrayfun(@(v) sprintf('%d: %.1f%%', v, 100*mean(dst == v)), u, 'UniformOutput', false), ', '));
                end
            end
        end
        % detached branch after the x/y moves: exceptions are folds beyond the window
        k = ismember(EX(:,4), [1 2 3]); det = ismember(EX(k,6), [2 4 5]);
        fprintf(f, '  after x/y moves: detached branch in the window in %.1f%% of %d sets; all others have the detached fold above 1e6: %d of %d\n', ...
            100*mean(det), sum(k), sum(EX(k(:) & ~ismember(EX(:,6), [2 4 5]), 7) > 1e6), sum(~det));
        % symmetry breaking
        % moves that break u* = r: every single parameter except k1, k11, k5, k55
        % (KM1 and KM5 do not enter u*, so those moves keep the branch point)
        brk = setdiff(1:14, [3 4 9 10]);
        k = SC(:,1) >= 3 & ismember(SC(:,2), [6 7 8]) & ismember(SC(:,4), brk);
        lost = ~ismember(SC(k,6), [6 7 8]);
        hasDet = ismember(SC(k,6), [2 4 5]); pred = SC(k,9) == 1;
        agree = (hasDet == pred);
        excWin = ~agree & pred & SC(k,10) > 1e6;
        kp = SC(:,1) >= 3 & ismember(SC(:,2), [6 7 8]) & ismember(SC(:,4), [3 4 9 10]);
        fprintf(f, '\nMoves of k1, k11, k5, k55 (symmetry kept) of Class 6-8 sets: branch point kept in %.1f%% of %d\n', ...
            100*mean(ismember(SC(kp,6), [6 7 8])), sum(kp));
        fprintf(f, ['Symmetry breaking (single moves of Class 6-8 sets that break u* = r, %d moves): branch point lost in %.1f%%; ', ...
            'among those, detached branch present iff x < 0 and y > 0 in %.1f%%; ', ...
            'disagreements with the detached fold above 1e6: %d of %d\n'], sum(k), 100*mean(lost), ...
            100*mean(agree(lost)), sum(excWin(lost)), sum(~agree(lost)));
        kk = find(k);
        for mv = brk
            m = SC(kk,4) == mv & lost;
            if ~any(m), continue; end
            u = unique(SC(kk(m), 6))';
            fprintf(f, '  %-5s: %s\n', pn{mv}, strjoin(arrayfun(@(v) sprintf('%d: %.0f%%', v, 100*mean(SC(kk(m),6) == v)), u, 'UniformOutput', false), ', '));
        end
    end
    fclose(fid);
end
