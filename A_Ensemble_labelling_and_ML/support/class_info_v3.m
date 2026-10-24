function C = class_info_v3()
%CLASS_INFO_V3  Names, plot labels and colours of the ten bifurcation classes (v3).
%   One place for the class definitions, used by every script and helper.
%
%     #   name                  group              rule (window ATOT in cfg.ATOT_window)
%     1   LP_S_ONLY             limit points       two folds on the main branch, S-shaped,
%                                                  no fold on a detached branch
%     2   LP_S_DETACHED         limit points       as 1, plus a fold on a detached branch
%     3   LP_Z_ONLY             limit points       two folds on the main branch, Z-shaped,
%                                                  no fold on a detached branch
%     4   LP_Z_DETACHED         limit points       as 3, plus a fold on a detached branch
%     5   LP_DETACHED_ONLY      limit points       fold only on a detached branch
%     6   BP_ONLY               branch points      branch point, no fold (pitchfork)
%     7   LP_BP_LPgreaterBP     branch points      branch point + fold, x2(LP) > x2(BP)
%     8   LP_BP_LPlessBP        branch points      branch point + fold, x2(LP) < x2(BP)
%     9   NO_SING_biphasic      no singularity     one interior maximum of x2
%    10   NO_SING_monotone      no singularity     monotone
%    11   OTHER                 -                  anything else
%     0   FAILED

C.names = {'LP_S_ONLY','LP_S_DETACHED','LP_Z_ONLY','LP_Z_DETACHED','LP_DETACHED_ONLY', ...
           'BP_ONLY','LP_BP_LPgreaterBP','LP_BP_LPlessBP','NO_SING_biphasic', ...
           'NO_SING_monotone','OTHER'};
C.labels = {'LP S-shaped only','LP S-shaped + detached','LP Z-shaped only', ...
            'LP Z-shaped + detached','LP on detached branch only','BP only (pitchfork)', ...
            'LP + BP, x_2^{LP} > x_2^{BP}','LP + BP, x_2^{LP} < x_2^{BP}', ...
            'No singularity, biphasic','No singularity, monotone','Other'};
% colours of the class wheel figure, used in every figure
C.colors = [0.545 0.820 0.549;   %  1 LP S only (green)  #8BD18C
            0.310 0.651 0.353;   %  2 LP S + detached (dark green)  #4FA65A
            0.549 0.686 0.910;   %  3 LP Z only (blue)  #8CAFE8
            0.310 0.498 0.816;   %  4 LP Z + detached (dark blue)  #4F7FD0
            0.788 0.635 0.478;   %  5 LP detached only (tan)  #C9A27A
            0.910 0.549 0.553;   %  6 BP only (red-pink)  #E88C8D
            0.965 0.741 0.431;   %  7 LP + BP, LP > BP (orange)  #F6BD6E
            0.843 0.557 0.878;   %  8 LP + BP, LP < BP (violet)  #D78EE0
            0.498 0.847 0.847;   %  9 biphasic (teal)  #7FD8D8
            0.741 0.741 0.741;   % 10 monotone (grey)  #BDBDBD
            0.478 0.478 0.478];  % 11 other (dark grey)  #7A7A7A
% v2 (9-class) names, kept for traceability
C.names_v2 = {'LP_ONLY_S-shaped','LP_ONLY_Z-shaped','BP_ONLY','LP_BP_LPgreaterBP', ...
              'LP_BP_LPlessBP','NO_SING_biphasic','LP_DETACHED','NO_SING_monotone','OTHER'};
% original 7 classes
C.names_7 = [C.names_v2(1:6), {'OTHER'}];
C.n = numel(C.names);
end
