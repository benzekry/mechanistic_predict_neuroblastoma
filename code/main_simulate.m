table_res         = readtable('mechanistic/df_fit_indiv.csv');
idx_patients      = table_res{:, 'id'};
table_birth_times = table();
for idx_iter = 1:length(idx_patients)
    idx_patient                                            = idx_patients(idx_iter);
    table_birth_times{idx_patient, 'id'}                   = idx_patient;
    birth_times                                            = simulate_f(idx_iter);
    if ~isempty(birth_times)
        table_birth_times{idx_patient, 'birth_time_met_1'} = birth_times(1);
    else
        table_birth_times{idx_patient, 'birth_time_met_1'} = NaN;
    end
end
writetable(table_birth_times, './simulations/table_birth_times.csv')
