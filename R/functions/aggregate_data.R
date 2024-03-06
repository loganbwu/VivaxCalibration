aggregate_data = function(data) {
  data_init(data_sim,
            ts = aggregate_time(data_sim$ts),
            cases = aggregate_cases(data_sim$cases),
            n_times = length(aggregate_time(data_sim$ts)))
}
