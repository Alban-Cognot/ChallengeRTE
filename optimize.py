from ortools.sat.python import cp_model
from collections import defaultdict

def optimize(data):
    model,output_data = create_model(data)
    print("Model created")
    status = solve_model(model,output_data)
    print(f"Model solved with status {status}")


def create_model(data) -> cp_model.CpModel:

    #Parameters
    resources = data["Resources"]
    seasons = data["Seasons"]
    interventions = data["Interventions"]
    exclusions = data["Exclusions"]
    T = data["T"]
    scenarios_number = data["Scenarios_number"]
    quantile = data["Quantile"]
    alpha = data["Alpha"]


    model = cp_model.CpModel()

    #start times necessary for the output
    start_times = {}
    end_times = {}
    #variables corresponding to Delta[start_time]
    durations = {}

    output_data = {}
    output_data["start_times"] = start_times
    output_data["end_times"] = end_times

    
    for (name,intervention) in interventions.items():


        start_time = model.NewIntVar(1,intervention["tmax"],f"start_time_{name}")

        possible_durations = intervention["Delta"]

        #All interventions must be completed at t = T+1
        possible_end_times = cp_model.Domain.FromValues([i+delta+1 for (i,delta) in enumerate(possible_durations) if i+delta+1 <= T+1])


        end_time = model.NewIntVarFromDomain(possible_end_times,f"end_time_{name}")
        

        duration = model.NewIntVarFromDomain(cp_model.Domain.FromValues(possible_durations),f"duration_{name}")

        # duration = possible_durations[start_time]
        model.AddElement(start_time-1,possible_durations,duration)

        model.Add(start_time + duration == end_time).WithName(f"interval_{name}")

        start_times[name] = start_time
        end_times[name] = end_time
        durations[name] = duration

    #Workload contraints and risk calculations

    total_resources = defaultdict(list)
    output_data["total_ressources"] = total_resources
    for t in range(1,T+1):
        time_string = str(t)
        for c_name,resource in resources.items():
            r_c_t = []
            for i_name,intervention in interventions.items():
                if c_name in intervention["workload"].keys() and time_string in intervention["workload"][c_name].keys():
                    workload_i_c_t = intervention["workload"][c_name][time_string]
                    #Create an array to have workload_i_c_t[t'] = workload
                    t_primes = sorted(list(int(k) for k in workload_i_c_t.keys()))
                    offset = min(t_primes)
                    
                    possible_wl_necessary = cp_model.Domain.FromValues(workload_i_c_t.values())
                    workload = model.NewIntVarFromDomain(possible_wl_necessary,f"workload_{i_name}_{c_name}_{t}")
                    model.AddElement(start_times[i_name]-offset,t_primes,workload)
                    r_c_t.append(workload)
            workload_sum = cp_model.LinearExpr.Sum(r_c_t)
            #Upper and Lower limit per resource constraint
            model.Add(workload_sum >= resource["min"][t-1])
            model.Add(workload_sum <= resource["max"][t-1])

            total_resources[c_name].append(workload_sum)

    interval_vars = []

    #Exclusion constraints
    for name,exclusion in exclusions.items():

        intervention_a = exclusion[0]
        intervention_b = exclusion[1]

        start_a = start_times[intervention_a]
        end_a = end_times[intervention_a]

        start_b = start_times[intervention_b]
        dur_b = durations[intervention_b]
        end_b = end_times[intervention_b]

        exclusion_season = exclusion[2]
        non_overlap_times = seasons[exclusion_season]
        
        forbidden_end = non_overlap_times[0]

        time_index = 0

        #loop on the different continious intervals 
        while(time_index != len(non_overlap_times)-1):
            forbidden_start = non_overlap_times[time_index]
            #get the forbidden interval
            while(time_index+1 < len(non_overlap_times) and non_overlap_times[time_index] == non_overlap_times[time_index+1]-1):
                time_index += 1
        forbidden_end = non_overlap_times[time_index]
        forbidden_dur = forbidden_end - forbidden_start


        # --- PROJECT A onto forbidden window ---
        overlap_start_a = model.NewIntVar(forbidden_start, forbidden_end, f"ov_start_a_{name}_from{forbidden_start}_to_{forbidden_end}")
        model.AddMaxEquality(overlap_start_a, [start_a, forbidden_start])

        overlap_end_a = model.NewIntVar(forbidden_start, forbidden_end, f"ov_end_a_{name}_from{forbidden_start}_to_{forbidden_end}")
        model.AddMinEquality(overlap_end_a, [end_a, forbidden_end])

        # overlap duration (>=0)
        overlap_dur_a = model.NewIntVar(0, forbidden_dur, f"ov_dur_a_{name}_from{forbidden_start}_to_{forbidden_end}")
        model.Add(overlap_end_a - overlap_start_a == overlap_dur_a)

        # # presence boolean: true iff overlap_dur_a >= 1 (i.e., A actually covers some part)
        # a_in_f = model.NewBoolVar(f"a_in_forbidden_{name}_from{forbidden_start}_to_{forbidden_end}")
        # model.Add(overlap_dur_a >= 1).OnlyEnforceIf(a_in_f)
        # model.Add(overlap_dur_a == 0).OnlyEnforceIf(a_in_f.Not())

        # optional interval representing A's projection inside the forbidden window
        interval_a_in_f = model.NewIntervalVar(
            overlap_start_a, overlap_dur_a, overlap_end_a, f"interval_a_in_f_{name}_from{forbidden_start}_to_{forbidden_end}"
        )

        
        interval_b = model.NewIntervalVar(
            start_b, dur_b, end_b, f"interval_b_in_f_{name}_from{forbidden_start}_to_{forbidden_end}"
        )

        interval_vars.append({"start" : overlap_start_a, "end" : overlap_end_a, "duration" : overlap_dur_a})
        interval_vars.append({"start" : start_b, "end" : end_b, "duration" : dur_b})

        # --- Forbid overlap inside the forbidden window only ---
        model.AddNoOverlap([interval_a_in_f, interval_b])
        

    # Objective expressions and function
    # obj1 = 1/T * sum_avg_risks

    # obj2 = 1/T * sum_Excess

    # obj = alpha*obj1 + (1-alpha)*obj2
    # model.Minimize(obj)

    return model,output_data

def solve_model(model:cp_model.CpModel,output_data) -> str:

    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    if status in [cp_model.FEASIBLE,cp_model.OPTIMAL]:
        for i in range(len(output_data["start_times"])):
            print(f"I{i+1} : {solver.Value(output_data["start_times"][i])} -> {solver.Value(output_data["end_times"][i])}")
        

    return solver.StatusName()
