
# Duration of immunity without immune boosting

using DrWatson
@quickactivate :ImmuneBoostingODEs

# Generate duration of immunity by assuming waning at rate 1 and no new infections
immunedurationmodel = let 
    p = SirnsParameters(.0, .0, .0, .0, .0, .0, 1.)
    u0 = sirns_u0(.0, .0; equalrs = false, p)
    sol = run_sirns(u0, p, ( 0., 2.5 ))
    modelcompartments(sol, p)
end
