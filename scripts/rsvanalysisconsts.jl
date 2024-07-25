
# This script is called by `rsvanalysis.jl` and `rsvanalysisplots.jl`. Relevant `using` 
# lines are in those scripts

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RSV data from Scotland
data = processrsvdata("respiratory_scot.csv", "rsv.csv")

# Age-specific data
agedata = processagedata("respiratory_age.csv", "rsv_age.csv")

# Data from Oxford Covid-19 Government Response Tracker
crgtdata = processcrgtvdata("OxCGRT_compact_subnational_v1.csv", "crgt.csv")

# To avoid splitting outbreaks, count cases from April each year 
let 
    april1value = MONTHDAYS[4] / 365
    offsetdate = data.Date .- april1value
    aprilyear = round.(Int, offsetdate, RoundDown)
    aprilfractiondate = offsetdate - aprilyear
    insertcols!(data, :AprilYear => aprilyear)
    insertcols!(data, :AprilFractionDate => aprilfractiondate)
    # insert cumulative cases since last April 
    cumulativecases = Vector{Float64}(undef, size(data, 1))
    cumulativecases[1] = data.Cases[1]
    for i ∈ axes(data, 1)
        i == 1 && continue
        if data.AprilYear[i] == data.AprilYear[i-1]
            cumulativecases[i] = data.Cases[i] + cumulativecases[i-1]
        else 
            cumulativecases[i] = data.Cases[i]
        end 
    end 
    insertcols!(data, :AprilCumulativeCases => cumulativecases)
end 
println("In the 12 months from 1 April each year")
for y ∈ 2016:2022 
    inds = findall(x -> y <= x < y + 1, data.AprilYear)
    println("    $(sum(data.Cases[inds])) cases in $y")
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load constants and callbacks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const contacts_alllocations = makecontactmatrix(
    datadir("exp_raw", "MUestimates_all_locations_2.csv")
)
const contacts_schools = makecontactmatrix(datadir("exp_raw", "MUestimates_school_2.csv"))
const contacts_notschools = contacts_alllocations .- contacts_schools
@assert minimum(contacts_notschools) >= 0
const contacts_home = makecontactmatrix(datadir("exp_raw", "MUestimates_home_2.csv"))
const contacts_other = makecontactmatrix(
    datadir("exp_raw", "MUestimates_other_locations_2.csv")
) 
const contacts_lockdownschoolopen = (contacts_home .+ contacts_other .+ contacts_schools) .* 0.8
const contacts_lockdownschoolclosed = (contacts_home .+ contacts_other) .* 0.8

const births::Vector{Float64} = [
    55_690,
    57_781,
    60_041,
    59_046,
    58_791,
    58_590,
    58_027,
    56_014,
    56_725,
    55_098,
    54_488,
    52_861,
    51_308,
    49_863,
    46_809,
    47_786,
    46_959,
]

const save_positions = ( false, false )
# callback to age children by one year 
function ageswitch!(integrator)
    for g ∈ 1:15, j ∈ 1:5 
        i = 16 - g 
        integrator.u[j, (i + 1)] += integrator.u[j, i]
        integrator.u[j, i] = 0.0 
    end
end

startholidays!(integrator) = integrator.p.c = contacts_notschools
endholidays!(integrator) = integrator.p.c = contacts_alllocations
lockdownschoolopen!(integrator) = integrator.p.c = contacts_lockdownschoolopen
lockdownschoolclosed!(integrator) = integrator.p.c = contacts_lockdownschoolclosed

function updatebirthrate!(integrator)
    i = round(Int, integrator.t - 2005, RoundDown) 
    integrator.p.ν = births[i] 
end

const agetimes = [ x + 0.581 for x ∈ 2007:2023 ]
const agecb = PresetTimeCallback(agetimes, ageswitch!; save_positions)

const holidayfractions = [ 0.104, 0.255, 0.389, 0.498, 0.759, 0.893, 0.973 ]
const holidaytimes = [ 
    vec([ x + f for x ∈ 2006:2019, f ∈ holidayfractions ]); 
    [ 2020.104, 2021.603, 2021.759, 2021.893, 2021.973 ]; 
    vec([ x + f for x ∈ 2022:2023, f ∈ holidayfractions ]) 
]
const holidaycb = PresetTimeCallback(holidaytimes, startholidays!; save_positions)

const schoolfractions = [ 0.016, 0.123, 0.304, 0.403, 0.619, 0.805, 0.904 ]
const schooltimes = [ 
    vec([ x + f for x ∈ 2006:2019, f ∈ schoolfractions ]); 
    [ 2020.016, 2021.619, 2021.805, 2021.904 ]; 
    vec([ x + f for x ∈ 2022:2023, f ∈ schoolfractions ]) 
]
const schoolcb = PresetTimeCallback(schooltimes, endholidays!; save_positions)

const lockdownschoolopentimes = [ 
    2020.219, 2020.628, 2020.805, 2020.904, 2021.203, 2021.304, 2021.403 
]
const lockdownschoolopencb = PresetTimeCallback(
    lockdownschoolopentimes, lockdownschoolopen!; 
    save_positions
)

const lockdownschoolclosedtimes = [ 
    2020.221, 2020.759, 2020.893, 2020.973, 2021.255, 2021.389, 2021.498 
]
const lockdownschoolclosedcb = PresetTimeCallback(
    lockdownschoolclosedtimes, lockdownschoolclosed!; 
    save_positions
)

const yearstime = [ x + 0.581 for x ∈ 2007:2022 ]
const yearscb = PresetTimeCallback(yearstime, updatebirthrate!; save_positions)

# for convenience, set cumulative cases to 0 at the start of the analysis 
function resetcumcases!(integrator)
    for i ∈ 1:28
        integrator.u[6, i] = 0.0 
    end
end
const resetcumcasescb = PresetTimeCallback(2006.581, resetcumcases!; save_positions)

cbs = CallbackSet(
    agecb, holidaycb, schoolcb, lockdownschoolopencb, lockdownschoolclosedcb, yearscb
)

# mortality estimated from mid-point of age ranges for women in 2020 from 
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/mortalityratesqxprincipalprojectionscotland
const mortality = [
    506.44,  # 0
    29.53,  # 1
    21.56,  # 2
    7.04,  # 3
    14.07,  # 4
    10.38,  # 5
    13.39,  # 6
    3.25,  # 7
    3.15,  # 8
    3.15,  # 9
    19.32,  # 10
    6.46,  # 11
    12.64,  # 12
    18.90,  # 13
    9.40,  # 14
    26.27,  # 17
    32.97,  # 22
    47.73,  # 27
    71.67,  # 32
    90.65,  # 37
    196.06,  # 42
    276.98,  # 47
    428.51,  # 52
    723.39,  # 57
    1042.47,  # 62
    1625.24,  # 67
    2761.07, # 72
    100_000 / 14  # life expetancy at 75 in UK is 14 years 
    # https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/articles/lifeexpectancycalculator/2019-06-07
] ./ 100_000

const pops = [ 
    # https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
    0,
    54_700,
    54_100,
    52_600,
    51_800,
    52_800,
    53_700,
    56_300,
    57_400,
    59_100,
    58_700,
    59_500,
    61_500,
    63_000,
    65_400,
    328_700,
    339_100,
    310_300,
    317_200,
    384_800,
    405_100,
    377_900,
    335_400,
    344_800,
    279_600,
    243_400,
    212_400,
    168_500 + 118_600 + 62_800 + 32_400
]

const saveat = collect(2016.7568:0.019165:2023.5233)

const intagedatamatrix = let 
    agedatamatrix = zeros(Int(size(agedata, 1) / 7), 7)
    for i ∈ axes(agedatamatrix, 1), j ∈ axes(agedatamatrix, 2)
        agedatamatrix[i, j] = agedata[(7 * (i - 1) + j), (j + 11)]
    end
    # this is per 100_000 -- convert to numbers using the population numbers used above 
    agedatamatrix[:, 1] .*= (55_690 / 100_000)
    agedatamatrix[:, 2] .*= (sum(@view pops[2:5]) / 100_000)
    agedatamatrix[:, 3] .*= (sum(@view pops[6:15]) / 100_000)
    agedatamatrix[:, 4] .*= (sum(@view pops[16:21]) / 100_000)
    agedatamatrix[:, 5] .*= (sum(@view pops[22:25]) / 100_000)
    agedatamatrix[:, 6] .*= (sum(@view pops[26:27]) / 100_000)
    agedatamatrix[:, 7] .*= (pops[28] / 100_000)
    for k ∈ axes(agedatamatrix, 1) 
        i = size(agedatamatrix, 1) + 1 - k 
        i == 1 && continue
        if agedatamatrix[i, 1] >= agedatamatrix[(i - 1), 1]
            for j ∈ axes(agedatamatrix, 2)
                agedatamatrix[i, j] += -agedatamatrix[(i - 1), j]
            end
        end
    end
    @assert minimum(agedatamatrix) >= 0

    round.(Int, agedatamatrix)
end
