# process csv data 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NB assuming 28 days in February
const MONTHDAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

const SCOTLANDREGIONS = [  # not exported
    "Aberdeen City", 
    "Aberdeenshire", 
    "Angus Council", 
    "Argyll and Bute Council", 
    "Clackmannanshire", 
    "Dumfries and Galloway", 
    "Dundee City Council", 
    "East Ayrshire Council", 
    "East Dunbartonshire Council", 
    "East Lothian Council", 
    "East Renfrewshire Council", 
    "Edinburgh", 
    "Falkirk", 
    "Fife", 
    "Glasgow City",  
    "Highland Council", 
    "Inverclyde", 
    "Midlothian", 
    "Moray", 
    "Na h-Eileanan an Iar", 
    "North Ayrshire Council", 
    "North Lanarkshire",  
    "Orkney", 
    "Perth and Kinross", 
    "Renfrewshire", 
    "Scottish Borders", 
    "Shetland Islands", 
    "South Ayrshire Council", 
    "South Lanarkshire",  
    "Stirling", 
    "West Dunbartonshire Council", 
    "West Lothian", 
]

"""
    processagedata(rawfilename[, processedfilename])

Process data from CSV file on age-stratified respiratory syncytial virus data.

If `processedfilename` is supplied and is a file in the `exp_pro` folder, this will 
    be loaded. Otherwise, will use `rawfilename` from the `exp_raw` folder. The 
    processed version will then be saved as `processedfilename`.

Variables should be strings of filenames, including `.csv` at the end.
"""
function processagedata(filename)
    data = CSV.read(filename, DataFrame)

    # Select rows that describe respiratory syncytial virus
    subset!(data, :Pathogen => x -> x .== "Respiratory syncytial virus")

    # Convert dates to a proportion of the year 
    processcsvdates!(data, :WeekBeginning)

    # Set dates to start in April
    insertcols!(data, :Offsetdate => data.Date .- (MONTHDAYS[4] / 365))
    insertcols!(data, :Year => round.(Int, data.Offsetdate, RoundDown))
    insertcols!(data, :FractionDate => data.Offsetdate .- data.Year)

    insertcols!(data, :_Row => axes(data, 1))
    for (i, a) ∈ enumerate(unique(data.AgeGroup))
        d = subset(data, :AgeGroup => x -> x .== a)
        cumulativecases = Vector{Float64}(undef, size(d, 1))
        cumulativecases[1] = d.RatePer100000[1]
        for i ∈ axes(d, 1)
            i == 1 && continue
            if d.Year[i] == d.Year[i-1]
                cumulativecases[i] = d.RatePer100000[i] + cumulativecases[i-1]
            else 
                cumulativecases[i] = d.RatePer100000[i]
            end 
        end 
        insertcols!(d, :CumulativeCases => cumulativecases)
        select!(d, :_Row, :CumulativeCases)
        rename!(d, :CumulativeCases => "Rate$a")
        leftjoin!(data, d, on = :_Row)
    end

    return data 
end

function processagedata(rawfilename, processedfilename)
    return processdata(processagedata, rawfilename, processedfilename)
end

"""
    processrsvdata(rawfilename[, processedfilename])

Process data from CSV file on respiratory syncytial virus data.

If `processedfilename` is supplied and is a file in the `exp_pro` folder, this will 
    be loaded. Otherwise, will use `rawfilename` from the `exp_raw` folder. The 
    processed version will then be saved as `processedfilename`.

Variables should be strings of filenames, including `.csv` at the end.
"""
function processrsvdata(filename)
    data = CSV.read(filename, DataFrame)

    # Select the rows that describe respiratory syncytial virus
    subset!(data, :Pathogen => x -> x .== "Respiratory syncytial virus")

    # Convert the dates to a proportion of the year 
    processcsvdates!(data, :WeekBeginning)

    # Rename cases
    rename!(data, :NumberCasesPerWeek => :Cases)

    # Select columns that are needed 
    select!(data, :Date, :gt, :Cases, :Pathogen)

    return data 
end 

function processrsvdata(rawfilename, processedfilename)
    return processdata(processrsvdata, rawfilename, processedfilename)
end

"""
    processrsvdata(rawfilename[, processedfilename])

Process data from the Oxford Covid-19 Government Response Tracker's stringency index.

If `processedfilename` is supplied and is a file in the `exp_pro` folder, this will 
    be loaded. Otherwise, will use `rawfilename` from the `exp_raw` folder. The 
    processed version will then be saved as `processedfilename`.

Variables should be strings of filenames, including `.csv` at the end.
"""
function processcrgtvdata(filename)
    data = CSV.read(filename, DataFrame)

    # Select rows describing Scotland 
    subset!(data, :RegionName => x -> x .== "Scotland")

    # Rename Date
    rename!(data, :Date => :RawDate)

    # Convert the dates to a proportion of the year 
    processcsvdates!(data, :RawDate)

    return data 
end 

function processcrgtvdata(rawfilename, processedfilename)
    return processdata(processcrgtvdata, rawfilename, processedfilename)
end

function processcsvdates!(df, datecolumn)
    startdates = getproperty(df, datecolumn)
    years = @. round(Int, startdates / 10000, RoundDown)
    months = @. round(Int, (startdates - years * 10000) / 100, RoundDown)
    days = @. startdates - years * 10000 - months * 100
    yeardays = @. MONTHDAYS[months] + days
    dates = @. years + yeardays / 365
    insertcols!(df, :Date => dates)

    # Days since 3 October 2016 (the first day in the RSV dataset)
    datadays = @. round(Int, 365 * (years - 2016.76) + yeardays)

    # Call this :gt for consistency with the model outputs 
    insertcols!(df, :gt => datadays)
end

function processdata(func, rawfilename, processedfilename)
    if isfile(datadir("exp_pro", processedfilename))
        data = CSV.read(datadir("exp_pro", processedfilename), DataFrame)
    else 
        data = func(datadir("exp_raw", rawfilename))
        CSV.write(datadir("exp_pro", processedfilename), data)
    end 
    return data
end

function printrawdate(rawdate::Int)
    stringdate = "$rawdate"
    printrawdate(stringdate)
end 

printrawdate(stringdate::String) = "$(stringdate[7:8])/$(stringdate[5:6])/$(stringdate[1:4])"

_elementwisebetareduction(::Missing) = missing 
_elementwisebetareduction(x::Int) = (100 + x) / 100

function processmobilitydata(filenames...)
    mobilitydata = CSV.read(datadir("exp_raw", filenames[1]), DataFrame)

    for (i, name) in enumerate(filenames) 
        i == 1 && continue 
        md = CSV.read(datadir("exp_raw", name), DataFrame)
        mobilitydata = vcat(mobilitydata, md)
    end

    select!(
        mobilitydata, 
        :date, 
        :sub_region_1, 
        :workplaces_percent_change_from_baseline,
        :transit_stations_percent_change_from_baseline,
        :retail_and_recreation_percent_change_from_baseline,
    )

    # areas of Scotland 
    filter!( 
        :sub_region_1 => x -> !ismissing(x) && x in SCOTLANDREGIONS,
        mobilitydata
    )

    mobilitydata.workplaces_percent_change_from_baseline .= 
        _elementwisebetareduction.(mobilitydata.workplaces_percent_change_from_baseline)
    mobilitydata.transit_stations_percent_change_from_baseline .= 
        _elementwisebetareduction.(
            mobilitydata.transit_stations_percent_change_from_baseline
        )
    mobilitydata.retail_and_recreation_percent_change_from_baseline .= 
        _elementwisebetareduction.(
            mobilitydata.retail_and_recreation_percent_change_from_baseline
        )

    combinedmobilitydata = DataFrame(:date => unique(mobilitydata.date))
    insertcols!(
        combinedmobilitydata,
        :betareduction => [
            (
                tdf = filter(:date => x -> x == combinedmobilitydata.date[i], mobilitydata);
                mean(
                    skipmissing(
                        [
                            tdf.workplaces_percent_change_from_baseline; 
                            tdf.transit_stations_percent_change_from_baseline; 
                            tdf.retail_and_recreation_percent_change_from_baseline
                        ]
                    )
                )
            )
            for i in axes(combinedmobilitydata, 1)
        ],
        :fractiondate => [
            (
                yr = year(combinedmobilitydata.date[i]);
                mth = month(combinedmobilitydata.date[i]);
                d = day(combinedmobilitydata.date[i]);
                yd = MONTHDAYS[mth] + d;
                yr + yd / 365
            )
            for i in axes(combinedmobilitydata, 1)
        ]
    )
    return combinedmobilitydata
end


