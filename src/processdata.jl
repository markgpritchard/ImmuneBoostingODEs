
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Process csv data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

function processmobilitydata(fn1, fn2, fn3)
    mobilitydata = CSV.read(datadir("exp_raw", fn1), DataFrame)
    append!(mobilitydata, CSV.read(datadir("exp_raw", fn2), DataFrame))
    append!(mobilitydata, CSV.read(datadir("exp_raw", fn3), DataFrame))
    filter!(:sub_region_1 => ismissing, mobilitydata)
    insertcols!(
        mobilitydata, 
        :gtdate => Dates.value.(mobilitydata.date .- Date("2016-10-03")) ./ 365 .+ 2016.76
    )
    # reduction is the mean of transit, workplace and retail, as used 
    # by https://doi.org/10.1371/journal.pcbi.1012452
    insertcols!(
        mobilitydata, 
        :rawreduction => [ 
            +(
                mobilitydata.transit_stations_percent_change_from_baseline[i],
                mobilitydata.workplaces_percent_change_from_baseline[i],
                mobilitydata.retail_and_recreation_percent_change_from_baseline[i],
            ) / 3
            for i ∈ axes(mobilitydata, 1)
        ]
    )
    insertcols!(
        mobilitydata, 
        :proportionreduction => @. 1 + mobilitydata.rawreduction * 0.01
    )
    insertcols!(
        mobilitydata, 
        :reduction => [ 
            mobilitydata.proportionreduction[1:6]; 
            rollmean(mobilitydata.proportionreduction, 7) 
        ]
    )
    select!(mobilitydata, :gtdate, :reduction)
    return mobilitydata
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
