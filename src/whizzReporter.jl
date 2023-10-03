# whizzReporter contains the reporting functions for GeoWhizzNC.

"""
    report_whizz(whizz_file; ...)

 Print a summary of the contents of a geoWhizz file.

    Parameters
    ----------
    whizz_file : String
        The name of the geoWhizz file.
    line : String, optional
        The line number of a particular survey line. DEFAULT="". If specified, the attributes
        of the line will be reported.
    channel_name : String, optional
        The name of a particular data channel. DEFAULT="". If specified, the attributes
        of the channel will be reported.
"""
function report_whizz(whizz_file; line="", channel_name="", detailed=false)
    NCDataset(whizz_file) do df
        report_whizz_hdr(whizz_file, df; detailed=detailed)

        lines_group = df.group["Lines"]
        lines = keys(lines_group.group)
        lines_header = @sprintf("%.0f Lines:", size(lines,1))
        distance_flown(df, lines_group)
        println("\nLines:")
        println(lines)
        
        if line == ""
            line = lines[1]
        end
        report_line = lines_group.group[line]
        println("\nLine attributes")
        println(report_line.attrib)   
        vars = keys(report_line)
        channel_header = @sprintf("Line %s; %.0f channels:", line, size(vars,1))
        println(channel_header)
        println(vars)

        if channel_name != ""
            report_channel = report_line[channel_name]
            println("Channel " * channel_name * " attributes")
            println(report_channel)   
        end
    end
end


"""
   report_flights(whizz_file; flight_channel="FLIGHT", lines=[], detailed=false)

 Reports a summary of the flight information in GeoWhizz file, `whizz_file`. 
"""
function report_flights(whizz_file; flight_channel="FLIGHT", lines=[], detailed=false)

    NCDataset(whizz_file) do df
        report_whizz_hdr(whizz_file, df; detailed=detailed)
        lines_group = df.group["Lines"]
        if lines == []
            lines = keys(lines_group.group)
        end
        flight_dict = Dict()
        for line in lines
            this_flight = lines_group.group[line][flight_channel][1]
            if this_flight in keys(flight_dict)
                append!(flight_dict[this_flight], [line])
            else
                flight_dict[this_flight] = [line]
            end
        end
        # sort the flight dict
        println("Flights")
        for flight in keys(flight_dict)
            println("    $flight")
            if detailed
                for line in flight_dict[flight]
                    print("L$line  ")
                end
                print("\n")
            end
        end
    end
end


"""
    report_sampling(whizz_file; time_chan="", n_chan="", e_chan="", lines=[], detailed=false)

 Reports a summary of the sample spacing of the GeoWhizz NCDataset, whizz_file.
"""
function report_sampling(whizz_file; time_chan="", n_chan="", e_chan="", lines=[], detailed=false)

    NCDataset(whizz_file) do df
        report_whizz_hdr(whizz_file, df; detailed=detailed)
        lines_group = df.group["Lines"]
        if lines == []
            lines = keys(lines_group.group)
        end
        num_lines = size(lines, 1)

        if n_chan == ""
            n_chan = df.attrib["northing"]
        end
        if e_chan == ""
            e_chan = df.attrib["easting"]
        end
        if time_chan == ""
            time_chan = df.attrib["time"]
        end

        time_deltas = Vector{Float64}()
        dd_deltas = Vector{Float64}()
        for line in lines
            t = lines_group.group[line][time_chan]
            n = lines_group.group[line][n_chan]
            e = lines_group.group[line][e_chan]
            append!(time_deltas, diff(t))
            dd = _distance.(n, e)
            append!(dd_deltas, abs.(diff(dd)))
        end
        mean_dt = mean(time_deltas)
        min_dt = minimum(time_deltas)
        max_dt = maximum(time_deltas)
        mean_dd = mean(dd_deltas)
        min_dd = minimum(dd_deltas)
        max_dd = maximum(dd_deltas)

        println("Sample time and distance statistics")
        println(@sprintf("  Min  = %.3f s, %.1f m", min_dt, min_dd))
        println(@sprintf("  Max  = %.3f s, %.1f m", max_dt, max_dd))
        println(@sprintf("  Mean = %.3f s, %.1f m", mean_dt, mean_dd))
    end
end


"""
    function report_whizz_hdr(whizz_file, df; detailed=true)
 
 Print out the project (global) attributes for the
 GeoWhizz NCDataSet, `whizz_file`, with fileID, `df`.
"""
function report_whizz_hdr(whizz_file, df; detailed=true)
    println("The geoWhizz Filename:")
    whizz_filename = splitpath(whizz_file)[end]
    println("    ", whizz_filename)
    println("Global Attributes:")
    if detailed
        println(df.attrib)
    else
        println("    Whizz Version: ", df.attrib["Whizz_Version"])
        println("    Project Name: ", df.attrib["project_name"])
        println("    Block Name: ", df.attrib["block_name"])
        println("    Customer: ", df.attrib["customer"])
        println("    Acquirer: ", df.attrib["acquirer"])
        println("    Customer: ", df.attrib["customer"])
    end
    print("\n")
end


"""
    distance(x, y)

 The length of the hypotenuse of the triangle with other side lengths, `x` and `y`.
"""
function _distance(x, y)
    return sqrt(x ^ 2.0 + y ^ 2.0)
end


"""
    linelength(x, y)

 Calculates the length of the line (x,y) where `x` and `y` are vectors of the same length.
"""
function linelength(x, y)
    return sum(sqrt.(diff(x) .^ 2 .+ diff(y) .^2))
end


"""
    distance_flown(dataset, lines_group)

 Calculates the total number of lines, and the cumulative sum of the line lengths in km.
"""
function distance_flown(dataset, lines_group)
    lineDistance = 0.0
    count = 0
    distance_calculable = false        

    if haskey(dataset.attrib, "northing") && haskey(dataset.attrib, "easting")
        distance_calculable = true
        n = dataset.attrib["northing"]
        e = dataset.attrib["easting"]
    end
    lines = keys(lines_group.group)
    for line in lines
        if distance_calculable
            nPos = lines_group.group[line][n]
            ePos = lines_group.group[line][e]
            lineDistance += linelength(nPos, ePos)
        end
        count += 1
    end
    if distance_calculable
        reportstr = @sprintf("%.0f lines: total distance flown = %.3f km.", count, lineDistance/1000.0)
    else
        reportstr = @sprintf("%.0f lines: total distance flown unknown (no easting/northing).", count)
    end            
    println(reportstr)
    return (count, lineDistance/1000.0)
end


"""
    title_string(df)

 Returns a string that can be used as a plot or report title, and built
 from the metadata in the GeoWhizz file, `df`.
"""
function title_string(df)
    my_string = ""
    join = ""
    if df.attrib["project_name"] != ""
        my_string = @sprintf("Project %s", df.attrib["project_name"])
        join = ", "
    end
    if df.attrib["block_name"] != ""
        my_string *= join * @sprintf("Block %s", df.attrib["block_name"])
    end
    join = "\n    "
    if df.attrib["acquirer"] != ""
        my_string *= join * @sprintf("Acquired by %s", df.attrib["acquirer"])
        join = ", "
    end
    if df.attrib["acquirer_projectID"] != ""
        my_string *= join * @sprintf("Acquirer project ID: %s", df.attrib["acquirer_projectID"])
    end
    return my_string
end
