# xyzReader creates a GeoWhizz file from a XYZ file.

"""
    xyz_to_whizz(xyz_filename, whizz_file; ...)

 Create a `geoWhizz` airborne gravity database in netCDF4 format using the data from a
 Geosoft `XYZ` columnar text file.

    Parameters
    ----------
    xyz_filename : String
        The name of the Geosoft XYZ file.
    whizz_file : String
        The name of the geoWhizz file.
    line_style : String, optional
        The name of the line numbering style used by the data acquirer. DEFAULT = "", (unknown)
    missing_value : String, optional
        The missing_value for the data. DEFAULT=-1.0E-64
"""
function xyz_to_whizz(xyz_filename::String; whizz_file="", line_style="", missing_value=-1.0E-64)
    # fspec = FormatSpec("E")
    missing_str = "-1.0E-64"#fmt(fspec, missing_value)

    # Default output file has same name has input, with "nc" suffix.
    if whizz_file == ""
        whizz_file = replace(normalize(xyz_filename, casefold=true), ".xyz" => ".nc")
    end

    # count number of header records, number of flight lines in xyz_file and number of channels
    (num_head_recs, num_lines, num_channels, field_precisions) = xyz_count(xyz_filename)
    
    # get channel names
    channelnames = xyz_channels(xyz_filename, num_head_recs, num_channels)

    # get line numbers and lengths
    (line_ids, num_fids) = xyz_lines(xyz_filename, num_lines)

    # create the whizz file to write the data to.
    create_whizz(whizz_file; line_style=line_style)

    # create the whizz line sub-groups and channels.
    NCDataset(whizz_file, "a") do df
        lines_group = df.group["Lines"]
        for ii = 1:num_lines
            add_line_whizz(lines_group, line_ids[ii], num_fids[ii], line_style=line_style)
            for jj = 1:num_channels
                add_channel_whizz(lines_group, line_ids[ii], channelnames[jj], precision=field_precisions[jj])
            end
        end
        println("\nAdded $num_lines lines, each of $num_channels channels.")
        println("\nLines:\n$line_ids")
        println("\nChannels:\n $channelnames")
    end

    num_lines_saved = 0
    num_fids_in_line = 0
    recs_of_line_read = 0
    curr_line_group = nothing
    data = Float64[]

    NCDataset(whizz_file, "a") do df
        lines_group = df.group["Lines"]
        open(xyz_filename) do xyzf
            for xyz_rec in eachline(xyzf)
                if startswith(xyz_rec, "/")
                    continue
                elseif start_upper(xyz_rec, "LINE") || start_upper(xyz_rec, "TIE")
                    line_number = split(xyz_rec)[2]
                    curr_line_group = lines_group.group[line_number]
                    ii = findfirst(isequal(line_number), line_ids)
                    num_fids_in_line = num_fids[ii]
                    data = zeros(num_fids_in_line, num_channels)
                    fill!(data, missing_value)
                    recs_of_line_read = 0
                else
                    # check & clean xyz_rec
                    rec_str = split(xyz_rec)
                    if length(rec_str) != num_channels
                        println("ERROR - Line $line_number - Number of channel names mis-match to size of data ")
                        println("Record has $(length(rec_str)) elements. Number of channels = $num_channels")
                        println("Error occurred after $recs_of_line_read records in line.")
                        break
                    end
                    for ii = 1:num_channels
                        if rec_str[ii] == "*"
                            rec_str[ii] = missing_str
                        end
                    end
                    rec_values = parse.([Float64], rec_str)
                    recs_of_line_read += 1
                    data[recs_of_line_read,:] = rec_values
                    if recs_of_line_read == num_fids_in_line
                        write_data_whizz(curr_line_group, data, channelnames)
                        num_lines_saved += 1
                        data = Float64[]
                    end
                    if num_lines_saved == num_lines
                        break
                    end
                end
            end
        end
    end
end


"""
    xyz_count(file_name::String)

 Count the number of header records, flight lines and channels in XYZ file.

 Also the decimal precision of each channel in the Geosoft XYZ file.

    Parameters
    ----------
    file_name : String
        The name of the Geosoft XYZ file.

    Returns
    -------
    num_head_recs : Int
        The .
    num_lines : Int
        The .
    num_channels : Int
        The .
    field_precisions : Int
        The .
"""
function xyz_count(file_name::String)
    num_lines = 0
    num_head_recs = 0
    num_channels = 0
    need_first_data_rec = true
    field_precisions = Int64[]
    xyz_filename = splitpath(file_name)[end]

    open(file_name) do fid
        println("Accessing XYZ data in %s.\n\nFirst few records are: $(xyz_filename)")
        for file_line in eachline(fid)
            if num_lines == 0
                # Always useful to see the first few records in the file.
                println(file_line)
            end
            if file_line[1] == '/'
                num_head_recs += 1
            elseif start_upper(file_line, "LINE") || start_upper(file_line, "TIE")
                num_lines += 1
            # Need every value in the row to be a decimal number (no '*' dummies)
            elseif occursin("*", file_line)
                continue
            elseif num_lines > 0 && need_first_data_rec
                first_data_rec = split(file_line)
                num_channels = length(first_data_rec)
                for ii = 1:length(first_data_rec)
                    bits = split(first_data_rec[ii], '.')
                    if bits[1] == first_data_rec[ii]
                        push!(field_precisions, 0)
                    else
                        push!(field_precisions,length(bits[2]))
                    end
                end
                need_first_data_rec = false
            else
                continue
            end
        end
    end
    println("\n  Found $num_head_recs header records")
    println("  Found $num_lines lines")
    println("  Found $num_channels channels\n")
    return (num_head_recs, num_lines, num_channels, field_precisions)
end


"""
    xyz_channels(file_name::String, num_head_recs::Integer, num_channels::Integer)

 Get the names of the channels in a Geosoft XYZ `file_name`. The algorithm checks
 `num_head_recs` header records. If it finds one with a number of words equal to
 the number of channels, `num_channels`, then it returns those words as an array
 of channel names.
"""
function xyz_channels(file_name::String, num_head_recs::Integer, num_channels::Integer; verbose=false)

    channelnames = repeat([""],num_channels)
    header_rec = 0
    open(file_name) do fileid
        for file_line in eachline(fileid)
            temp_names = split(lstrip(file_line, '/'))
            header_rec += 1
            if header_rec > num_head_recs
                channelnames = [""]
                println("Error - can't find header record with $num_channels channel names.")
                break
            elseif length(temp_names) == num_channels
                channelnames = temp_names
                if verbose
                    for ii = 1:length(channelnames)
                        println("$(channelnames[ii])")
                    end
                end
                break
            end
        end
    end
    return channelnames
end


start_upper(my_string::String, start_string::String) = startswith(uppercase(lstrip(my_string)), start_string)


"""
    xyz_lines(file_name::String, num_lines)

 Returns the line numbers (`line_ids`), and the number of fiducials in each line.
 A helper funtion for `xyz_to_whizz`.

    Parameters
    ----------
    whizz_file : String
        The name of the geoWhizz file.

    Returns
    -------
    None
"""
function xyz_lines(file_name::String, num_lines)

    # initialise counter for the number of fids per line and storage for lines, flights and dates
    num_fids = zeros(Int64, num_lines)
    line_ids = String[]
    line_ctr = 1
    # flight_nos = zeros(Int64, num_lines)
    # flight_dates = [[1, 1, 1980] for k in range(num_lines)]
    # flight_no = 0
    # flight_date = [1, 1, 1980]

    open(file_name) do fid
        for file_line in eachline(fid)
            if length(file_line) < 4
                continue
            elseif file_line[1] == '/'
                if file_line[2] != '/'  # already got channel names
                    continue
                elseif start_upper(file_line, "//FLIGHT")
                    flight_no = parse(Int64, split(file_line)[2])
                elseif start_upper(file_line, "//DATE")
                    test = split(file_line)[2]
                    y, m, d = split(test, '/')
                    flight_date = [parse(Int64, d), parse(Int64, m), parse(Int64, y)]
                end
            elseif start_upper(file_line, "LINE") || start_upper(file_line, "TIE")
                push!(line_ids, split(file_line)[2])
                # flight_nos[line_ctr] = flight_no
                # flight_dates[line_ctr] = flight_date
                line_ctr += 1
            # elseif file_line.count('*') != 0
            #     continue  # print('Skip record for dummy')
            else
                num_fids[line_ctr-1] += 1
            end
        end
    end
    return line_ids, num_fids
end


