module Measurement

using CSV
using MAT
using FileIO
using JLD2

export data 

function load_data()
    directory = raw"Z:\Daten\Netzselbstregeleffekt\UW_Osterburken"  
    files = readdir(directory)
    paths = String[]
    for f in files 
        path = joinpath(directory,f)
       push!(paths,path)
    end
    file = matread(paths[1])
    return file
end
 
const data = load_data()

end

