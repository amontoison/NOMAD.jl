const nomad_archive = joinpath(@__DIR__, "downloads", "NOMAD.zip")

builddir = @__DIR__

nomad_path = joinpath(builddir,"nomad.3.9.1")

if isdir(nomad_path)
    rm(nomad_path, recursive=true)
end

run(`unzip $nomad_archive -d $builddir`)

nomad_path = joinpath(builddir,"nomad.3.9.1")

ENV["NOMAD_HOME"] = nomad_path

cd(nomad_path)

run(`./configure`)

run(`make`)

hpp_path = joinpath(nomad_path,"hpp")

mkpath(hpp_path)

sgt_src_path = joinpath(nomad_path,"ext/sgtelib/src")

sgtelib_src = readdir(sgt_src_path)

for i=1:length(sgtelib_src)
    file = sgtelib_src[i]
    n = length(file)
    if file[n-3:n]==".hpp" && file!="Exception.hpp"
        cp(joinpath(sgt_src_path,file),joinpath(hpp_path,file))
    end
end

nmd_src_path = joinpath(nomad_path,"src")

nomad_src = readdir(nmd_src_path)

for i=1:length(nomad_src)
    file = nomad_src[i]
    if file=="Defines.hpp"
        println("found strange file")
        file=="defines.hpp"
    end
    n = length(file)
    if file[n-3:n]==".hpp"
        cp(joinpath(nmd_src_path,file),joinpath(hpp_path,file))
    end
end

cd(joinpath(nomad_path,"src"))

run(`make all`)
