
using Plots
using DataFrames
using Dates
using CSV
using NetCDF
using Chain
using LsqFit
using Statistics
using JSON

function load_aossmps(ss, ts)
    dp = "../level3dataall/aossmps/"

    files = @chain begin
        readdir(dp)
        filter(f -> split(f, ".")[end] == "nc", _)
        filter(_) do f
            d = Date(split.(f, ".")[3], dateformat"yyyymmdd")
            return (d >= ss) & (d <= ss .+ Day(2))
        end
    end

    ff = readdir(dp)
    S = mapfoldl(f -> ncread(dp * f, "dN_dlogDp") ./ log(10), hcat, files)
    Dp = ncread(dp * ff[1], "diameter_mobility")
    Des = ncread(dp * ff[1], "diameter_mobility_bounds")
    Nt = mapfoldl(f -> ncread(dp * f, "total_N_conc"), vcat, files)
    bt = mapfoldl(f -> ncread(dp* f, "base_time"), vcat, files)
    tnix =
        mapfoldl(i -> ncread(dp* ff[i], "time_offset") .+ bt[i], vcat, 1:length(files))
    t = @chain map(unix2datetime, tnix) round.(_, Dates.Second)

    return (Dp=Dp, t=t, S=S, Nt=Nt)
end

function load_cpcu(ss, ts)
    dp = "../level3dataall/cpcu/"

    files = @chain begin
        readdir(dp)
        filter(f -> split(f, ".")[end] == "nc", _)
        filter(_) do f
            d = Date(split.(f, ".")[3], dateformat"yyyymmdd")
            return (d >= ss) & (d <= ss)
        end
    end

    ff = readdir(dp)
    N = mapfoldl(f -> ncread(dp * f, "concentration"), hcat, files)
    bt = mapfoldl(f -> ncread(dp * f, "base_time"), vcat, files)
    tnix =
        mapfoldl(i -> ncread(dp * files[i], "time_offset") .+ bt[i], vcat, 1:length(files))
    t = @chain map(unix2datetime, tnix) round.(_, Dates.Second)
    df =  DataFrame(t = t, N = N) 

    return (df = df,)
end

function gf2kappa(gf, RH, Dd)
    aw = RH / 100.0 / exp(2.1 / (Dd * gf))
    return (gf^3.0 - 1.0) * (1.0 - aw) / aw
end

# df = CSV.read("../level3dataall/hou120ncsuhtdmainvertedM1.b1.20220720.csv", DataFrame; header=false)
df = CSV.read("../flaggedlevel3/invertedflaggedHTDMA.csv", DataFrame)
# allowmissing!(df)

Dds = unique(df[!, 3])

function get_mode(myDd)
    df1 = @chain df begin
        filter(:Column3 => D -> D .== myDd, _)
        filter(:Column1 => x -> x .== "invertTDMA", _)
    end

    gfs = @chain df begin
        filter(:Column3 => D -> D .== Dds[1], _)
        filter(:Column1 => x -> x .== "Midpoints", _)
        _[1, 8:end]
        Vector
        reverse
    end

    tgf = df1[!, :Column2] .- Hour(5)
    gfM = df1[!, :8:end] |> Matrix |> x -> reverse(x; dims=2) |> transpose
    n, m = size(gfM)
    modegf = map(i -> gfs[findmax(gfM[:, i])[2]], 1:m)
    modek = map(g -> gf2kappa(g, 70.3, myDd .* 1e9), modegf)

    return tgf, modek, modegf, gfM
end

ss = Date(2022, 8, 17)
ts = DateTime(ss):Minute(1):(DateTime(ss)+Day(2)-Second(1))
smps = load_aossmps(ss, ts)
cpcu = load_cpcu(ss, ts)

CPCUT1 =  cpcu.df.t - Hour(5)
CPCUTii = (CPCUT1 .> DateTime(2022,08,17,12,03,00)) .& (CPCUT1 .< DateTime(2022,08,17,12,23,00))
cpcu.df.N[CPCUTii] .= NaN

t = Dates.format.(smps.t .- Hour(5), "yyyy-mm-dd HH:MM:SS") 

ii = (smps.Dp .> 8) .& (smps.Dp .< 110)

jam1 = smps.t .- Hour(5)
# kk = findall((jam1 .> DateTime(2022,08,17,12,03,00)) .& (jam1 .< DateTime(2022,08,17,12,23,00)))

myS = smps.S[ii, :]'
myS[myS.<0] .= NaN

# smps.Nt[kk] .= NaN
# myS[kk,:] .=NaN
tgf50, modek50, modegf50, gM50 = get_mode(Dds[5])
tgf40, modek40, modegf40, gM40 = get_mode(Dds[4])
tgf30, modek30, modegf30, gM30 = get_mode(Dds[3])
tgf20, modek20, modegf20, gM20 = get_mode(Dds[2])

modek20[modek20 .< 0] .= NaN
modek30[modek30 .< 0] .= NaN
modek40[modek40 .< 0] .= NaN
modek50[modek50 .< 0] .= NaN

# event 1 k----
event1ii50 = (tgf50 .>= DateTime(2022,08,17,10,14,0)) .& (tgf50.<= DateTime(2022,08,17,13,0,0) )
event1ii40 = (tgf40 .>= DateTime(2022,08,17,10,14,0)) .& (tgf40.<= DateTime(2022,08,17,13,0,0) )
event1ii30 = (tgf30 .>= DateTime(2022,08,17,10,14,0)) .& (tgf30.<= DateTime(2022,08,17,13,0,0) )
event1ii20 = (tgf20 .>= DateTime(2022,08,17,10,14,0)) .& (tgf20.<= DateTime(2022,08,17,13,0,0) )


modekall = vcat(modek20[event1ii20],modek30[event1ii30],modek40[event1ii40],modek50[event1ii50])
events = filter(!isnan, modekall)
eventk_allnm = mean(events)
eventk_allstdnm = std(events)
# event 1 prior k----
pevent1ii50 = (tgf50 .> DateTime(2022,08,17,08,14,0)) .& (tgf50.< DateTime(2022,08,17,10,14,0) )
pevent1ii40 = (tgf40 .> DateTime(2022,08,17,08,14,0)) .& (tgf40.< DateTime(2022,08,17,10,14,0) )
pevent1ii30 = (tgf30 .> DateTime(2022,08,17,08,14,0)) .& (tgf30.< DateTime(2022,08,17,10,14,0) )
pevent1ii20 = (tgf20 .> DateTime(2022,08,17,08,14,0)) .& (tgf20.< DateTime(2022,08,17,10,14,0) )

modekpall = vcat(modek20[pevent1ii20],modek30[pevent1ii30],modek40[pevent1ii40],modek50[pevent1ii50])
pevents = filter(!isnan, modekpall)
peventk_allnm = mean(pevents)
peventk_allstdnm = std(pevents)


grt = DateTime(2022,8,17,10,14,0):Minute(1):DateTime(2022,8,17,13,0,0) |> collect
len = Dates.value(grt[end] .- grt[1])  / 1000 
gr = 19.0./60
grD = 3.0 .+ gr*(0:1:len) |> collect


data = Dict(("D" => smps.Dp[ii]), ("S" => (myS)),
    ("t" => t), ("Nt" => smps.Nt), 
    ("tgf20" => tgf20), ("k20" => modek20), 
    ("tgf30" => tgf30), ("k30" => modek30),
    ("tgf40" => tgf40), ("k40" => modek40),
    ("tgf50" => tgf50), ("k50" => modek50),
    # ("tcpc" => [!,:t].-Hour(5)), ("Ncpc" => d[!,:N]),
    ("tcpc" => cpcu.df[!,:t].-Hour(5)), ("Ncpc" => cpcu.df[!,:N]),
    ("grt" => grt), ("grD" => grD))
jdata = JSON.json(data)

open("data.js", "w") do file
    write(file, "const x = ")
    JSON.print(file, jdata)
end