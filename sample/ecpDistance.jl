import Ecp
import NPZ
import Serialization

# Read ECPs
ecp_1 = Serialization.deserialize(ARGS[1])
ecp_2 = Serialization.deserialize(ARGS[2])

filt_params = 0
fitlt_types = Int64
for k in keys(ecp_1)
    println(k)
    global filt_params = length(k)
    global filt_type = typeof(k[1])
    break
end

fmins = Tuple(parse(Int64, ARGS[i]) for i in 3:(3+filt_params - 1))
fmaxs = Tuple(parse(Int64, ARGS[i]) for i in (3+filt_params):(2 + 2 * filt_params))

println("Distance: ", Ecp.distance_ecp(ecp_1, ecp_2, fmins, fmaxs))