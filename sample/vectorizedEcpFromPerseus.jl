import Ecp
import NPZ

# Read perseus-like file
input_data = Ecp.read_input_from_perseus_file(ARGS[1])
# Parse number of steps for each filtration parameter
number_of_steps = Tuple(parse(Int64, ARGS[i]) for i in 3:length(ARGS))

# Compute and vectorize ECP
contributions = length(size(input_data)) == 2 ? Ecp.compute_contributions_2d(input_data) : Ecp.compute_contributions_3d(input_data)
vectorised_ecp = Ecp.vectorize_ecp(contributions, number_of_steps)

# Save numpy array to file (.npy format)
NPZ.npzwrite(ARGS[2], vectorised_ecp)