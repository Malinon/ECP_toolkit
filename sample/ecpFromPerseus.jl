import Ecp
import NPZ
import Serialization

# Read perseus-like file
input_data = Ecp.read_input_from_perseus_file(ARGS[1])

# Compute and vectorize ECP
contributions = length(size(input_data)) == 2 ? Ecp.compute_contributions_2d(input_data) : Ecp.compute_contributions_3d(input_data)

# Save ECP contributions to JDL file
Serialization.serialize(ARGS[2], contributions)