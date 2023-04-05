import math
import copy
import csv

# Represents the costs of something with circuits
# optimizing for low width, low depth, etc.
class Cost:
	def __init__(self, low_depth, low_T, low_width):
		self.low_depth = low_depth
		self.low_T = low_T
		self.low_width = low_width

	def add(self, cost2):
		return Cost(
			low_depth = self.low_depth.add(cost2.low_depth), 
			low_T = self.low_T.add(cost2.low_T),
			low_width = self.low_width.add(cost2.low_width)
		)

	def subtract(self, cost2):
		return Cost(
			low_depth = self.low_depth.subtract(cost2.low_depth), 
			low_T = self.low_T.subtract(cost2.low_T),
			low_width = self.low_width.subtract(cost2.low_width)
		)
	def multiply(self, n):
		return Cost(
			low_depth = self.low_depth.multiply(n), 
			low_T = self.low_T.multiply(n),
			low_width = self.low_width.multiply(n)
		)
	def message(self):
		 message = ""
		 message += "Width-optimal: \n" +self.low_width.message()
		 message += "T-optimal: \n" + self.low_T.message()
		 message += "Depth-optimal: \n" + self.low_depth.message() + "\n"
		 return message

	# Computes the surface code costs for all strategies
	# epsilon is the physical error rate
	# target is the target error rate of the overall algorithm
	# n is the bit-length of the prime
	def surface_code_costs(self, n, epsilon = 0.001, target = 0.01):
		t_surface = self.low_T.surface_code_costs(n, False, epsilon, target)
		depth_surface = self.low_depth.surface_code_costs(n, True, epsilon, target)
		width_surface = self.low_width.surface_code_costs(n, False, epsilon, target)
		return {'T': t_surface, 'Depth' : depth_surface, 'Width': width_surface}



# Contains all relevant cost metrics for a single circuit
class SingleCost:
	def __init__(self, width, T_depth, full_depth, measure, T_count, single_qubit, CNOT):
		self.width = width
		self.T_depth = T_depth
		self.T_count = T_count
		self.full_depth = full_depth
		self.measure = measure
		self.single_qubit = single_qubit
		self.CNOT = CNOT

	# Costs of two sequential circuits
	# Uses the maximum width (assumes circuits are run sequentially)
	def add(self, cost2):
		return SingleCost(
			width = max(self.width, cost2.width),
			T_depth = self.T_depth + cost2.T_depth,
			T_count = self.T_count + cost2.T_count,
			full_depth = self.full_depth + cost2.full_depth,
			measure = self.measure + cost2.measure,
			single_qubit = self.single_qubit + cost2.single_qubit,
			CNOT = self.CNOT + cost2.CNOT,
		)

	# Subtracts
	# Assumes the width actually stays the same (does not decrease!)
	def subtract(self, cost2):
		return SingleCost(
			width = self.width,
			T_depth = self.T_depth - cost2.T_depth,
			T_count = self.T_count - cost2.T_count,
			full_depth = self.full_depth - cost2.full_depth,
			measure = self.measure - cost2.measure,
			single_qubit = self.single_qubit - cost2.single_qubit,
			CNOT = self.CNOT - cost2.CNOT,
		)
	def multiply(self, n):
		return SingleCost(
			width = self.width,
			T_depth = self.T_depth * n,
			T_count = self.T_count * n,
			full_depth = self.full_depth * n,
			measure = self.measure * n,
			single_qubit = self.single_qubit * n,
			CNOT = self.CNOT * n
		)

	# Computes costs in a surface code architecture
	# n = bit-length of the prime
	# low_depth = boolean on whether this is a low-depth implementation
	# epsilon = physical error rate
	# target = target error rate of final algorithm
	def surface_code_costs(self, n, low_depth = False, epsilon=0.001, target=0.01):
		# For negative/0 depth, this ensures the program doesn't crash
		# (happens with blank rows in the data spreadsheet)
		if self.full_depth <= 0:
			distance = 1
		else:
			# Uses the sam formula as Gidney and Ekera 2019
			distance = (math.log(target/0.1)*2 - 2*math.log(self.full_depth) - 2*math.log(self.width) - math.log(100*epsilon))/math.log(100*epsilon)
			distance = math.ceil(distance)
		# Each CNOT takes depth 3 
		# This determines the average rate of CNOTs
		cnot_rate = self.CNOT*3/(self.full_depth)
		# Guesses at how many simultaneous CNOTs might be necessary
		# The factor 4 is a magic constant
		cnot_space = math.ceil(4*cnot_rate)
		# The rate of T gates, and hence the number of factories
		t_factories = math.ceil(self.T_count*6.5/(self.full_depth))
		# The number of buckets to store T-gates
		# The factor 4 is a magic constant
		t_buckets = t_factories*4
		# In the low-depth regime, we know we need at least 3n simultaneous T-gates
		if low_depth:
			t_buckets = max(t_buckets, n*3)
		# Add up how many physical qubits are necessary for different purposes
		physical_qubits = self.width * distance * distance # for logical qubits
		physical_qubits += cnot_space * (2*distance*distance + 4*distance) # for CNOTs
		physical_qubits += t_factories * 72 * distance * distance # for t factories
		physical_qubits += t_buckets * (3*distance*distance + 4*distance) # for t-buckets
		cycles = self.full_depth * distance
		return {'distance': distance, 'cycles': cycles, 'cnot_rate' : cnot_rate, 't_factories' : t_factories, 't_buckets': t_buckets, 'physical_qubits': physical_qubits}


	#Outputs the costs to a string
	def message(self):
		message = ""
		message += "    CNOT: " + str(self.CNOT) + "\n"
		message += "    Single-qubit: " + str(self.single_qubit) + "\n"
		message += "    Measurements: " + str(self.measure) + "\n"
		message += "    T gates: " + str(self.T_count) + "\n"
		message += "    T-depth: " + str(self.T_depth) + "\n"
		message += "    Full depth: " + str(self.full_depth) + "\n"
		message += "    Width: " +str( self.width) + "\n"
		return message

	#A string which can act as a header for a similar csv as Q# outputs
	@classmethod
	def CSV_header(Class):
		return "Size, CNOT, Single Qubit, T gates, R gates, Measurements, T-depth, Initial Width, Extra Width, Full Depth, Window Size, Code Distance, Cycles, CNOT Rate, T Factories, T Buckets, Physical Qubits\n"

	# Outputs the costs in a row that matches the header
	def csv_row(self):
		message = ""
		message += str(self.CNOT) + ", "
		message += str(self.single_qubit) + ", "
		message += str(self.T_count) + ", ,"
		message += str(self.measure) + ", "
		message += str(self.T_depth) + ", ,"
		message += str(self.width) + ","
		message += str(self.full_depth)
		return message

def surface_code_message(surface_code_cost):
	message = ""
	message += str(surface_code_cost['distance']) + ", "
	message += str(surface_code_cost['cycles']) + ", "
	message += str(surface_code_cost['cnot_rate']) + ", "
	message += str(surface_code_cost['t_factories']) + ", "
	message += str(surface_code_cost['t_buckets']) + ", "
	message += str(surface_code_cost['physical_qubits']) + ", "
	return message


# Returns the cost of a lookup for an n-bit elliptic curve point 
# among a table of 2^window_size points
# Extrapolations based on output of Q#
def Lookup_Cost(n, window_size):
	main_exponent = 2**window_size;
	costs = Cost(
		low_T = SingleCost(
			width = 4.07*window_size + 2.54*n - 2.16,
			T_depth = 0, #incorrect
			full_depth = 112*main_exponent* math.log(n,2) + 20.3*main_exponent + 514.2,
			measure = 1.00*main_exponent + 4.00 + 2.0*n,
			T_count = 4.0*main_exponent + 24.0,
			single_qubit = 7.74 * main_exponent + 10.68 + 2.0*n,
			CNOT = main_exponent * (7.40 + 3.0*n) + 835.6
		),
		low_width = SingleCost(
			width = 2.0*window_size + 2.0 + 2.0*n,
			T_depth = 0, # incorrect
			full_depth = 105.1*n + 6.0 * n * main_exponent + 68 * main_exponent + 245,
			measure = main_exponent + 4.0 + 2.0*n,
			T_count = 4.0*main_exponent + 7.07*n - 2.9,
			single_qubit = 7.793 * main_exponent + 5.75 + 2.0*n,
			CNOT = main_exponent*(3.45+1.015*n) + 2665
		),
		low_depth = SingleCost(
			width = 3.84*window_size + 1.34 + 2.534*n,
			T_depth = 0, # incorrect
			full_depth = 108.48*main_exponent + 0.16*main_exponent*window_size + 20.55*main_exponent*math.log(n,2)+534.4,
			measure =  main_exponent + window_size +2.0*n,
			T_count = 4.0*main_exponent+24.0,
			single_qubit = 23*main_exponent + window_size + 2.0*n - 12,
			CNOT = main_exponent*(6.152 + 3.006 * n) + 817.113
		)
	)	
	return costs

# Returns the costs of a single point addition with 
# window size of 8
# Based on estimates from Q#
def point_addition_cost(n):
	costs = {}
	nsquared = n*n
	lgn = math.log(n)/math.log(2.0)
	nlgn = n*lgn
	n2lgn = nsquared*lgn
	costs = Cost(
		low_T = SingleCost(
			width = 10.0*n + 1.5*math.floor(lgn) + 18.9,
			T_depth = 431.6*nsquared + 17572,
			full_depth = 1562 * nsquared + 120830,
			measure = 85*nsquared + 19465,
			T_count = 1182*nsquared + 92166,
			single_qubit = 648*nsquared + 101890,
			CNOT = 2391*nsquared + 473340
		),
		low_width = SingleCost(
			width = 7.99*n + 3.81*math.floor(lgn) + 17.1,
			T_depth = 144.5 * n2lgn + 626302,
			full_depth = 464.6 * n2lgn + 2074976,
			measure = 753.7*n2lgn - 21095,
			T_count = 503.4*n2lgn + 1318387,
			single_qubit = 167.7 * n2lgn + 544865,
			CNOT = 751.2*n2lgn + 2296571
		),
		low_depth = SingleCost(
			width = 11.0*n + 28.6,
			T_depth = 226.1*nlgn + 14469,
			full_depth = 1485*nlgn + 52413,
			measure = 202.5*nsquared - 14509,
			T_count = 2745*nsquared - 85878,
			single_qubit = 1462*nsquared - 35830,
			CNOT = 6481*nsquared+44882
		)
	)	

	return costs

# Load costs from a CSV
# Assumes the data is formatted as surface-code specific data
def load_from_csv(csv_file_name, n):
	csv.register_dialect('cost_csv_dialect', skipinitialspace = True)
	with open(csv_file_name, newline="\n") as csvfile:
		csvCosts = csv.DictReader(csvfile, dialect='cost_csv_dialect')
		costs = SingleCost(0,0,0,0,0,0,0)
		for row in csvCosts:
			if (row['size'] == str(n)):
				costs.CNOT = int(row['CNOT count'])
				costs.single_qubit = int(row['1-qubit Clifford count'])
				costs.T_count = int(row['T count'])
				costs.measure = int(row['M count'])
				costs.T_depth = 0
				costs.width = int(row['extra width'])
				costs.full_depth = int(row['Full depth'])/2 # halved since Q# doubles
	return costs

# Returns the costs of a single point addition with 
# window size of 8 for fixed-modulus curves
# Based on estimates from Q#
def fixed_modulus_point_addition_cost(n):
	low_T_costs = load_from_csv('EllipticCurveEstimates/LowT/Fixed-modulus-signed-all-gates.csv', n)

	low_width_costs = load_from_csv('EllipticCurveEstimates/LowWidth/Fixed-modulus-signed-all-gates.csv', n)

	low_depth_costs = load_from_csv('EllipticCurveEstimates/LowDepth/Fixed-modulus-signed-all-gates.csv', n)

	return Cost(low_T = low_T_costs, low_width = low_width_costs, low_depth = low_depth_costs)



def get_optimal_shor(addition_cost, n, epsilon = 0.001):
	#addition_cost = point_addition_cost(n)
	eight_lookup_cost = Lookup_Cost(n, 8).multiply(6)
	# The cost of an addition without any lookups
	# We also want to remove the qubits, too
	blank_addition_cost = addition_cost.subtract(eight_lookup_cost)
	blank_addition_cost.low_depth.width -= eight_lookup_cost.low_depth.width
	blank_addition_cost.low_T.width -= eight_lookup_cost.low_T.width
	blank_addition_cost.low_width.width -= eight_lookup_cost.low_width.width
	best_T_size = 8
	best_T = addition_cost.multiply(2*n)
	best_depth = addition_cost.multiply(2*n)
	best_depth_size = 8
	best_width = addition_cost.multiply(2*n)
	best_width_size = 8
	# Check all window sizes up to n/2
	for i in range(int(n/2)):
		# number of windows
		num_windows = math.floor( n / (i+1))
		# size of remainder window
		remainder_window = max(n - num_windows * (i+1), 0)
		main_lookup_costs = Lookup_Cost(n,i).multiply(6)
		# Add in the cost of doing that many lookups
		main_addition_cost = blank_addition_cost.add(main_lookup_costs)
		# The number of point additions that need to be done
		total_cost = main_addition_cost.multiply(2*num_windows)
		
		# If there is a "remainder window" (a window smaller than the 
		# others to finish the remaining bits), add that cost
		if remainder_window > 0:
			second_lookup_costs = Lookup_Cost(n,remainder_window).multiply(6)
			second_addition_cost = blank_addition_cost.add(second_lookup_costs)
			total_cost = total_cost.add(second_addition_cost.multiply(2))
			#Here we add whichever width is greater
			total_cost.low_depth.width += max(second_lookup_costs.low_depth.width, main_lookup_costs.low_depth.width)
			total_cost.low_T.width += max(second_lookup_costs.low_T.width, main_lookup_costs.low_T.width)
			total_cost.low_width.width += max(second_lookup_costs.low_width.width, main_lookup_costs.low_width.width)
		else:
			#With no remainder window, we add just the main lookup width
			total_cost.low_depth.width += main_lookup_costs.low_depth.width
			total_cost.low_T.width += main_lookup_costs.low_T.width
			total_cost.low_width.width += main_lookup_costs.low_width.width
		# Compare to previous best counts, update as needed
		# Choose by lowest T count
		if total_cost.low_T.T_count < best_T.low_T.T_count:	
			best_T = copy.deepcopy(total_cost)
			best_T_size = i
		# Choose by lowest physical width
		if total_cost.low_width.surface_code_costs(n,epsilon)['physical_qubits'] < best_width.low_width.surface_code_costs(n)['physical_qubits']:
			best_width = copy.deepcopy(total_cost)
			best_width_size = i
		# Choose by lowest depth
		if total_cost.low_depth.full_depth < best_depth.low_depth.full_depth:
			best_depth = copy.deepcopy(total_cost)
			best_depth_size = i

	return {"T": best_T, "depth": best_depth, "width" : best_width, "T-window" : best_T_size, "depth-window" : best_depth_size, "width-window" : best_width_size}


# Physical error rate
physical_err = 0.001

# Check fixed modulus sizes
T_CSV = SingleCost.CSV_header()
depth_CSV = SingleCost.CSV_header()
width_CSV = SingleCost.CSV_header()
# for i in {110,224, 256, 384, 521}:
for i in {256, 384, 521}:
	costs = get_optimal_shor(fixed_modulus_point_addition_cost(i), i, physical_err)
	t_surface = costs['T'].surface_code_costs(i, physical_err)
	depth_surface = costs["depth"].surface_code_costs(i, physical_err)
	width_surface = costs["width"].surface_code_costs(i, physical_err)
	T_CSV += str(i) + ", " + costs["T"].low_T.csv_row() + "," + str(costs["T-window"]) + "," + surface_code_message(t_surface['T']) + "\n"
	depth_CSV += str(i) + ", " + costs["depth"].low_depth.csv_row() + "," + str(costs["depth-window"]) + "," + surface_code_message(depth_surface['Depth']) + "\n"
	width_CSV += str(i) + ", " + costs["width"].low_width.csv_row() + "," + str(costs["width-window"]) + "," + surface_code_message(width_surface['Width']) + "\n"


t_file = open('shor_low_t_fixed.csv', 'a')
t_file.write(T_CSV)
t_file.close()

depth_file = open('shor_low_depth_fixed.csv', 'a')
depth_file.write(depth_CSV)
depth_file.close()

width_file = open('shor_low_width_fixed.csv', 'a')
width_file.write(width_CSV)
width_file.close()
