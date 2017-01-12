### Create a dictionary of physical constants

pcs = {}

# Seconds per day
pcs['spd'] = 60.0 * 60.0 * 24.0
# Seconds in a year
pcs['spy'] = 60.0 * 60.0 * 24.0 * 365.0              
# Seconds per month 
pcs['spm'] = pcs['spy'] / 12.0 
# Density of water (kg / m^3)
pcs['rho_w'] = 1000.0  
# Density of ice (kg / m^3)
pcs['rho_i'] = 910.0
# Gravitational acceleration (m / s^2)
pcs['g'] = 9.8 
# Flow rate factor of ice (1 / Pa^3 * s) 
pcs['A'] = 2.5e-25
# Average bump height (m)
pcs['h_r'] = 0.1
# Typical spacing between bumps (m)
pcs['l_r'] = 2.0      
# Sheet conductivity (m^(7/4) / kg^(1/2))
pcs['k'] = 5e-3
# Exponents 
pcs['alpha'] = 5. / 4.
pcs['beta'] = 3. / 2.
pcs['delta'] = pcs['beta'] - 2.0