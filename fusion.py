

#input settings
n0 = 4.6e13                 #electron background density in #/m^3
phi_Cathode = -100000        #cathode potential
phi0 = 0                    #reference potential
Te = np.abs(phi_Cathode)    #electron temperature in eV
Ti = 0.1                    #ion velocity in eV (not used yet)
vth = np.sqrt(2*QE*Ti/m_ion)   #thermal velocity with Ti in eV
Operating_Pressure   = 7      # Pa (Not used yet)

# Calculate specific weight and prepair to insert particles
np_insert = 400                      #insert 2 particles per anode cell.
flux = 4.6e20                          #Flux of entering ions [ions per second]
npt = flux*dt
spwt = npt/np_insert                    #specific weight, real particles per macroparticle
mp_q = 1                                #macroparticle charge
max_part = 200000                       #buffer size
#allocate particle array
part_x = np.zeros([max_part,2]) #particle positions
part_v = np.zeros([max_part,2]) #particle velocities
sourse_storage = np.zeros([max_part,2]) #particle sourses

#calculate maximum expected velocity and timestep
E_av = (np.abs(phi_Cathode) - 0) / (R2 - R1)
a_av = E_av*QE / m_ion
v_max = np.sqrt(2*a_av*R2)
dt = 1e-10     #time step size, at vmax move 0.10dx

# Define the potential solver function

PHI_B = np.zeros([nx, ny])
idx = np.round(nx/2).astype(int)
idy = np.round(ny/2).astype(int)

#INITIALIZE
# PHI_B[INDEX_X_Cathode + idx, INDEX_Y_Cathode + idy] = -20
# PHI_B[INDEX_X_Anode + idx, INDEX_Y_Anode + idy] = 0
# ni=0
#n0=0
# EPS0=1
# kb = 1
# iterz = 300
#den = 0*PHI_B
def get_Potential(PHI_B, den, nx, ny, iters):
    for k in range(iters):
        PHI_OLD = PHI_B
        for i in range(1,nx-2):
            for j in range(1, ny-2):
                ni = den[i,j]
                rho = QE*(ni - n0*np.exp((PHI_OLD[i,j] - phi0)/(kb*Te))) / EPS0
                chrg = -rho*dx**2
                PHI_B[i,j] = (chrg - PHI_B[i+1, j] - PHI_B[i-1, j] - PHI_B[i, j+1] - PHI_B[i, j-1])/(-4)
        PHI_B[INDEX_X_Cathode + idx, INDEX_Y_Cathode + idy] = phi_Cathode
        PHI_B[INDEX_X_Anode + idx, INDEX_Y_Anode + idy] = phi0
    return PHI_B
#PHI_B = get_Potential(PHI_B, den, nx, ny, iterz)

def sample_Source(nump, R_S):
    xv = np.zeros([nump])
    yv = np.zeros([nump])
    for i in range(nump):
        theta = np.random.rand(1)*2*np.pi   # Generate random polar angle
        x = R_S*np.cos(theta)                # Get x position
        y = R_S*np.sin(theta)                # Get y position
        xv[i] = x
        yv[i] = y
    return np.array([xv, yv])

def fusion_Cross_Section(vx, vy):
    # Takes in velocity componants in m/s and returns a cross section in barns
    E = .5*m_ion*(vx**2 + vy**2)
    E = E*6.242e15  # convert J to KeV
    A1 = 46.097
    A2 = 372
    A3 = 4.36e-4
    A4 = 1.22
    A5 = 0
    AA1 = 47.88
    AA2 = 482
    AA3 = 3.08e-4
    AA4 = 1.177
    AA5 = 0
    term1 = A5 + A2/((A4 - A3*E)**2 + 1)
    term2 = E*(np.exp(A1/np.sqrt(E)) - 1)
    term3 = AA5 + AA2/((AA4 - AA3*E)**2 + 1)
    term4 = E*(np.exp(AA1/np.sqrt(E)) - 1)
    sig1 = term1/term2
    sig2 = term3/term4
    return sig1 + sig2

# print(fusion_Cross_Section(100))
# EVE = np.arange(0, 1000, 1)
# sigma = fusion_Cross_Section(EVE)
# EVEev = EVE
# plt.plot(EVEev, sigma)



# nUm = 85000
# plt.axis('equal')
# vec = sample_Source(nUm, R2)
# plt.scatter(vec[0,:] + R_Chamber, vec[1,:]+ H_Chamber/2)

#############
# MAIN LOOP
############

#INITIALIZE
PHI_M = np.zeros([nx, ny])
idx = np.round(nx/2).astype(int)
idy = np.round(ny/2).astype(int)
PHI_M[INDEX_X_Cathode + idx, INDEX_Y_Cathode + idy] = phi_Cathode
PHI_M[INDEX_X_Anode + idx, INDEX_Y_Anode + idy] = phi0
fuse_pos_x = np.array([])
fuse_pos_y = np.array([])
fuse_time = np.array([])
col_counter = 0
top_counter = 0
bot_counter = 0
left_counter = 0
right_counter = 0
anode_counter = 0
cathode_counter = 0
fuse_counter = 0
num_p = 0                         #Clear number of particles
iters = 600                       #Number of iterations used in the potential solver
print('Beginning Main Loop. This could take a while. \n')

for it in range(ts):
    print('Time Step ', it, 'Particles ', num_p)
    #reset field quantities
    den = np.zeros([nx,ny])              #number density
    efx = np.zeros([nx,ny])              #electric field, x-component
    efy = np.zeros([nx,ny])              #electric field, y-component
    chg = np.zeros([nx,ny])              #charge distribution
    col_counter = 0

    # *** 1. Calculate Charge Density ***
    # deposit charge to nodes
    for p in range(num_p):                       #loop over particles
        fi = (part_x[p, 0] + R_Chamber-dx)/dx    #real i index of particle's cell
        i = np.floor(fi).astype(int)             #integral part
        hx = fi - i                              #the remainder
        fj = (part_x[p,1] + (H_Chamber/2)-dx)/dx #real i index of particle's cell
        j = np.floor(fj).astype(int)             #integral part
        hy = fj - j                              #the remainder
        #interpolate charge to nodes
        chg[i, j] = chg[i, j] + (1-hx)*(1-hy)
        chg[i+1, j] = chg[i+1, j] + hx*(1-hy)
        chg[i, j+1] = chg[i, j+1] + (1-hx)*hy
        chg[i+1, j+1] = chg[i+1, j+1] + hx*hy

    # Calculate the Density
    den = spwt*mp_q*chg / (dx**2)
    den[0,:] = 2*den[0,:]              # Double density since only half volume contributing
    den[nx-1,:] = 2*den[nx-1,:]
    den[:,0] = 2*den[:,0]
    den[:,ny-1] = 2*den[:,ny-1]
    den = den + 1e3                    # Add density floor to help the solver

    # *** 2. Calculate the Potential
    PHI_M = get_Potential(PHI_M, den, nx, ny, iters)
    print('Potential Solution Complete.')
    # *** 3. Calculate the Electric Field ***
    efx[1:nx-2,:] = PHI_M[0:nx-3,:] - PHI_M[2:nx-1,:]   #central difference on internal nodes
    efy[:,1:ny-2] = PHI_M[:,0:ny-3] - PHI_M[:,2:ny-1]   #central difference on internal nodes
    efx[0,:] = 2*(PHI_M[0,:] - PHI_M[1,:])             #forward difference on x=0
    efx[nx-1,:] = 2*(PHI_M[nx-2,:] - PHI_M[nx-1,:])     #backward difference on x=Lx
    efy[:,0] = 2*(PHI_M[:,0] - PHI_M[:,1])              #forward difference on y=0
    efy[:,ny-1] = 2*(PHI_M[:,ny-2] - PHI_M[:,ny-1])     #forward difference on y=Ly
    efx = efx / (dx**2)                                 #divide by denominator
    efy = efy / (dx**2)

    # *** 4. Generate New Particles
    print('Generating Particles')
    #insert particles randomly distributed in y and in the first cell
    Posv = sample_Source(np_insert, R_Sourse)
    part_x[num_p:num_p+np_insert, 0] = Posv[0,:]   #x position
    part_x[num_p:num_p+np_insert, 1] = Posv[1,:]   #y position
    #sourse_storage[num_p:num_p+np_insert, 0] = Posv[0,:]
    #sourse_storage[num_p:num_p+np_insert, 1] = Posv[1,:]
    #sample maxwellian in x and y
    pt1 = np.random.rand(np_insert)
    pt2 = np.random.rand(np_insert)
    pt3 = np.random.rand(np_insert)
    pt11 = np.random.rand(np_insert)
    pt12 = np.random.rand(np_insert)
    pt13 = np.random.rand(np_insert)
    part_v[num_p:num_p+np_insert,0] = (-1.5 + pt1 + pt2 + pt3)*vth     # x velocity
    part_v[num_p:num_p+np_insert,1] = (-1.5 + pt11 + pt12 + pt13)*vth  # y velcoity
    num_p = num_p + np_insert                                            #increment particle counter
    # *** Move Particles ***
    print('Moving Particles...   DO NOT INTERRUPT')
    p=0
    while p < num_p:        # Loop over particles
        fi = (part_x[p, 0] + R_Chamber-dx)/dx   # i index of particle's cell. Taking into account
        i = np.floor(fi).astype(int)         # that the physical domain is centered at (0,0)
        hx = fi - i
        fj = (part_x[p,1] + (H_Chamber/2)-dx)/dx
        j = np.floor(fj).astype(int)
        hy = fj-j
        #print('Fi: ', fi, 'Fj: ', fj)
        # gather electric field
        E = np.array([0,0])
        E = np.array([efx[i,j], efy[i,j]])*(1-hx)*(1-hy)      #contribution from (i,j)
        E = E + np.array([efx[i+1,j], efy[i+1,j]])*hx*(1-hy)            #(i+1,j)
        E = E + np.array([efx[i,j+1], efy[i,j+1]])*(1-hx)*hy            #(i,j+1)
        E = E + np.array([efx[i+1,j+1], efy[i+1,j+1]])*hx*hy            #(i+1,j+1)
        #print('E:', E)
        # Update Velocity and Position
        F = QE*E      # Lorenz force F = qE
        a = F/m_ion   # Acceleration
        part_v[p,:] = part_v[p,:] + a*dt

        part_x[p,:] = part_x[p,:] + part_v[p,:]*dt
        #print(part_v[p,:])
        #print(part_x[p,:])
        # Get Fusion probability
        vx = part_v[p, 0]
        vy = part_v[p, 1]
        delx = vx * dt
        dely = vy * dt
        path_len = np.sqrt(delx**2 + dely**2)
        sigma = fusion_Cross_Section(vx, vy)*1e-28 # microscopic cross section [m^2]

        # gather density at particle position
        Rho = den[i,j]*(1-hx)*(1-hy)      #contribution from (i,j)
        Rho = Rho + den[i+1,j]*hx*(1-hy)            #(i+1,j)
        Rho = Rho + den[i,j+1]*(1-hx)*hy            #(i,j+1)
        Rho = Rho + den[i+1,j+1]*hx*hy            #(i+1,j+1)

        # Calculate macroscopic cross section
        SIGMA = Rho * sigma
        Prob_fusion = SIGMA * path_len * spwt

        # Prepare to check for fusion
        random = np.random.rand(1)

        # Process Boundries and Sinks
        R = np.sqrt(part_x[p,0]**2 + part_x[p,1]**2)
        #print('Finished Position and Velocity Update')
        #print('Checking Top Wall')
        # Top Wall
        if part_x[p,1] > H_Chamber/2:
            part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
            part_v[p,:] = part_v[num_p-1,:]
            part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
            part_v[num_p-1,:] = np.array([0,0])
            num_p = num_p - 1              # Reduce particle count
            p = p - 1                      # Reduce particle index
            col_counter = col_counter + 1
            top_counter = top_counter + 1

        # Bottom Wall
        elif part_x[p,1] < -H_Chamber/2:
            #print('Checking Bottom Wall')
            part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
            part_v[p,:] = part_v[num_p-1,:]
            part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
            part_v[num_p-1,:] = np.array([0,0])
            num_p = num_p - 1              # Reduce particle count
            p = p - 1                      # Reduce particle index
            col_counter = col_counter + 1
            bot_counter = bot_counter + 1

        # Right Wall
        elif part_x[p,0] > R_Chamber:
            #print('Checking Right Wall')
            part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
            part_v[p,:] = part_v[num_p-1,:]
            part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
            part_v[num_p-1,:] = np.array([0,0])
            num_p = num_p - 1              # Reduce particle count
            p = p - 1                      # Reduce particle index
            col_counter = col_counter + 1
            right_counter = right_counter + 1

        # Left Wall
        elif part_x[p,0] < -R_Chamber:
            #print('Checking Left Wall')
            part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
            part_v[p,:] = part_v[num_p-1,:]
            part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
            part_v[num_p-1,:] = np.array([0,0])
            num_p = num_p - 1              # Reduce particle count
            p = p - 1                      # Reduce particle index
            col_counter = col_counter + 1
            left_counter = left_counter + 1

        # Process grids

        # Anode
        elif (R < R2 + r_wire) and (R > R2 - r_wire):
            #print('Check Anode')
            prob = np.random.rand(1)
            if prob > GEO_A:
                # Delete particle if it collides with the grid
                part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
                part_v[p,:] = part_v[num_p-1,:]
                part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
                part_v[num_p-1,:] = np.array([0,0])
                num_p = num_p - 1              # Reduce particle count
                p = p - 1                      # Reduce particle index
                col_counter = col_counter + 1
                anode_counter = anode_counter + 1

        # Cathode
        elif (R < R1 + r_wire) and (R > R1 - r_wire):
            #print('Check Cathode')
            prob = np.random.rand(1)
            if prob > GEO_A:
                # Delete particle if it collides with the grid
                part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
                part_v[p,:] = part_v[num_p-1,:]
                part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
                part_v[num_p-1,:] = np.array([0,0])
                num_p = num_p - 1              # Reduce particle count
                p = p - 1                      # Reduce particle index
                col_counter = col_counter + 1
                cathode_counter = cathode_counter + 1

        elif random <= Prob_fusion:
            print('FUSION!\n')
            fuse_pos_x = np.append(fuse_pos_x, part_x[p, 0])
            fuse_pos_y = np.append(fuse_pos_y, part_x[p, 1])
            fuse_time = np.append(fuse_time, dt*it)
            fuse_counter = fuse_counter + 1
            # Delete particle if it fused
            part_x[p,:] = part_x[num_p-1,:]  # Kill particle by replacing with last particle
            part_v[p,:] = part_v[num_p-1,:]
            part_x[num_p-1,:] = np.array([0,0]) # Reset last particle to 0
            part_v[num_p-1,:] = np.array([0,0])
            num_p = num_p - 1              # Reduce particle count
            p = p - 1                      # Reduce particle index

        p = p + 1                              # Move to the next particle
    print('Finished Moving Particles.')
    print('Net Change in Ion Population: ', np_insert - col_counter)
    print(col_counter, ' particles lost.\n')
#     if (it >= 0) and (num_p != 0):
#         clear_output()
#     else:
#         break
