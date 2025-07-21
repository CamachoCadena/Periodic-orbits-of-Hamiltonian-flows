import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
############ Trace of figure 8 in paracolic locus #############
# function in coordinates
def traceFig8Par(s1, t1):
    numerator = (1 + t1) ** 2 + 3 * s1 ** 2 * (1 + t1) ** 2 + s1 ** 3 * (1 + t1) ** 2 + 3 * s1 * (
                    1 + 3 * t1 + t1 ** 2)
    denominator = s1 * t1
    return numerator / denominator

# ODE in coordinates
def ODETraceFig8Par(y, t):
    s1, t1 = y  # Extract the current values of s1 and t1

    # Define the system of differential equations
    ds1dt = (2 * (1 + s1) ** 3 * (-1 + t1 ** 2)) / (s1 * t1)
    dt1dt = - (2 * (-2 + s1) * (1 + s1) ** 2 * (1 + t1) ** 2) / s1 ** 2

    return [ds1dt, dt1dt]
# Plot level sets
x = np.linspace(0.05, 8, 1000)  # Avoid x=0 and y=0 to prevent division by zero
y = np.linspace(0.05, 50, 1000)
X, Y = np.meshgrid(x, y)
# Calculate Z values based on f(x, y)
Z = traceFig8Par(X, Y)

# Plot the level set using contour
levels = [35,40,50,60,70,80,90,100,150,200,250,300,350]
plt.figure(figsize=(6, 6))
cp = plt.contour(X, Y, Z, levels=levels, cmap='viridis')
plt.colorbar(cp)
plt.xlabel(r'$\sigma_1$')
plt.ylabel(r'$\tau_1$')
plt.savefig('LevelSetsParabolicFig8.eps', format='eps')
plt.show()
# Plot numerical solution
initial_conditions = [1, 1]  # Fuchsian point in parabolic locus
t = np.linspace(0, 0.5, 100000)
solution = odeint(ODETraceFig8Par, initial_conditions, t)
s1_solution = solution[:, 0]
t1_solution = solution[:, 1]
plt.figure(figsize=(10, 6))
plt.plot(t, s1_solution, label=r'$\sigma_1(t)$', color='teal', linestyle='-', linewidth=2)
plt.plot(t, t1_solution, label=r'$\tau_1(t)$', color='violet', linestyle='-', linewidth=2)
plt.xlabel('Time (t)')
plt.ylabel(r'Hamiltonian flow, $\sigma_1$ and $\tau_1$')
plt.legend()
plt.savefig('solutionHamiltonianFig8ParFuchs.eps', format='eps')
plt.tight_layout()
plt.show()

############## Trace of the figure 8 in symplectic leaf determined by F = (3,1/3,6/1,6,8,1/8) ##########
# function in coordinates
def traceFig8Fuchs(t1, s1):
    numerator = (4 + 580 * t1 + 576 * t1**2 +
                 6 * s1**3 * (1 + t1)**2 +
                 s1**2 * (17 + 107 * t1 + 90 * t1**2) +
                 s1 * (15 + 949 * t1 + 408 * t1**2))
    denominator = 24 * s1 * t1
    return numerator / denominator

# ODE in coordinates
def ODEFig8Fuchs(y, t):
    s1, t1 = y  # Extract the current values of s1 and t1

    # Define the system of differential equations
    ds1dt =(4 - 576 * t1**2 + s1 * (15 - 408 * t1**2) + s1**2 * (17 - 90 * t1**2) - 6 * s1**3 * (-1 + t1**2)) / (12 * t1)
    dt1dt = (1 + t1) * (-4 - 576 * t1 + 12 * s1**3 * (1 + t1) + s1**2 * (17 + 90 * t1)) / (12 * s1)

    return [ds1dt, dt1dt]
# Plot level sets
x = np.linspace(0.05, 10, 1000)  # Avoid x=0 and y=0 to prevent division by zero
y = np.linspace(0.05, 16, 1000)
X, Y = np.meshgrid(x, y)
# Calculate Z values based on f(x, y)
Z = traceFig8Fuchs(X, Y)

# Plot the level set using contour
levels = [35,40,50,60,70,80,90,100,150,200,210,230,240,250,300,310,350,400,420]
plt.figure(figsize=(6, 6))
cp = plt.contour(X, Y, Z, levels=levels, cmap='viridis')
plt.colorbar(cp)
plt.xlabel(r'$\sigma_1$')
plt.ylabel(r'$\tau_1$')
plt.savefig('LevelSetsParabolicFig368.eps', format='eps')
plt.show()
# Plot numerical solution for initial condition the Fuchsian point
initial_conditions = [2, 1]  # Fuchsian point in the symplectic leaf determined by F = (3,1/3,6/1,6,8,1/8)
t = np.linspace(0, 0.5, 100000)
solution = odeint(ODEFig8Fuchs, initial_conditions, t)
s1_solution = solution[:, 0]
t1_solution = solution[:, 1]
plt.figure(figsize=(10, 6))
plt.plot(t, s1_solution, label=r'$\sigma_1(t)$', color='teal', linestyle='-', linewidth=2)
plt.plot(t, t1_solution, label=r'$\tau_1(t)$', color='violet', linestyle='-', linewidth=2)
plt.xlabel('Time (t)')
plt.ylabel(r'Hamiltonian flow, $\sigma_1$ and $\tau_1$')
plt.legend()
plt.savefig('solutionHamiltonianFig8368Fuchs.eps', format='eps')
plt.tight_layout()
plt.show()
# Plot numerical solution for initial condition (4,3)
initial_conditions = [4, 3]
t = np.linspace(0, 0.5, 100000)
solution = odeint(ODEFig8Fuchs, initial_conditions, t)
s1_solution = solution[:, 0]
t1_solution = solution[:, 1]
plt.figure(figsize=(10, 6))
plt.plot(t, s1_solution, label=r'$\sigma_1(t)$', color='teal', linestyle='-', linewidth=2)
plt.plot(t, t1_solution, label=r'$\tau_1(t)$', color='violet', linestyle='-', linewidth=2)
plt.xlabel('Time (t)')
plt.ylabel(r'Hamiltonian flow, $\sigma_1$ and $\tau_1$')
plt.legend()
plt.savefig('solutionHamiltonianFig8368IC43.eps', format='eps')
plt.tight_layout()
plt.show()

################ Trace of the commutator in parabolic locus ##########
# function in coordinates
def traceCommPar(s1, t1):
    numerator =  ((1 + t1)**3 + s1**6 * (1 + t1)**3 +
                           3 * s1 * (1 + t1)**2 * (1 + 2 * t1) +
                           3 * s1**5 * (1 + t1)**2 * (1 + 2 * t1) +
                           3 * s1**2 * (1 + t1)**2 * (1 + 5 * t1) +
                           3 * s1**4 * (1 + t1)**2 * (1 + 5 * t1) +
                           s1**3 * (2 + 27 * t1 + 42 * t1**2 + 20 * t1**3))
    denominator = (s1**3 * t1)
    return numerator / denominator
# ODE in coordinates
def ODEcommPar(y, t):
    s1, t1 = y  # Extract the current values of s1 and t1

    # Define the system of differential equations
    ds1dt = -((2 * (1 + s1)**4 * (1 + t1) * (-1 + t1 + 2 * t1**2 + s1**2 * (-1 + t1 + 2 * t1**2) +
                                          s1 * (1 - t1 + 4 * t1**2))) / (s1**2 * t1))
    dt1dt = (6 * (-1 + s1) * (1 + s1)**3 * (1 + t1)**2 * (1 + t1 + 2 * s1 * t1 +
                                                     s1**2 * (1 + t1))) / s1**3

    return [ds1dt, dt1dt]

# Plot level sets
x = np.linspace(0.05, 3, 1000)  # Avoid x=0 and y=0 to prevent division by zero
y = np.linspace(0.05, 1.9, 1000)
X, Y = np.meshgrid(x, y)
# Calculate Z values based on f(x, y)
Z = traceCommPar(X, Y)

# Plot the level set using contour
levels = [35,40,50,60,70,80,90,100,150,200,210,230,240,250,300,310,350,400,420,450,470,490,500,520,550,600]
plt.figure(figsize=(6, 6))
cp = plt.contour(X, Y, Z, levels=levels, cmap='viridis')
plt.colorbar(cp)
plt.xlabel(r'$\sigma_1$')
plt.ylabel(r'$\tau_1$')
plt.savefig('LevelSetsParabolicComm.eps', format='eps')
plt.show()
# Plot numerical solution for initial condition the Fuchsian point
initial_conditions = [1, 1]  # Fuchsian point in parabolic locus
t = np.linspace(0, 0.5, 100000)
solution = odeint(ODEcommPar, initial_conditions, t)
s1_solution = solution[:, 0]
t1_solution = solution[:, 1]
plt.figure(figsize=(10, 6))
plt.plot(t, s1_solution, label=r'$\sigma_1(t)$', color='teal', linestyle='-', linewidth=2)
plt.plot(t, t1_solution, label=r'$\tau_1(t)$', color='violet', linestyle='-', linewidth=2)
plt.xlabel('Time (t)')
plt.ylabel(r'Hamiltonian flow, $\sigma_1$ and $\tau_1$')
plt.legend()
plt.savefig('solutionHamiltonianParComm.eps', format='eps')
plt.tight_layout()
plt.show()
####################### Trace of a^kc^(-1) in parabolic locus #################
# function in coordinates
def tracePowPar(s1, t1, N):
    numerator =  (6 * s1 * t1 + N * (1 + s1)**3 * (1 + t1)**2 +
        N**2 * (1 + s1)**3 * (1 + t1)**2)
    denominator = (2 * s1 * t1)
    return numerator / denominator
#ODE in coordinates
def ODEPowPar(y, t):
    s1, t1 = y  # Extract the current values of s1 and t1

    # Define the system of differential equations
    ds1dt = -((3 * (1 + 3) * (1 + s1)**3 * (-1 + t1**2)) / t1)
    dt1dt =  (3 * (1 + 3) * (1 + s1)**2 * (-1 + 2 * s1) * (1 + t1)**2) / s1

    return [ds1dt, dt1dt]
# Plot level sets
x = np.linspace(0.05, 2.5, 1000)  # Avoid x=0 and y=0 to prevent division by zero
y = np.linspace(0.05, 8, 1000)
X, Y = np.meshgrid(x, y)
# Calculate Z values based on f(x, y)
Z = tracePowPar(X, Y,3)

# Plot the level set using contour
levels = [165,170,180,200,210,220,230,240,250,270,280,300,310,320,330,350,400]
plt.figure(figsize=(6, 6))
cp = plt.contour(X, Y, Z, levels=levels, cmap='viridis')
plt.colorbar(cp)
plt.xlabel(r'$\sigma_1$')
plt.ylabel(r'$\tau_1$')
plt.savefig('LevelSetsParabolicPowers.eps', format='eps')
plt.show()
# Plot numerical solution for initial condition the Fuchsian point
initial_conditions = [1, 1]  # Fuchsian point in parabolic locus
t = np.linspace(0, 0.5, 100000)
solution = odeint(ODEPowPar, initial_conditions, t)
s1_solution = solution[:, 0]
t1_solution = solution[:, 1]
plt.figure(figsize=(10, 6))
plt.plot(t, s1_solution, label=r'$\sigma_1(t)$', color='teal', linestyle='-', linewidth=2)
plt.plot(t, t1_solution, label=r'$\tau_1(t)$', color='violet', linestyle='-', linewidth=2)
plt.xlabel('Time (t)')
plt.ylabel(r'Hamiltonian flow, $\sigma_1$ and $\tau_1$')
plt.legend()
plt.savefig('solutionHamiltonianParPowers.eps', format='eps')
plt.tight_layout()
plt.show()