import numpy as np
from scipy.integrate import odeint
import math
from matplotlib import pyplot as plt
import matplotlib as mpl

TIME_INTERVAL = 30.0
STEPS = 1000000.0
TIME_PER_FRAME = .001

#meters
init_r = 600001.0
r = init_r
#meters/sec
init_r_vel = 0.0
r_vel = init_r_vel
#meters/sec^2
init_r_acc = 0.0
r_acc = init_r_acc
#seconds
init_t = 0.0
t = init_t

#meters
init_theta_x = 0.0
theta_x = init_theta_x
init_theta_y = 0.0
theta_y = init_theta_y

#meters/sec
init_theta_x_vel = 0.0
theta_x_vel = init_theta_x_vel
init_theta_y_vel = 0.0
theta_y_vel = init_theta_y_vel

#meters/sec^2
init_theta_x_acc = 0.0
theta_x_acc = init_theta_x_acc
init_theta_y_acc = 0.0
theta_y_acc = init_theta_y_acc

#kgs
init_mass = 6690

mass = init_mass

init_fuel_mass = 4000

#mass of Kerbin (kg) 
k = 5.2915158e+22

#gravitational constant (m^2/s^2)
G = 3.5316000e+12

#atmosphere scaling for Kerbin (m)
H = 5600.0

#pressure in pascals (Pa)(J/m^3)(kg/(m*s))(N/m^2)
p_init = 101325.0
p = p_init

#distance of sea level from the center of the planet (radius)
#m
s = 600000.0

#cross-section (m^2) Default is the flat surface of a small part (1.25m in radius)
#using 8.0% or slightly less of the flat area cross-section seems to produce accurate results
A = .4

#specific gas constant (J/kg*K)
R = 287.053

#temp (K) average temp in atmosphere
T = 240.0

#dimensionless
drag_coeff = .15

#specific impulse (N/(kg/s))
isp = 265.0

def mass_func(t):
    sol_space = np.linspace(0.0, t, 1000)
    y0 = [0.0]
    #print(sol_space)
    #print(t)
    y_r = odeint(thrust_func_r, y0, sol_space)
    y_theta_x = odeint(thrust_func_theta_x, y0, sol_space)
    #print(len(y))
    #print(init_mass - (y[len(y)-1])*isp)
    if (y_r[len(y_r)-1] + y_theta_x[len(y_theta_x)-1])*isp < init_fuel_mass:
        return init_mass - (y_r[len(y_r)-1] + y_theta_x[len(y_theta_x)-1])*isp
    else:
        return init_mass - init_fuel_mass

#thrust (N)
def thrust_func_r(t):
    if t < 20.0:
        return 197897.0
    #elif t < 40.0:
    #    return 197897.0/2.0
    else:
        return 0.0
    
def thrust_func_theta_x(t):
    if t >= 20.0 and t < 24.0:
        return 197897.0
    else:
        return 0.0

    
#data about r, time, then init time and step length
def accel_eq_r():
    global t, k, G, r, p, R, T, H, s, r_acc, mass
    new_r_acc = (1.0/(mass))*(((-1.0*mass*G)/(r**(2.0)))+thrust_func_r(t)-.5*drag_coeff*(p/(R*T))*A*(r_vel)**2)
    #print("drag: " + str(-1.0*.5*drag_coeff*(p/(R*T))*A*(r_vel)**2))
    #print("gravity: " + str(-1.0*mass_func(t)*G/(r**(2.0))))
    return new_r_acc

def accel_eq_theta_x():
    global t, k, G, r, p, R, T, H, s, theta_x_acc, mass
    new_theta_x_acc = (1.0/(mass))*(thrust_func_theta_x(t)-.5*drag_coeff*(p/(R*T))*A*(theta_x_vel)**2)
    #print("drag: " + str(-1.0*.5*drag_coeff*(p/(R*T))*A*(r_vel)**2))
    #print("gravity: " + str(-1.0*mass_func(t)*G/(r**(2.0))))
    return new_theta_x_acc


def update_approximate(interval):
    global r_acc, r_vel, r, t, mass, p, theta_x, theta_x_vel, theta_x_acc
    mass = mass - interval*(thrust_func_r(t)+thrust_func_theta_x(t))/(isp*9.80665)
    p = p_init*math.exp((1/H)*(r-s)*-1.0)
    #print(mass)
    if mass < init_mass-init_fuel_mass:
        #exit()
        mass = init_mass-init_fuel_mass
    if r > s:
        r_acc = accel_eq_r()
        r_vel = r_vel + r_acc*interval
        r = r + r_vel*interval
        theta_x_acc = accel_eq_theta_x()
        theta_x_vel = theta_x_vel + theta_x_acc*interval
        theta_x = theta_x + theta_x_vel*interval
        #print(theta_x_vel)
    else:
        r = s
        r_acc = 0
        r_vel = 0
        theta_x_acc = 0
        theta_x_vel = 0

def calc_delta_v():
    delta_v = np.log(init_mass/(init_mass-init_fuel_mass))*isp*9.81
    return delta_v
        
x_col = []
z_col = []
max_r_vel = 0.0
max_theta_x_vel = 0.0
for x in xrange(0, int(STEPS)):
    #print(mass_func(t))
    #print("height: " + str(r) + " vel: " + str(r_vel) + " acc: " + str(r_acc))
    if x%(int((STEPS/TIME_INTERVAL)*(TIME_PER_FRAME))) == 0:
        z_col.append(r*np.sin(theta_x))
        x_col.append(r*np.cos(theta_x))
    update_approximate(TIME_INTERVAL/STEPS)
    if r_vel > max_r_vel:
        max_r_vel = r_vel
    if theta_x_vel > max_theta_x_vel:
        max_theta_x_vel = theta_x_vel
    #print(mass_func(t))
    t = t+(TIME_INTERVAL/STEPS)
    
print("delta v: " + str(calc_delta_v()))
print("final radius: " + str(r))
print("final x coord: " + str(x_col[len(x_col)-1]))
print("final z coord: " + str(z_col[len(z_col)-1]))
print("max radial velocity: " + str(max_r_vel))
print("max angular velocity: " + str(max_theta_x_vel))


#print(mpl.rcParams['agg.path.chunksize'])
mpl.rcParams['agg.path.chunksize'] = STEPS*TIME_PER_FRAME

ax = plt.subplot(111)
ax.plot(x_col, z_col)
ax.grid(True)
plt.show()
