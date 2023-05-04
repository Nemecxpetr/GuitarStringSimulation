"""
MPC-FYA 2023 - Individual Project

Author: Bc. Petr Nemec
Student ID: 221480
Email: xnemec93@vutbr.cz

Assignment: The creation of wave function with transients on string 

Script: 

General one dimensional wave function:
d2y/dt^2 = c^2 d2y/dx^2 - b * dy/dt

"""
import time
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

"""
# Set the initial conditions
"""

# Set up the simulation parameters
dt = 0.0002 # time step
fps = 25 # frames per second
animation_length = 10**4
interval = 1./fps

# Initial string length in cm
string_length = 40

T = 80 # N tension on ends
mu = 6.8 # longitude density

# initial string speed
c = math.sqrt(T/mu)*100

# damping factor
damp = dt*6*10**-5 # set to 0 only viscosity dampening for now
# EDIT: viscosity ensures that now 
# TO DO: change damping factor to a function dependent on frequency of the wave function (higher frequencies will be dumped more
viscosity = 5000

class String():    
    def __init__(self, x, y0, c):
        # Initialize the string's state
        self.y = self.pad_array(y0) 
        self.x = self.pad_array(x)
        self.y0 = self.pad_array(y0)
        self.c = c
        self.viscosity = viscosity
        self.y_prev = np.copy(self.y0)

    # Pad the array with ghost points to avoid boundary conditions
    def pad_array(self, arr):
        padded_arr = np.concatenate(([arr[0]] * (len(arr)-1), arr, [arr[-1]] * (len(arr)-1)))
        return padded_arr[len(arr)-1 : -len(arr)+1]

    def increment(self, dt):
        """Increment shape of string by dt with viscosity damping"""
        # Calculate the time step based on the string speed and grid spacing
        r = (self.c * dt / np.gradient(self.x))**2

        # Update the string shape using the wave equation with viscosity damping
        temp = np.copy(self.y)
        # first delete old line (2y-yold) and then increment with diff. eq 
        # (here done by multiplying r with y - indexes set the initial condition of being mounted at both ends)
        self.y[1:-1] = 2 * self.y[1:-1] - self.y_prev[1:-1] + r[1:-1] * (self.y[2:] - 2 * self.y[1:-1] + self.y[:-2]) 
        # then add viscosity 
        # If we comment out the line before we will so what it does to the initial state alone - it flatens it out 
        # when in motion it basically damps the higher freqencies the higher the viscosity the bigger damping
        self.y[1:-1] += self.viscosity * dt * r[1:-1] * (self.y[2:] + self.y[:-2] - 2 * self.y[1:-1])
        # amplitude damping done from diff. eq.
        self.y[1:-1] -= damp * (self.y[1:-1] - self.y_prev[1:-1]) / dt

        # Update the previous string shape for the next time step
        self.y_prev = temp
"""
    def increment(self, dt):
        """"""Increment shape of string by dt with frequency-dependent damping""""""
        # Calculate the time step based on the string speed and grid spacing
        r = (self.c * dt / np.gradient(self.x))**2

        # Calculate the frequency at each point on the string
        f = self.c / (2 * self.x[-1] * np.sqrt(r))

        # Update the string shape using the wave equation and frequency-dependent damping factor
        temp = np.copy(self.y)
        self.y[1:-1] = 2 * self.y[1:-1] - self.y_prev[1:-1] + r[1:-1] * (self.y[2:] - 2 * self.y[1:-1] + self.y[:-2])
        self.y[1:-1] -= dt * self.freq_damping(f[1:-1]) * (self.y[1:-1] - self.y_prev[1:-1]) / dt

        # Update the previous string shape for the next time step
        self.y_prev = temp
"""

"""
Initial wave shape 
"""
# create empty matrix to be filled with shape
x = np.linspace(1, string_length, 2**8)
y = np.empty_like(x)

# x is percentage length coordinate on string (so it should be between 0 and 1)
# y is y coordinate of deflection in cm
#v[x, y]
#v[0, 1]
v0 = [.3, 2]
v1 = [.3, 1]
v2 = [.3, 2]
v3 = [.4, 5]
v4 = [.5, 1]

v = [[0, 0], v0, v2, v3, v4, [1, 0]]

"""
# Shapes:
"""


"""
# 'CONECTED DOTS'
# Sort the points v0 and v1 based on their x-coordinates
#if v0[0] > v1[0]:
#    v0, v1 = v1, v0

# x[-1] refers to the last point in x 
# - this sets the lenght of the string,
#  hence v is just percentage of lenght
#condition = np.where(x <= v[0][0]*x[-1])
#y[condition] = v[1][1] / (v[1][0] * x[-1]) * x[condition]
"""


""" CONNECTED VERTEXES
xold, yold = x[0], y[0]

for [xnew, ynew] in v:

    condition = np.logical_and(x > xold*x[-1], x <= xnew*x[-1])
    # Compute the slope and intercept of the line connecting the middle points
    slope = (ynew - yold) / ((xnew - xold)*x[-1])
    intercept = yold - slope * xold
    
    y[condition] = slope * (x[condition] - xold*x[-1]) + intercept
    
    xold, yold = xnew, ynew
"""

"""
#condition = np.where(x > xold*x[-1])
#y[condition] = -v[-1][1] / ((1.0 - v[-1][0]) * x[-1]) * (x[condition] - x[-1])
"""
# 'GAUSS'
sigma = .4 # d0 * (x[-1] - x[0]) # quality factor the smaller the number the steeper the G. function
y = v2[1] * np.exp(-(x - v2[0]*x[-1])**2 / (2*sigma**2)) # Gaussian function shape
# can be used with
# 'DOUBLE GAUSS'
#sigma = 3
#y[np.where(x > v3[0]*x[-1])] = v4[1] * np.exp(-(x[np.where(x > v3[0]*x[-1])] - v4[0]*x[-1])**2 / (2*sigma**2)) # Gaussian function shape

# 'SINE'
y = v0[1]*np.sin(2*np.pi/(x[-1] - x[0]) * x)

# 'NOISE'
y[1:-1] = np.random.uniform(-2,2,254)
# OR
#y[1:-1] = np.random.randint(-2,2,254)

# 'SAW'
y[np.where(x <= v0[0]*x[-1])] = v0[1]/(v0[0]*x[-1]) * x[np.where(x <= v0[0]*x[-1])]
y[np.where(x > v0[0]*x[-1])] = -v0[1]/((1.0-v0[0])*x[-1]) * (x[np.where(x > v0[0]*x[-1])] - x[-1])


# Initialize the string object with the chosen initial shape
string = String(x, y, c)

class MyAnimation:
    def __init__(self):
        """
        Plotting
        """
        global line, string, fig, ax
        # Create ínstance of the string object
        self.string = string

        # Set up the plot
        #fig, ax = plt.subplots()
        line, = ax.plot(self.string.x, self.string.y)
        ax.set_ylim(-1.1*v0[1], 1.1*v0[1])
        #ax.set_xlim()
        ax.set_title("String wave function animation")
        ax.set_xlabel("String lenght (cm)")
        ax.set_ylabel("String Deflection (cm)")

        self.paused = False

        self.animation = FuncAnimation(fig, self.update, init_func=self.init, frames=range(animation_length), blit=True, interval=interval, repeat=True)
        self.animation.pause()

        # Connect event handlers
        fig.canvas.mpl_connect('button_press_event', self.toggle_pause)

    def toggle_pause(self, *args, **kwargs):
        if self.paused:
            self.animation.resume()
        else:
            self.animation.pause()
        self.paused = not self.paused

    def init(self):
        global line
        return line,

    def update(self, frame_no):
        global string
        # Real-time plotting parameters
        t0 = time.time()
        t = t0
        t_next = time.time() + interval
        # while the time of the simulation has started
        while t_next <= t0:
             # Increment the string shape by the time step and damping factor
             string.increment(dt)
             t+=dt
        # Increment the string shape for the current frame with a damping factor
        string.increment(dt)
        line.set_ydata(string.y)
        return line,

# Set up the plot
fig, ax = plt.subplots()
line, = ax.plot(string.x, string.y)
ax.set_title("Initial string deflexion shape")
ax.set_xlabel("String lenght (cm)")
ax.set_ylabel("String Deflection (cm)")

pa = MyAnimation()
    
#ani = MyAnimation()

plt.show() 
