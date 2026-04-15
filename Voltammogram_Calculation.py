import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt  
import pandas as pd

# Define Grünwald–Letnikov fractional derivative function for evenly spaced data. This function returns
# an array with the values of the fractional derivative. It also works for fractional integrals.

def gl_fractional_derivative(x, y, alpha=0.5): # X is the variable and Y is the discrete function to be analyzed. 
                                               # Alpha is the order of the operation, and, if a value is not given 
                                               # when the function is used, a semi-derivative (alpha = 0.5) will be calculated.

    # Turn the X and Y vectors of data into a Numpy array (if they are not already np arrays)
    x = np.array(x)
    y = np.array(y)
    # Calculate the step in the X array
    h = x[1] - x[0]
    # Check if the values in the X array are evenly spaced
    if not np.allclose(np.diff(x), h):
        raise ValueError("X values must be evenly spaced.")

    # Define n as the length of the X array
    n = len(x)
    # Create an n-sized array for the coefficients that will be calculated. They will be placed in the array later
    coeffs = np.zeros(n)
    
    # Binomial coefficient for k=0
    coeffs[0] = 1.0
    # Compute binomial coefficient for order alpha
    for k in range(1, n):
        coeffs[k] = coeffs[k-1] * (alpha - k + 1) / k
    

    # Create an n-sized array for the calculated fractional derivative. They will be placed in the array later
    gl_frac_derv = np.zeros(n)
    # Calculate sum
    for i in range(n):
        s = 0.0
        for k in range(i + 1):
            s += ((-1)**k) * coeffs[k] * y[i - k]
        gl_frac_derv[i] = s / (h**alpha)

    return gl_frac_derv

# Define Riemann-Liouville fractional derivative function for evenly spaced data.
def rl_fractional_derivative(x, y, alpha=0.5): # X is the variable and Y is the discrete function to be analyzed. 
                                               # Alpha is the order of the operation, and, if a value is not given 
                                               # when the function is used, a semi-derivative (alpha = 0.5) will be calculated.

    # Turn the X and Y vectors of data into a Numpy array (if they are not already np arrays)
    x = np.array(x)
    y = np.array(y)
    # Calculate the step in the X array
    h = x[1] - x[0]
    # Check if the values in the X array are evenly spaced
    if not np.allclose(np.diff(x), h):
        raise ValueError("X values must be evenly spaced.")
        
    # Define n as the length of the X array
    n = len(x)
    
    # Create an n-sized array for the calculated integral term. They will be placed in the array later
    I = np.zeros(n)

    # Approximate the integral using rectangular integration
    for j in range(n):
        s = 0
        for k in range(j + 1):
            s += y[k] * (x[j] - x[k])**(alpha-1) if j != k else 0
        I[j] =  s * h / gamma(alpha)

    # Differentiate I using central differences
    rl_frac_deriv = np.gradient(I, h)
    return rl_frac_deriv

# Define equally spaced values for theta function
n = 1
T = 298 # K
F = 96485 # C/mol
R = 8.314 # J/(K*mol)
delta_E = 0.5 # Delta_E = Ei - E° in V
scan_rate = 0.05 # V/s
a = (n*F*scan_rate)/(R*T) # 
b = (n*F*delta_E)/(R*T) # For Dr = Dox
time = np.linspace(0, 20, 20001) 
theta = -a*time+b

# Calculate concentration profile sigmoidal curve
y = 1/(1+np.exp((theta)))

# Obtain semi-derivative by numerical calculations
y_half_deriv = gl_fractional_derivative(a*time, y, alpha=0.5)
y_half_rl = rl_fractional_derivative(a*time, y, alpha=0.5)

# Plot results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(theta, y)
ax1.set_title('Concentration Profile', fontsize=18, fontweight='bold')
ax1.set_xlabel(r'$\theta$',fontsize=16, fontweight='bold')
ax1.set_ylabel(r'$\frac{1}{1 + \exp(\theta)}$',fontsize=16, fontweight='bold')
ax1.set_xlim([-20, 20])
ax1.set_ylim([-0.05, 1.1])
ax1.set_xticks([-20,-15, -10, -5, 0, 5, 10, 15,20])
ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

ax2.plot(theta, -y_half_rl, label='Riemann-Liouville', linestyle='--')
ax2.plot(theta, -y_half_deriv, label='Grünwald–Letnikov', linestyle='--')
ax2.set_title('Semi-Derivatives', fontsize=18, fontweight='bold')
ax2.set_xlabel(r'$\theta$',fontsize=16, fontweight='bold')
ax2.set_ylabel(r'$-\frac{d^{\frac{1}{2}}}{d(at)^{\frac{1}{2}}}\left[\frac{1}{1 + \exp(\theta)}\right]$',fontsize=16, fontweight='bold')
ax2.legend()
ax2.set_xlim([-20, 20])
ax2.set_ylim([-0.5, 0.05])
ax2.set_xticks([-20,-15, -10, -5, 0, 5, 10, 15,20])
ax2.set_yticks([0.0, -0.1, -0.2, -0.3, -0.4, -0.5])

plt.tight_layout()
plt.show()
