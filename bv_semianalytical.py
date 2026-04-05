# Simulation of proton couple electron transfer in electrocatalysis:

# 1 electron and 1 proton transfers in buffered and non-buffered systems according to Costentin, C., 2020. Available at 'doi.org/10.1021/acscatal.0c02532

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

options = ["1. ETPT buffered", "2. ETPT non-buffered", "3. PTET buffered", "4. PTET non-buffered", "5. CPET buffered", "6. CPET non-buffered",
           "7. Quit program"]
plt.figure()
clear_plot = False

# Define general parameters of simulation. Should only be edited by the developer.
h = 0.005
n = 3000
a = (4/3)*np.sqrt(h/np.pi)


# List the possible cases.
print("Simulate voltammograms for electrocatalytic proton-coupled electrons transfer reactions, under heterogeneous catalysis conditions")
print("The possible cases are:")
print()

for i in range(len(options)):
    print(options[i])
    
print()

# Define a function to avoid invalid inputs.
def get_valid_number(prompt):
    while True:
        try:
            return int(input(prompt))
        except ValueError:
            print()
            print("Invalid choice. Choose an integer between 1 and 7")
            print()

choice = get_valid_number("Choose number from 1 to 6 to simulate or 7 to quit:")           
# Loop until quit, generating one voltammogram at each loop
while True:
    g = 25
    g_min = -25
    dimless_current = np.zeros(2*n)
    dimless_pot = np.zeros(2*n)
    conv_prev_list = np.zeros(2*n)
    forward_scan = True

    if clear_plot:
        plt.clf()
        clear_plot = False
    
    if choice == 1:
        print()
        print("You choose ETPT buffered simulation. Enter the parameters:")
        print()
        l = float(input("lambda:")) 
        k = float(input("lambda_a:"))
        m = float(input("lamba_-a:"))

                  
        for i in range(2*n):
            if i == 0:
                conv = 0
                psi = 0
            else:
                if forward_scan:
                    g -= h  # Forward scan (decreasing g)
                    if g <= g_min:
                        forward_scan = False  # Reverse at minimum potential
                else:
                    g += h  # Reverse scan (increasing g)
                conv = np.sqrt(h/np.pi)*sum((2/3)*(dimless_current[j-1]-dimless_current[j])*((i-j+1)**1.5-(i-j)**1.5)+2*(dimless_current[j]*(i-j+1)
                                                                    -dimless_current[j-1]*(i-j))*(np.sqrt(i-j+1)-np.sqrt(i-j)) for j in range(i-1))
                f = (a/2)*psi + conv
                psi = (1/(4 * a * (np.exp(g) * k + k - m)))*(-np.sqrt((2 * a * l * m - 2 * f * np.exp(g) * k - 2 * f * k + 2 * f * m + np.exp(g) * k + 2 * np.exp(g) * l + k + 2 * l + m)**2
                - 4 * (2 * f * l * m + l * m) * (-2 * a * np.exp(g) * k - 2 * a * k + 2 * a * m)) + 2 * a * l * m - 2 * f * np.exp(g) * k - 2 * f * k + 2 * f * m + np.exp(g) * k + 2 * np.exp(g) * l + k + 2 * l + m)

                dimless_current.append(psi)
                dimless_pot.append(g)
                conv_prev_list.append(f)
        # print(dimless_current) # Verify the output dimless_current
                

        plt.plot(dimless_pot,dimless_current)
        #plt.plot(dimless_pot,conv_prev_list)
        plt.show(block=False)
        plt.pause(0.1)
        print()

    elif choice == 2:
        print()
        print("You choose ETPT non-buffered simulation. Enter the parameters:")
        print()
        l = float(input("lamba_ah:"))
        w = float(input("lambda_w:"))
        k = float(input("lamda_-ah:"))

        for i in range(2*n):
            if i == 0:
                conv = 0
                psi = 0
            else:
                if forward_scan:
                    g -= h  # Forward scan (decreasing g)
                    if g <= g_min:
                        forward_scan = False  # Reverse at minimum potential
                else:
                    g += h  # Reverse scan (increasing g)

                conv = np.sqrt(h/np.pi)*sum((2/3)*(dimless_current[j-1]-dimless_current[j])*((i-j+1)**1.5-(i-j)**1.5)+2*(dimless_current[j]*(i-j+1)
                                                                    -dimless_current[j-1]*(i-j))*(np.sqrt(i-j+1)-np.sqrt(i-j)) for j in range(i-1)) 
                f = (a/2)*psi + conv
                psi = (1/(2*a*k))*(+np.sqrt((a*k*w+f*k+np.exp(g)*l+np.exp(g)*w+k+l+w)**2 -4*a*k*(f*k*w+k*w))-a*k*w-f*k-np.exp(g)*l-np.exp(g)*w-k-l-w) 
                dimless_current.append(psi)
                dimless_pot.append(g)
                conv_prev_list.append(f)

        #print(dimless_pot)
        plt.plot(dimless_pot,dimless_current)
        #plt.plot(dimless_pot,conv_prev_list)
        plt.show(block=False)
        plt.pause(0.1)
        print()  

    elif choice == 3:
        print()
        print("You choose PTET buffered simulation. Enter the parameters:")
        print()
        l = float(input("lambda:"))
        k = float(input("lambda_a:"))
        m = float(input("lambda_-a:"))
      
        for i in range(2*n):
            if i == 0:
                conv = 0
                psi = 0
            else:
                if forward_scan:
                    g -= h  # Forward scan (decreasing g)
                    if g <= g_min:
                        forward_scan = False  # Reverse at minimum potential
                else:
                    g += h  # Reverse scan (increasing g)

                conv = np.sqrt(h/np.pi)*sum((2/3)*(dimless_current[j-1]-dimless_current[j])*((i-j+1)**1.5-(i-j)**1.5)+2*(dimless_current[j]*(i-j+1)
                                                                    -dimless_current[j-1]*(i-j))*(np.sqrt(i-j+1)-np.sqrt(i-j)) for j in range(i-1))
                f = (a/2)*psi + conv
                psi = (1/(4 * a * (-np.exp(g) * k + np.exp(g) * m + m)))*(+np.sqrt((2 * a * l * m - 2 * f * np.exp(g) * k + 2 * f * np.exp(g) * m + 2 * f * m + np.exp(g) * k + np.exp(g) * m + 2 * l + m)**2
                - 4 * (2 * f * l * m + l * m) * (-2 * a * np.exp(g) * k + 2 * a * np.exp(g) * m + 2 * a * m)) - 2 * a * l * m + 2 * f * np.exp(g) * k - 
                2 * f * np.exp(g) * m - 2 * f * m - np.exp(g) * k - np.exp(g) * m - 2 * l - m)
                
                dimless_current[i] = psi
                dimless_pot[i] = g
                conv_prev_list[i] = f
                
        plt.plot(dimless_pot,dimless_current)
        #plt.plot(dimless_pot,conv_prev_list)
        plt.show(block=False)
        plt.pause(0.1)
        print()
        
    elif choice == 4:
        print()
        print("You choose PTET non-buffered simulation. Enter the parameters:")
        print()
        L = float(input("Lambda:"))
        alpha = float(input("alpha:"))
        # k = float(input("lamda_-a:"))
        
        for i in range(2*n):
            if i == 0:
                conv = 0
                psi = 0
            else:
                if forward_scan:
                    g -= h  # Forward scan (decreasing g)
                    if g <= g_min:
                        forward_scan = False  # Reverse at minimum potential
                else:
                    g += h  # Reverse scan (increasing g)
                conv = np.sqrt(h/np.pi)*sum((2/3)*(dimless_current[j-1]-dimless_current[j])*((i-j+1)**1.5-(i-j)**1.5)+2*(dimless_current[j]*(i-j+1)
                                                                    -dimless_current[j-1]*(i-j))*(np.sqrt(i-j+1)-np.sqrt(i-j)) for j in range(i-1))
                f = (a/2)*psi + conv
                psi =-((-1 + f + np.exp(g) * f) * L * alpha) / (np.exp(g * alpha) + a * L * alpha + a * np.exp(g) * L * alpha)
                
                dimless_current[i] = psi
                dimless_pot[i] = g
                conv_prev_list[i] = f
                
        #print(dimless_pot)
        plt.plot(dimless_pot,dimless_current)
        #plt.plot(dimless_pot,conv_prev_list)
        plt.show(block=False)
        plt.pause(0.1)
        print()  

    elif choice == 5:
        print()
        print("You choose CPET buffered simulation. Enter the parameters:")
        print()
        l = float(input("lambda:"))
        L = float(input("LAMBDA:"))
        c = float(input("alpha:"))
        
        for i in range(2*n):
            if i == 0:
                conv = 0
                psi = 0
            else:
                if forward_scan:
                    g -= h  # Forward scan (decreasing g)
                    if g <= g_min:
                        forward_scan = False  # Reverse at minimum potential
                else:
                    g += h  # Reverse scan (increasing g)
                conv = np.sqrt(h/np.pi)*sum((2/3)*(dimless_current[j-1]-dimless_current[j])*((i-j+1)**1.5-(i-j)**1.5)+2*(dimless_current[j]*(i-j+1)
                                                                    -dimless_current[j-1]*(i-j))*(np.sqrt(i-j+1)-np.sqrt(i-j)) for j in range(i-1))
                f = (a/2)*psi + conv
                psi = (1/(4 * a * L))*(+np.sqrt((2 * a * l * L + 2 * l * np.exp(c * g) + 2 * f * L + 2 * np.exp(g) * L + L)**2
                - 8 * a * L * (2 * f * l * L + l * L)) - 2 * a * l * L - 2 * l * np.exp(c * g) - 2 * f * L - 2 * np.exp(g) * L - L)
                
                dimless_current.append(psi)
                dimless_pot.append(g)
                conv_prev_list.append(f)
                
        plt.plot(dimless_pot,dimless_current)
        #plt.plot(dimless_pot,conv_prev_list)
        plt.show(block=False)
        plt.pause(0.1)
        print()

    elif choice == 6:
        print()
        print("You choose CPET buffered simulation. Enter the parameters:")
        print()
        w = float(input("lambda_w:"))
        L = float(input("LAMBDA:"))
        c = float(input("alpha:"))
        
        for i in range(2*n):
            if i == 0:
                conv = 0
                psi = 0
            else:
                if forward_scan:
                    g -= h  # Forward scan (decreasing g)
                    if g <= g_min:
                        forward_scan = False  # Reverse at minimum potential
                else:
                    g += h  # Reverse scan (increasing g)
                conv = np.sqrt(h/np.pi)*sum((2/3)*(dimless_current[j-1]-dimless_current[j])*((i-j+1)**1.5-(i-j)**1.5)+2*(dimless_current[j]*(i-j+1)
                                                                    -dimless_current[j-1]*(i-j))*(np.sqrt(i-j+1)-np.sqrt(i-j)) for j in range(i-1))
                f = (a/2)*psi + conv
                psi = (1/(a*L))* (0.5 * (+np.sqrt(a**2 * L**2 * w**2 + 2 * a * L * w**2 * np.exp((c - 1) * g + g) - 2 * a * f * L**2 * w + 2 * a * np.exp(g) * L**2 * w
                 - 2 * a * L**2 * w + 2 * f * L * w * np.exp((c - 1) * g + g) + 2 * L * w * np.exp((c - 1) * g + g) + 2 * L * w * np.exp((c - 1) * g + 2 * g)
                + w**2 * np.exp(2 * (c - 1) * g + 2 * g) + f**2 * L**2 + 2 * f * np.exp(g) * L**2 + 2 * f * L**2 + 2 * np.exp(g) * L**2 + np.exp(2 * g) * L**2
                + L**2) - a * L * w - w * np.exp((c - 1) * g + g) - f * L - np.exp(g) * L - L))

                dimless_current.append(psi)
                dimless_pot.append(g)
                conv_prev_list.append(f)


       
        plt.plot(dimless_pot,dimless_current)
        #plt.plot(dimless_pot,conv_prev_list)
        plt.show(block=False)
        plt.pause(0.1)
        print()

    elif choice == 7:
        
        break # Quit if choice is 7.

        
    else: # Alert input invalidity.
        print()
        print("Invalid choice. Choose an integer between 1 and 7.")
        print()
    
    if input("Do you want to export the last voltammogram as .csv file? (y/n): ").strip().lower() == 'y':
        file_name = str(input("Name the file:"))

        df = pd.DataFrame({"Potential (g)": dimless_pot, "Current (psi)": dimless_current})
        df.to_csv(f"{file_name}.csv", index=False)

    if input("Do you want to clear the plot before the next simulation? (y/n): ").strip().lower() == 'y':
        clear_plot = True