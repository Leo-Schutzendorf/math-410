import numpy as np

def Adaptive(a, b, TOL, f, N):
    APP = 0
    
    # Initialize data for the entire interval [a, b]
    a_i = a
    h_i = (b - a) / 2
    FA_i = f(a)
    FC_i = f(a + h_i)
    FB_i = f(b)
    S_i = h_i * (FA_i + (4 * FC_i) + FB_i) / 3
    TOL_i = 10 * TOL  # Initial tolerance scaled up to account for the 1/15 factor
    L_i = 1

    # We use a standard Python list as our LIFO stack for clean append/pop operations
    stack = []
    
    # Pack the initial interval state into a dictionary
    initial_interval = {
        'a': a_i, 'FA': FA_i, 'FC': FC_i, 'FB': FB_i,
        'h': h_i, 'TOL': TOL_i, 'S': S_i, 'L': L_i
    }
    stack.append(initial_interval)

    # Process intervals while there are items remaining on the stack
    while len(stack) > 0:
        # Retrieve and unpack the next subinterval from the stack
        current = stack.pop()
        a_i = current['a']
        FA_i = current['FA']
        FC_i = current['FC']
        FB_i = current['FB']
        h_i = current['h']
        TOL_i = current['TOL']
        S_i = current['S']
        L_i = current['L']

        # Midpoints of the left and right halves
        FD = f(a_i + h_i / 2)
        FE = f(a_i + 3 * h_i / 2)

        # Simpson's approximations for both halves
        S1 = h_i * (FA_i + 4 * FD + FC_i) / 6
        S2 = h_i * (FC_i + 4 * FE + FB_i) / 6

        # Map current variables to your v_ placeholders
        v_1 = a_i
        v_2 = FA_i
        v_3 = FC_i
        v_4 = FB_i
        v_5 = h_i
        v_6 = TOL_i
        v_7 = S_i
        v_8 = L_i

        # Check accuracy constraint (multiplied by 15 or compared to v_6 directly)
        if abs(S1 + S2 - v_7) < v_6:
            APP += S1 + S2
        else:
            if v_8 >= N:
                raise ValueError("Level exceeded")
            else:
                # Prepare data for the right-half subinterval
                right_half = {
                    'a': v_1 + v_5,
                    'FA': v_3,
                    'FC': FE,
                    'FB': v_4,
                    'h': v_5 / 2,
                    'TOL': v_6 / 2,
                    'S': S2,
                    'L': v_8 + 1
                }

                # Prepare data for the left-half subinterval
                left_half = {
                    'a': v_1,
                    'FA': v_2,
                    'FC': FD,
                    'FB': v_3,
                    'h': right_half['h'],        
                    'TOL': right_half['TOL'],    
                    'S': S1,
                    'L': right_half['L']         
                }

                # Push both onto the stack (Right first, then Left so Left is processed first)
                stack.append(right_half)
                stack.append(left_half)
                
    return APP

# Define limits and parameters
a = 0
b = np.pi*0.5

f = lambda x: (6*np.cos(4*x)+4*np.sin(6*x))*np.exp(x)

print("Approximated Integral:", Adaptive(a, b, 0.00001, f, 15))
