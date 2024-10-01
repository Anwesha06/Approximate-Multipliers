import itertools
from operator import*
import math
import numpy as np

def decimal_to_bin(n, bits=8):
    return list(map(int, bin(n)[2:].zfill(bits)))

def bin_to_decimal(bin_list):
    return int(''.join(map(str, bin_list)), 2)

def approximate_multiplier(A, B, debug=False):
    A_bin = decimal_to_bin(A)
    B_bin = decimal_to_bin(B)
    
    partial_products = []
    
    for i in range(7,-1,-1):
        partial_product = [0] * (7-i+1) + [A_bin[7-i] & B_bin[j] for j in range(8)] + [0] * i
        partial_product = (partial_product + [0] * (16 - len(partial_product)))[:16]
        partial_products.append(partial_product)
        
        if debug:
            print(f"Partial product row {i}: {''.join(map(str, partial_product))} ({bin_to_decimal(partial_product)})")
    

    result = [0] * 16
   
    for i in range(15,7,-1):
        result[i]=partial_products[0][i]
        for j in range(0,8):
            result[i] |= partial_products[j][i]
            
        if debug:
            print(result[i])
    
    carry=0
    for i in range(7, -1,-1):
    
        sum_=carry
        for j in range(0,8):
            
            sum_ +=  partial_products[j][i]
        carry = sum_ //2
        result[i]=sum_ %2
        
        
        if debug:
            print(result[i])
    
    
    return bin_to_decimal(result)


def exact_multiplier(A, B):
    return A * B

# Parameters
in_size = 7
out_size = 15
zero = 0


totalSamples = 65536
emax = 32768  # (-32768 to 32767)


error = 0
error_distance = 0
error_distance_abs = 0
error_distance_square = 0
sum_ED = 0
sum_ED_abs = 0
sum_SED = 0
sum_RE = 0
max_error = 0


results = []


for i in range(256):
    A = i
    for j in range(256):
        B = j


        exact = exact_multiplier(A,B)
        apprx= approximate_multiplier(A,B)
        
        if exact > apprx:
            error_distance_abs = exact - apprx
        else:
            error_distance_abs = apprx - exact
            
        if exact != apprx:
            error += 1
            error_distance = error_distance_abs
            error_distance_square = ((error_distance_abs) ** 2)


        ER = error / totalSamples

        


        if error_distance_abs > max_error:
            max_error = error_distance_abs


        sum_ED += error_distance
        sum_ED_abs += error_distance_abs
        sum_SED += error_distance_square


        MED = sum_ED_abs / totalSamples
        MSED = sum_SED / totalSamples
        NMED = MED / emax
        NMSED = MSED / emax
        WCE = max_error
        NoEB = math.sqrt(MSED)
        NM = MED / emax

        if exact != 0:
            RE = error_distance_abs / exact
            sum_RE += RE
            MRED = sum_RE / totalSamples
        else:
            MRED = 0


        results.append({
            "A": A,
            "B": B,
            "exact": exact,
            "apprx": apprx,
            "error": error,
            "MED": MED,
            "MRED": MRED,
            "MNED": NMED,
            "ER": ER,
            "max": max_error,
            "MSED": MSED,
            "NMED": NMED,
            "NMSED": NMSED,
        })


with open("multiplier_1.txt", "w") as f:
    for r in results:
        f.write(
            f"A={r['A']} B={r['B']} exact={r['exact']} apprx={r['apprx']} #error={r['error']} MED={r['MED']} "
            f"MRED={r['MRED']} MNED={r['MNED']} ER={r['ER']} max={r['max']} MSED={r['MSED']} NMED={r['NMED']} NMSED={r['NMSED']}\n"
        )

print("Simulation complete. Check multiplier_1.txt for results.")
