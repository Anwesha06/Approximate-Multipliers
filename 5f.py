import itertools
from operator import *
import math
import numpy as np


def dec_to_bin_list(n, bits=8):
    return list(map(int, bin(n)[2:].zfill(bits)))


def bin_list_to_decimal(bin_list):
    return int(''.join(map(str, bin_list)), 2)


def approximate_multiply_L(A, B, debug=False):
    A_bin = dec_to_bin_list(A)
    B_bin = dec_to_bin_list(B)

    partial_products = []

    for i in range(7, 3, -1):
        partial_product = [0] * (i - 3) + [B_bin[i] & A_bin[j] for j in range(8)] + [0] * (7 - i)
        partial_product = (partial_product + [0] * (12 - len(partial_product)))[:12]
        partial_products.append(partial_product)

        if debug:
            print(f"Partial product row {i}: {''.join(map(str, partial_product))} ({bin_list_to_decimal(partial_product)})")

    result = [0] * 12
    result[11] = partial_products[0][11]
    for i in range(10, 4, -1):
        result[i] = partial_products[0][i]
        for j in range(0, 4):
            result[i] |= partial_products[j][i]

        if debug:
            print(result[i])

    sum_L = 0
    carry_L = 0
    for i in range(4,3,-1):
        a = partial_products[0][i]
        b = partial_products[1][i]
        c = partial_products[2][i]
        d = partial_products[3][i]
        sum_L = a | b | c | d | carry_L
        carry_L = (a & b & c) | (b & c & d) | (a & b & d)
        result[i]=sum_L


    carry= carry_L
    for i in range(3, -1, -1):
        sum_ = carry
        for j in range(4):
            sum_ += partial_products[j][i]
        carry = sum_ // 2
        result[i] = sum_ % 2

        if debug:
            print(result[i])

    return bin_list_to_decimal(result)



def approximate_multiply_M(A, B, debug=False):
    A_bin = dec_to_bin_list(A)
    B_bin = dec_to_bin_list(B)

    partial_products = []

    for i in range(3, -1, -1):
        partial_product = [0] * (i + 1) + [B_bin[i] & A_bin[j] for j in range(8)] + [0] * (3 - i)
        partial_product = (partial_product + [0] * (12 - len(partial_product)))[:12]
        partial_products.append(partial_product)

        if debug:
            print(
                f"Partial product row {i}: {''.join(map(str, partial_product))} ({bin_list_to_decimal(partial_product)})")

    result = [0] * 12
    result[11] = partial_products[0][11]
    for i in range(10, 8, -1):
        result[i] = partial_products[0][i]
        for j in range(0, 4):
            result[i] |= partial_products[j][i]

        if debug:
            print(result[i])

    sum_M=0
    carry_M = 0
    for i in range (8,3,-1):
        a=partial_products[0][i]
        b=partial_products[1][i]
        c=partial_products[2][i]
        d=partial_products[3][i]
        sum_M=a^b^c^d^carry_M
        carry_M=(a&b)|(b&c)| (c&d)| (d&a)| (a&c)| (b&d)
        result[i]=sum_M


    carry=carry_M
    for i in range(3, -1, -1):
        sum_ = carry
        for j in range(4):
            sum_ += partial_products[j][i]
        carry = sum_ // 2
        result[i] = sum_ % 2

        if debug:
            print(result[i])

    return bin_list_to_decimal(result)



def approximate_multiply(A, B, debug=False):
    L_ = approximate_multiply_L(A, B, debug) #decimal output from lsb
    M_ = approximate_multiply_M(A, B, debug) #decimal output from msb

    #converting into binary(16-bit)
    L_1= dec_to_bin_list(L_,bits=16)
    M_1 =dec_to_bin_list(M_ <<4, bits=16)


    result = [0]* 16
    for i in range(15,11,-1):
        result[i]=L_1[i]

    for i in range (11,7,-1):
            result[i] = L_1[i] | M_1[i]

    sum_2=0
    carry_2=0
    for i in range(7,4,-1):
        L_2 =L_1[i]
        M_2 =M_1[i]

        sum_2= L_2 ^ M_2 ^ carry_2
        carry_2=(L_2 & M_2)|(M_2 & carry_2)|(carry_2 & L_2)
        result[i]=sum_2

    carry_f = carry_2

    for i in range(4,-1,-1):
        sum_f = carry_f
        sum_f += L_1[i] + M_1[i]
        carry_f = sum_f // 2
        result[i] = sum_f % 2

        if debug:
            print(result[i])



    if debug:
        print(f"Multiplied {A} and {B} -> L: {L}, M: {M}, result: {result}")

    return bin_list_to_decimal(result)


def exact_multiply(A, B):
    return A * B



totalSamples = 0
for i in range(256):
    for j in range(256):
        totalSamples += 1

emax = 32768

# Initializing registers and variables
error = 0
error_distance = 0
error_distance_abs = 0
error_distance_square = 0
sum_ED = 0
sum_ED_abs = 0
sum_SED = 0
sum_RE = 0
max_error = 0
total = 0


results = []


for i in range(256):
    A = i
    for j in range(256):
        B = j
        total += 1


        exact = exact_multiply(A, B)
        apprx = approximate_multiply(A, B)

        if exact > apprx:
            error_distance_abs = exact - apprx
        else:
            error_distance_abs = apprx - exact

        if exact != apprx:
            error += 1
            error_distance = error_distance_abs
            error_distance_square = (error_distance_abs ** 2)


        ER = error / totalSamples


        if error_distance_abs > max_error:
            max_error = error_distance_abs


        sum_ED += error_distance
        sum_ED_abs += error_distance_abs
        if exact != apprx:
            sum_SED += error_distance_square


        MED = sum_ED_abs / totalSamples
        MSED = sum_SED / total
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


with open("5f.txt", "w") as f:
    for r in results:
        f.write(
            f"A={r['A']} B={r['B']} exact={r['exact']} apprx={r['apprx']} #error={r['error']} MED={r['MED']} "
            f"MRED={r['MRED']} MNED={r['MNED']} ER={r['ER']} max={r['max']} MSED={r['MSED']} NMED={r['NMED']} NMSED={r['NMSED']}\n"
        )

print("Simulation complete. Check 5f.txt for results.")