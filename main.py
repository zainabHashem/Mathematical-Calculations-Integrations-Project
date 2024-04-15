#libraries
import math
import cmath
import numpy as np
from scipy.integrate import quad
import sympy as sp
i = 1j

#choose function
#f(x) = cos(x)+6*i*sin(2*x)
#f1(x) = sin(x)*e^(x)
#f2(x) = e^(x)*(x^2)
def f(x):
  return math.cos(x)+6*i*math.sin(2*x)
def f1(x):
  return math.sin(x)*cmath.exp(x)
def f2(x):
  return cmath.exp(x)*x**2

#complete Integration equation
#[f(x) or f1(x) or f2(x)] * e^(iwx)
def integ_func(x,chosen_function,w):
    if chosen_function == 'f':
        return f(x) * cmath.exp(i * w * x)
    elif chosen_function == 'f1':
        return f1(x) * cmath.exp(i * w * x)
    elif chosen_function == 'f2':
        return f2(x) * cmath.exp(i * w * x)
    else:
        raise ValueError("invalid function choice")

#Calculating the integral
def integral(a,b,n,chosen_function,w):
    delta_x = (b-a)/n
    integral_sum = 0
    for i in range(n+1):
        x_i=a+(i)*delta_x
        f_x = integ_func(x_i,chosen_function,w)
        if i == 0 or i == n:  # For the first and last terms
            integral_sum += 0.5*f_x*delta_x
        else: 
            integral_sum += f_x*delta_x
    return integral_sum

#enter choose function[f,f1,f2] & a , b , w , n
chosen_function=input("choose a function (f, f1, f2): ")
if chosen_function not in ['f','f1','f2']:
    raise ValueError("invalid function choice")
a = float(input("enter the lower limit (a): "))
b = float(input("enter the upper limit (b): "))
w = float(input("enter the the value of (w): "))
n = int(input("enter the Number of intervals (n): "))


#without use built-in function
#print result in Cartesian form
result = integral(a,b,n,chosen_function,w)
real = result.real
imagine = result.imag
print(f"(Real part): {real:.6f}")
print(f"(Imaginary part): {imagine:.6f}")
print(f"integral in Cartesian form: {result:.6f}")

#use built-in function quad
result_real, _ = quad(lambda x: integ_func(x,chosen_function,w).real,a,b)
result_imag, _ = quad(lambda x: integ_func(x,chosen_function,w).imag,a,b)
result = result_real+result_imag*i
print(f"integral using built in function: {result:.6f}")

#use analytical form with f(x) = cos(x)+6*i*sin(2*x)
x = sp.symbols('x')
f = (sp.cos(x) + 6* i * sp.sin(2*x)) * sp.exp(1j * 4 * x)
analytical_integral = sp.integrate(f,x)
print("analytical form: ",sp.N(analytical_integral))
analytical_result = analytical_integral.subs(x,b) - analytical_integral.subs(x,a)
print("result = ",sp.N(analytical_result))

#print result in exponential & trigonometric polar form
magnitude = abs(result)
phase = math.atan2(result.imag,result.real)
print(f"phase: {phase:.6f}")
angle_degrees = phase*180/np.pi
#angle_degrees = math.degrees(phase)
print(f"exponential polar form: {magnitude:.6f} * e^(i * {angle_degrees:.6f}°)")
print(f"trigonometric Polar form: {magnitude:.6f}(cos({angle_degrees:.6f})+i*sin({angle_degrees:.6f}))")
polar = cmath.polar(result)
print(f"polar using built in function: {polar[0]:.6f}*e^(i*{math.degrees(polar[1]):.6f}°)")

#Trigonometric Fourier series
#function calculate a0,ak,bk
#f(x) = cos(x)+6*i*sin(2*x)
#f1(x) = x^2
def f(x):
    return math.cos(x)+6*i*math.sin(2*x)
def f1(x):
    return x**2

def a_0(x,w):
    return f(x) * cmath.exp(i * w * x)
def a_k(x,w,w_0,k):
    return f(x) * cmath.exp(i * w * x)*cmath.cos(k*w_0*x)
def b_k(x,w,w_0,k):
    return f(x) * cmath.exp(i * w * x)*cmath.sin(k*w_0*x)


#Calculating the integral
def integral(a,b,n,w,w_0,k,flag):
    delta_x = (b - a) / n
    integral_sum = 0
    for i in range(n + 1):
        x = a + i * delta_x
        if flag ==0 :
            f = a_0(x,w)
        if flag ==1 :
            f = a_k(x,w,w_0,k)
        if flag ==2 :
            f = b_k(x,w,w_0,k)

        if i == 0 or i == n:  # For the first and last terms
            integral_sum += 0.5 * f * delta_x
        else:
            integral_sum += f* delta_x
    return integral_sum


# function calculate real & imagine parts in a0,ak,bk
def a_0_real(x, w):
    return np.real(f(x) * np.exp(i * w * x))

def a_0_imag(x, w):
    return np.imag(f(x) * np.exp(i * w * x))

def a_k_real(x, w, w_0, k):
    return np.real(f(x) * np.exp(i * w * x)*np.cos(k*w_0*x))

def a_k_imag(x, w, w_0, k):
    return np.imag(f(x) * np.exp(i * w * x)*np.cos(k*w_0*x))

def b_k_real(x, w, w_0, k):
    return np.real(f(x) * np.exp(i * w * x)*np.sin(k*w_0*x))

def b_k_imag(x, w, w_0, k):
    return np.imag(f(x) * np.exp(i * w * x)*np.sin(k*w_0*x))
  

#enter T,w,n 
T = float(input("enter the period (T): "))
w = float(input("enter the the value of (w): "))
n = int(input("enter the Number of intervals (n): "))

# without use built-in function
T *= math.pi
w_0 = (2*math.pi)/T
k=3
a0 = (1/T) *integral(0,T,n,w,w_0,k,0)
ak =  (2/T) *integral(0,T,n,w,w_0,k,1)
bk =  (2/T) *integral(0,T,n,w,w_0,k,2)
print(f"coefficient a0: {a0}")
print(f"coefficient ak: {ak}")
print(f"coefficient bk: {bk}")

#use built-in function quad
real_a0, _ = quad(a_0_real, 0, T,args=(w))
imag_a0, _ = quad(a_0_imag, 0, T,args=(w))
real_ak, _ = quad(a_k_real, 0, T,args=(w,k,w_0))
imag_ak, _ = quad(a_k_imag, 0, T,args=(w,k,w_0))
real_bk, _ = quad(b_k_real, 0, T,args=(w,k,w_0))
imag_bk, _ = quad(b_k_imag, 0, T,args=(w,k,w_0))

a0 = (1/T) *(real_a0 + i * imag_a0)
ak =(2/T) * (real_ak + i * imag_ak)
bk =(2/T) * (real_bk + i * imag_bk)

print(f"results using built in function: ")
print(f"Coefficient a0: {a0}")
print(F"Coefficient ak: {ak}")
print(F"Coefficient bk: {bk}")