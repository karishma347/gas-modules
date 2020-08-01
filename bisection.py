# a python module to calculate bisection method
# Karishma Bharti, Srisha Rao M V, IISc, July 2020
# Global conventions and variables
# a,b are guessing values
# func is the function

def bis(a,b,func):
    if (func(a) * func(b)>=0):
        print("wrong assumption")
        return
    c=a
    while((b-a)>=0.01):
        c= (a+b)/2
        if(func(c) == 0.0):
            break
        if(func(c)*func(a)<0):
            b=c
        else:
            a=c
    return (a+b)/2