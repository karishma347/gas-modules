
def bisection_NSW(a,b,x,g,func):   #a and b are guessing values for bisection method
    if (func(a,x,g) * func(b,x,g)>=0):
        print("wrong assumption")
        return
    c=a
    while((b-a)>=0.01):
        c= (a+b)/2
        if(func(c,x,g) == 0.0):
            break
        if(func(c,x,g)*func(a,x,g)<0):
            b=c
        else:
            a=c
    return (a+b)/2


def bisection_OSW(a,b,x,g,func,beta):   #a and b are guessing values for bisection method for OBW relations
    if (func(a,x,g,beta) * func(b,x,g,beta)>=0):
        print("wrong assumption")
        return
    c=a
    while((b-a)>=0.01):
        c= (a+b)/2
        if(func(c,x,g,beta) == 0.0):
            break
        if(func(c,x,g,beta)*func(a,x,g,beta)<0):
            b=c
        else:
            a=c
    return (a+b)/2