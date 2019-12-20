import numpy
import scipy

"""
def function1(x,y):
    a,b=scipy.Symbol('x y')
    return y-2*x
"""
#Bounderies for the plot function
x_min=float(-1)
x_max=float(1)
y_min=float(-1)
y_max=float(1)
class fractals:
    epsilon=float(1e-7)
    delta=0.001
    iterations=1000
    #Takeas a startiing argument, being a function and/or it's derivitive
    def __init__(self,function,prime=None):
        if prime==None:
            self.prime=scipy.diff(function)
            self.function=function
        else:
            self.prime=prime
            self.function=function
        
    
    #Might not need repr.
    def __repr__(self):
        pass

    solutions=[]
    colours=[]
    
    def newton_method(self,guess):
        guess=complex(guess)
        x=guess
        for i in range(1000):#Should be iterations
            x=x-self.function(x)/self.prime(x)
            if abs(x-guess)<epsilon:
                print("We have convergeance")
                solutions.append(x)
                break
            return x
    

    def draw(self):
        pass


def p(z):
    x=z[0]
    y=[1]
