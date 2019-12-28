import scipy
import sympy
import math
import cmath
import matplotlib.pyplot as plt
import numpy as np
import numdifftools


class fractal2d:
    def __init__(self,function,div=None,epsilon=1e-8,delta=0.001,iterations=1000):
        """Intalizes a script with a function and possibly it's derivate"""
        if not callable(function):
            raise TypeError("The function isn't callable")
        if div!=None:
            if not callable(div):
                raise ValueError("The derivaite isn't callable.")
            else:
                self.function=function
                self.div=div
                self.epsilon=epsilon
                self.delta=delta
                self.iterations=iterations
        if div==None:
            self.function=function
            self.div=self.numerical_div
            self.epsilon=epsilon
            self.delta=delta
            self.iterations=iterations
    """
    def partialdiv(self,x,y):
        #Returns a jacobian matrix
        J = numdifftools.Jacobian(lambda z: self.function(z.reshape(x.shape), y).ravel())
        return J(x.ravel()).reshape(x.shape)
    """
    def divf(self,x,y):
        """This returns the jacobian matrix """
        x,y=sympy.symbols('x y')
        a=sympy.diff(self.function(x,y)[0],x)
        b=sympy.diff(self.function(x,y)[0],y)
        c=sympy.diff(self.function(x,y)[1],x)
        d=sympy.diff(self.function(x,y)[1],y)
        return np.array([[a,b],[c,d]]) 
    
    def numerical_div(self,guess):
        """This will return a numerical derivision of the self.function at a point 'guess'.
        -----------------------------------------------------------------------------------
        .-This is done by a numerical esimation of guess + some delta """
        if not isinstance(guess,(tuple,list,np.array)):
            raise TypeError("The derivision cannot be made due to wrong type of format")
        f1_x=(self.function(guess +np.array([self.delta,0]))[0])/(self.delta)
        f1_y=(self.function(guess +np.array([0,self.delta]))[0])/(self.delta)
        f2_x=(self.function(guess +np.array([self.delta,0]))[1])/(self.delta)
        f2_y=(self.function(guess +np.array([0,self.delta]))[1])/(self.delta)
        return np.array([[f1_x,f1_y],[f2_x,f2_y]]).T


    def newtonmethod(self,guess):
        if not isinstance(guess,(tuple,list,np.array)):
            raise TypeError("The input is of the wrong format.")
        else:
            guess=np.array(guess)
            if np.linalg.det(self.div(guess[0],guess[1]))==0:
                return None, self.iterations
            for i in range(self.iterations):
                xy=guess
                guess=guess-np.linalg.solve(self.div(guess),-self.function(guess))
                x=guess[0]
                y=guess[1]
                if abs(x-xy[0])<self.epsilon and abs(y-xy[1])<self.epsilon:#Yay, it converged
                    return guess,i
            else:
                return None,self.iterations#It didn't converge
    
    def simple_newtonmethod(self,guess):
        if not isinstance(guess,(tuple,list,np.array)):
            raise TypeError("The input is not of correct type")
        if np.linalg.det(self.div(guess[0],guess[1]))==0:
            return None, self.iterations
        Jacobian_inverse = np.linalg.inv(self.div(guess[0],guess[1]))
	    for i in range(self.iterations):
            try:
				guess = np.matmul(-Jacobian_inverse,self.function(guess[0],guess[1]) + guess
				if np.linalg.norm(self.function(guess)) < self.epsilon:
					return guess, i
            except:
                return None, self.iterations
        return None, self.iterations








def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])


"""
Partial derivitive might have to be worked on again
"""

