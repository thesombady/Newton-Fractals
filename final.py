import scipy
import sympy
import math
import matplotlib.pyplot as plt
import numpy as np


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
                self.zeros=np.array([[ ]])
        if div==None:
            self.function=function
            self.div=self.numerical_div
            self.epsilon=epsilon
            self.delta=delta
            self.iterations=iterations
            self.zeros=np.array([[ ]])
    
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
            raise TypeError("This is not of correct type")
        if np.linalg.det(self.div(guess[0],guess[1]))==0:
            return None, self.iterations
        try:
            jacobian_inverse=np.linalg.inv(self.div(guess[0],guess[1]))
            for i in range(self.iterations):
                guess=guess+np.matmul(-jacobian_inverse,self.function(guess[0],guess[1]))
                if np.linalg.norm(self.function(guess[0],guess[1]))<self.delta:
                    return guess, i
        except:
            return None
    
    def find_zeros(self,guess,which_newton=0):
        value,i=self.newtonmethod(guess)
        tolerance=1e-5
        if value==None:
            return -1,i#This is because the value didn't converge
        if self.zeros.size==0:
            self.zeros=np.array(value)
        t=(self.zeros-value)**2
        distance=t[:,0] + t[:,1]
        if (distance-tolerance).any():
            return np.where(distance<tolerance)[0][0] , i
        else:
            return np.reshape(np.append(self.zeros,value,value),(-1,2))#This line adds the zeros
    
    def call_find_zeros(self,x,y,simple=False):
        """This calls for the function 'find_zeros'"""
        return (self.find_zeros((x,y),simple))

    def plot(self,N,a,b,c,d,which_newton=None):
        if which_newton==None:
            self.newton=self.newtonmethod
        elif which_newton==1:
            self.newton=self.simple_newtonmethod
        else:
            self.newton=self.newton
        """
        Now just make a plot function with newton() instead of any function, this is just so that the thing cheacks out.
        """
        X,Y=np.meshgrid(np.linspace(a,b,N),np.linspace(c,d,N))
        Z=self.newton(X,Y)
        cs=plt.contour(X,Y,Z,np.logspace(0,3.5,7,base=10),cmap='gray')
        rosen=lambda x:Z(x[0],x[1])
        solution,iterates=scipy.optimize.fmin_powell(rosen,x0=np.array([0,-0.7]),retall=True)
        x,y=zip(*iterates)
        plt.plot(x,y,'ko') # plot black bullets
        plt.plot(x,y,'k:',linewidth=1) # plot black dotted lines
        plt.title("Newton Fractals")
        plt.clabel(cs)

    def number(self,guess):
        if not isinstance(guess,(tuple,list,np.array)):
            raise TypeError("Iterations method cannot be made!")
        value, i=self.simple_newtonmethod(guess)
        if isinstance(value,np.array):
            return value,i
        else:
            return self.newtonmethod(guess)

def f(x,y):
    return np.array([[x**3-3*x*y**2-1],[3*x**2*y-y**3]])

print(fractal2d(f).plot)