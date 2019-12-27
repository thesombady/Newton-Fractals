import scipy
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

class fractal2d:
    """Intalizes a 'script' that plots a fractal for the given function"""
    def __init__(self,function,div=None,epsilon=1e-8,delta=0.001,iterations=100,tol=1e-4):
        """Just the intalization function, there are more variables due to them being inside the
        class, this is to avoid having no constant inside the class named epsilon and so on."""
        if not callable(function) or callable(div):
            raise KeyError("The function is not callable")
        if div==None:
            self.function=function
            #Make a simple derivation for the f(x,y)
            self.iterations=iterations
            self.delta=delta
            self.epsilon=epsilon
            self.tol=tol
        else:
            self.function=function
            self.div=div
            self.iterations=iterations
            self.delta=delta
            self.epsilon=epsilon
            self.tol=tol
            self.zeros = np.array([[ ]])
    
    def __repr__(self):
        pass

    def partial_div(self,x,y):
        if isinstance(x,(float,int)) and isinstance(y,(float,int)):
            pass
        else:
            raise ValueError("The value is not of right type, only int and floats are allowed.")
    def newtonmethod(self,guess):
        if isinstance(guess,(tuple,list,np.array)):
            xy=np.array(0,0)
        else:
            raise KeyError("Wrong type of format")
        
        for i in range(self.iterations):
            
            jacobian_matrix=self.div(guess[0],guess[1])
            if np.linalg.det(jacobian_matrix)==0:
                return None
            else:
                jacobian_matrix_inv=np.linalg.inv(jacobian_matrix)
            guess_new=xy-np.matmul(jacobian_matrix_inv,self.function(xy[0],xy[1]))
            x=guess_new[0]
            y=guess_new[1]
            if abs(x-xy[0])<self.epsilon and abs(y-xy[1])<self.epsilon:
                return guess_new,i
            #It converged
            xy=guess_new
        else:
            return None, self.iterations
    
    def find_zeros(self,guess):
        value,i=self.newtonmethod(guess)
        if value ==None:
            raise ValueError("The solution didn't converge")
        if self.zeros.size==0:
            self.zeros=np.array([value])
        t = (self.zeros-value)**2
        dist = t[:,0] + t[:,1]
        if (dist<self.tol).any() : #zero already exists
            return np.where(dist<self.tol)[0][0] , i
        else:  # value dose not exist yet
            self.zeros=np.reshape(np.append(self.zeros,value),(-1,2)) #add value to zeros
            return self.zeros.size/2-1 , i 

        pass
