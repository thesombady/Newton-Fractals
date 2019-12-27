import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
            #Make a method that does self.div=div
            self.iterations=iterations
            self.delta=delta
            self.epsilon=epsilon
            self.tol=tol
            self.zeros=np.array([[ ]])
        else:
            self.function=function
            self.div=div
            self.iterations=iterations
            self.delta=delta
            self.epsilon=epsilon
            self.tol=tol
            self.zeros = np.array([[ ]])
       
    def partialdiv(self,x,y):
        if isinstance(x,(float,int)) and isinstance(y,(float,int)):
            f=self.function
            jacobian=scipy.zeros(2)
            jacobian=1/self.epsilon*np.array([f(x+self.epsilon,y)-f(x,y),f(x,y+self.epsilon)-f(x,y)])
            return jacobian
        else:
            raise ValueError("The partial derivitive cannot be calculated")
    
    def __repr__(self):
        pass

    def partial_div(self,x,y):
        if isinstance(x,(float,int)) and isinstance(y,(float,int)):
            pass
        else:
            raise ValueError("The value is not of right type, only int and floats are allowed.")
    def newtonmethod(self,guess):
        if isinstance(guess,(tuple,list,np.array)):
            guess=np.array(guess)
        else:
            raise KeyError("Wrong type of format")
        
        for i in range(self.iterations):
            
            jacobian_matrix=self.div(guess[0],guess[1])
            if np.linalg.det(jacobian_matrix)==0:
                return None, self.iterations
            else:
                jacobian_matrix_inv=np.linalg.inv(jacobian_matrix)
            guess_new=guess-np.matmul(jacobian_matrix_inv,self.function(guess[0],guess[1]))
            x=guess_new[0]
            y=guess_new[1]
            if abs(x-guess[0])<self.epsilon and abs(y-guess[1])<self.epsilon:
                return guess_new,i
            #It converged
            guess=guess_new
        else:
            return None, self.iterations
    
    def locate_zeros(self,guess,which_newton=0):
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

    def simple_newtonmethod(self,guess):
        if isinstance(guess,(list,np.array,tuple)):
            guess=np.array(guess)
        else:
            raise KeyError("Wrong type of format")
        guess_new=np.array(0,0)
        jacobian_matirx=self.div(guess[0],guess[1])
        if np.linalg.det(jacobian_matirx)==0:
            return None
        else:
            jacobian_matirx_inv=np.linalg.inv(jacobian_matirx)
        for i in range(self.iterations):
            guess_new=guess-np.matmul(jacobian_matirx_inv,self.function(guess[0],guess[1]))
            x=guess_new[0]
            y=guess_new[1]
            if abs(guess_new)>1/self.epsilon:
                return None, self.iterations
            if abs(x-guess[0])<self.epsilon and abs(y-guess[1])<self.epsilon:
                return guess_new,i
            guess=guess_new
        else:
            return None, self.iterations

    def optimized(self,guess):
        value,i=self.simple_newtonmethod(guess)
        if not isinstance(value,(np.array)):
            return self.newtonmethod(guess)
        else:
            return value,i
    
    def call_locace_zeros(self,a,b,simple=False):
        return self.locate_zeros((a,b),simple)
    
    def plot(self,N,a,b,c,d,which_newton=0,showplot=True):
        """
        if which_newton==0:
            self.newton=self.simple_newtonmethod
        elif which_newton==1:
            self.newton=self.optimized
        else:
            self.newton=self.newtonmethod
        """
        Y,X=np.meshgrid(np.linspace(a,b,N),np.linspace(c,d,N),indexing='ij')
        vector_zeros=scipy.vectorize(self.call_locace_zeros)

        A,B=vector_zeros(X,Y,self.simple_newtonmethod)
        while showplot==True:
            colour_map=plt.cm.get_cmap('plasma')
            colour_map.set_under('0')
            pcol=plt.pcolormesh(X,Y,A,colour_map,edgecolor='face',linewidth=0,rasterized=True,vmin=0,vmax=A.max())
            plt.colorbar()
            custom = [(0,(0, 0, 0, 0)), (1,(0, 0, 0, 0.5))]
            new_cmap=colors.LinearSegmentedColormap.from_list('something',custom)
            pcol=plt.pcolormesh(X,Y,B,edgecolor='face',cmap=new_cmap,linewidth=0,rasterized=True)
            plt.show()