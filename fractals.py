import sympy#Didn't use this because not in the standard liberary.
import numpy as np
import matplotlib.pyplot as plt
import scipy

class fractal2d:
    """A class"""
    def __init__(self,function,derivitive=None,epsilon=1e-8,delta=0.001,iterations=100):#Author: Andreas Evensen
        """The intalization. Using input a two dimensional function if wanted it's derivate. Could also change the iterations, the tolerance and increment. """
        if not callable(function):#If not the function you inputed is callable, do this
            raise KeyError("The function isn't Callable")
        if derivitive!=None:#If we have an input in derivitive, we first test it, if not true, we assign the method partialdiv as the derivate
            if not callable(derivitive):
                self.derivitive=self.partialdiv
                print("The derivitive isn't Callable, using the partialdiv function instead.")
        if derivitive==None:#If no input was detected we assign the method partialdiv as self.derivitive
            self.function=function
            self.derivitive=self.partialdiv
        else:#If we have both input we just assign nem
            self.function=function
            self.derivitive=derivitive
        #This will run anyway. Here we fix some values
        self.epsilon=epsilon
        self.delta=delta
        self.iterations=iterations
        self.zeros=np.array([[]])
    
    def partialdiv(self,x,y):#Author: Kasper Lindqvist
        """Using two input found in guess, to compute the partial derivitive:The jacobian matrix """
        if not isinstance(x,(float,int)) and isinstance(y,(int,float)):
            raise TypeError("the input is not of correct type in partialdiv")
        h=self.epsilon
        Jacobian_matrix=(1/h)*np.array([self.function(x+h,y)-self.function(x,y),self.function(x,y+h)-self.function(x,y)])
        return np.transpose(Jacobian_matrix)#The transpose is because of definition

    def newtonmethod(self,guess):#Author: Andreas Evensen
        """Using a 2-d array, list or tuple as input, as a startpoint for the newton-ralphson Method """
        if not isinstance(guess,(list,tuple,np.ndarray,np.generic)):#Cheacking the inputs-type
            raise TypeError("The input in newtonmethod isn't of correct type!")
        guess=np.array(guess)#Making the input as an array, being necessary as it can be a tuple or list.
        for i in range(self.iterations):
            if np.linalg.det(self.derivitive(guess[0],guess[1]))==0:#Cheacking if the determenant is zero, if so the inverse don't exist, which we use; Morveover if the determenant is zero it has no change to converge
                return None,self.iterations#Diverge
            xy=guess
            guess=guess-np.matmul(np.linalg.inv(self.derivitive(guess[0],guess[1])),self.function(guess[0],guess[1]))#could implement a sigma(self.delta) to improve solution,but would slow calculations down.
            if abs(xy[0]-guess[0])<self.epsilon and abs(xy[1]-guess[1])<self.epsilon:#Testing if it has converged underneath our epsilon
                return guess, i#Converge
        return None,self.iterations #Didn't converge, but not necissary diverged. Might have needed more repititions/iterations

    def simple_newtonmethod(self,guess):#Author: Zebastian
        """Using a 2-darray, list or tuple as input for the simple newton-ralphson Method.
        This method has jacobian matrix calculated only once and using that value throughout. """
        if not isinstance(guess,(tuple,list,np.ndarray,np.generic)):    #Checking if the input is of correct type
            raise TypeError("The input is of the wrong type in simple_newtonmethod!")
        guess=np.array(guess)   #Making the input to a 2-d array
        Jacobian=self.derivitive(guess[0],guess[1]) #Just the jacobian_matrix at point (guess)
        if np.linalg.det(Jacobian)==0:  #This is outside the while loop because it's only has to be computed once
            return None, self.iterations    #Diverges
        Jacobian_inv=np.linalg.inv(Jacobian)
        i=0
        while i<self.iterations:    #Using a while loop as a for loop to indicate that both works just fine
            xy=guess
            guess=guess-np.matmul(Jacobian_inv,self.function(guess[0],guess[1]))
            if abs(guess[0]).any()>(1/(self.epsilon)):
                return None,self.iterations
            if abs(xy[0]-guess[0])<self.epsilon and abs(xy[1]-guess[1])<self.epsilon:
                return guess, i     #Yay converged.
            i+=1
        return None,self.iterations #Didn't converge but not neccissary diverged.
        
    def find_zeros(self,X,Y,which_newton=None):#Author: Colab-Andreas Evensen & Kasper Lindqvist
        """Method that takes two inputs, guess (being a list, tuple or 2-d array) and which_newton which can be either False or something. """
        if not isinstance(X,(list,tuple,np.ndarray,np.generic)) and isinstance(Y,(list,tuple,np.ndarray,np.generic)):#Cheacking the input
            raise TypeError("The input is of wrong type!")
        if which_newton==False:     #Determining which newton-ralphson method to use
            self.newton=self.simple_newtonmethod
        else:
            self.newton=self.newtonmethod
        guess=(X,Y)
        guess=np.array(guess)
        value,i=self.newton(guess)
        if value is None:           #From either newton or simple newton, means no convergence
            return -1,i             #i being the number of iterations
        
        if self.zeros.size == 0:    # If not root existed before this, add this value, "value" as a root
            self.zeros = scipy.array([value])

        q = (self.zeros-value)**2
        distance = q[:,0] + q[:,1]  #Columns and arrays, The indicies represent different columns and arrays
        if (distance<self.delta).any() : #zero already exists
            return np.where(distance<self.delta)[0][0] , i
        else:                       # This will run if the value hasn't been found yet
            self.zeros=np.reshape(np.append(self.zeros,value),(-1,2))   #add value to zeros, resize index (-1,2)
            return self.zeros.size/2-1 , i

    def plot(self,N,a,b,c,d,which_newton=None):#Author: Ansgar
        X,Y=np.meshgrid(np.linspace(a,b,N),np.linspace(c,d,N),indexing='ij')#This is standard format for meshgrid
        vectorized=scipy.vectorize(self.find_zeros)         #This is done for ambigious values, this fixes it.
        A,B=vectorized(X,Y,which_newton)    #B is containg Both cordintes for and iterations to plot both position and color 
        plt.pcolor(A,cmap='viridis')
        plt.colorbar()
        plt.show()
        return A,B
    
    def plot2(self,N,a,b,c,d,which_newton=None):#Author: Ansgar
        """Plot function, taking six inputs. The resolution, x limits, y limits and which newton-ralphson method to use; Either False or something """
        X,Y=np.meshgrid(np.linspace(a,b,N),np.linspace(c,d,N),indexing='ij')
        vectorized=scipy.vectorize(self.find_zeros)#This is done for ambigious values, this fixes it.
        A,B=vectorized(X,Y,which_newton)#
        plt.pcolormesh(B,cmap='viridis')
        plt.colorbar()
        plt.show()
        return A,B


def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])
#print(fractal2d(f).partialdiv(0,0))#Print the partial derivitive of f at point (0,0)
#print(fractal2d(f).newtonmethod((20,10)))#Prints either a convering or divering number
#print(fractal2d(f).plot(100,-10,10,-10,10))
#fractal2d(f).plot(200,-1,1,-1,1)
def h(x,y):
    return np.array([x**8-28*x**6*y**2+70*x**4*y**4+15*x**4-28*x**2*y**6-90*x**2*y**2+y**8+15*y**4-16,8*x**7*y-56*x**5*y**3+56*x**3*y**5+60*x**3*y-8*x*y**7-60*x*y**3])
#fractal2d(h).plot(200,-1,1,-1,1)
def g(x,y):
    return np.array([x**3-3*x*y**2-2*x-2,3*x**2*y-y**3-2*y])

#fractal2d(h).plot2(1000,-10,10,-10,10)
#print(fractal2d(h).plot2(200,-10,10,-10,10,True))
fractal2d(g).plot(100,-10,10,-10,10)
fractal2d(g).plot2(100,-1,1,-1,1)