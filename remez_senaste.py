from numpy import *
from matplotlib.pyplot import *
from numpy.polynomial import Polynomial
from scipy import optimize

n=10  #graden på polynomet
a=2 #startvärde
b=5  #slutvärde

x=linspace(a,b,n+2) # vårt intervall mellan a och b, n+2 punkter
t=linspace(a,b,int(1e4))# xaxel med tätare xvärden för att testa funktionen

def f(x): #definera en testfunktion
    return cos(x)
def fprim(x):
    delta=1e-7
    return (f(x+delta)-f(x))/delta

def createPolynom(x):
    M=zeros((n+2,n+2)) #Skapa en n+2 x n+2 matris fylld med nollo
    for i in range(n+2): #Första kolumnen bara ettor? Fyll matrisen rad för rad, förutom sista kolumnen.
            for j in range(n+1):
                M[i][j]=x[i]**j

    for i in range(n+2): #Sista kolumnen, fylld med 1 och -1.
        M[i][-1]=(-1)**i
    F=f(x) #Kända funktionsvärden för testfunktionen i värdena x0,x1,x2 ... i vårt intervall.

    global Pcoeff
    #Pcoeff= linalg.inv(M)@(F) #Löser linjära ekvationsystemet. Får ut koefficienterna för polynomet, där sista elementet är h. 
    Pcoeff= linalg.solve(M,F) 
    global p
    p=Polynomial(Pcoeff[:-1]) #Polynomial skapar ett polynom objekt med värdena i Pcoeff förutom sista värdet 'h'.
    return p
createPolynom(x)

x=linspace(a,b,n+2)
counter=0
tol=1e-9
while abs(Pcoeff[-1]-max(abs(p(t)-f(t))))>tol and counter<50:
    counter+=1
    createPolynom(x) #Skapar nytt polynom med efter bytt ut en punkt i referensen.
    
    L = [(abs(f(i)-p(i))) for i in t] #Används för att hitta ett max i ett tätare intervall istället för derivera. 
    max_value = max(L)
    my = t[L.index(max_value)]
    
    #my=optimize.newton(lambda x:abs((p.deriv())(x)-fprim(x)),) #Alternativ metod för my? 
    #Behöver derivera testfunktionen. Gissa nollställe som maxvärdet?
    
#Steg 4.
#
     
    
    if my<x[0]:
        if sign(f(my)-p(my)) == sign(f(x[0])-p(x[0])):
            x[0]=my
        else: 
            x[-1]=my
    
    if my>x[-1]:
        if sign(f(my)-p(my)) == sign(f(x[-1])-p(x[-1])):
            x[-1]=my
        else:
            x[0]=my
    else:
        for i in range(len(x)):
            if (my<=x[i]) :
                x_i=i-1
                break
        if sign(f(my)-p(my)) == sign(f(x[x_i])-p(x[x_i])):
            x[x_i] = my
        else:
            x[x_i+1] = my 
    
    #x.sort()
    
    print(max(abs(p(t)-f(t))))
    #print(f'Opposite signs:{sign(f(x[0])-p(x[0]))!=sign(f(x[1])-p(x[1]))}') 
    #print(f'h={abs(Pcoeff[-1])}')
    #print(f'my={my}')

for i in range(len(x)):
    if (my<=x[i]) :
        print(f'i={i}')
        x_i=i
        break
if sign(f(my)-p(my)) == sign(f(x[x_i])-p(x[x_i])):
    x[x_i] = my
else:
    x[x_i+1] = my
print(x)

L = [(abs(f(i)-p(i))) for i in t] #Används för att hitta ett max i ett tätare intervall istället för derivera. 
max_value = max(L)
my = t[L.index(max_value)]
my

try:
    createPolynom(x)
except:    
    print(x)
L = [(abs(f(i)-p(i))) for i in t] #Används för att hitta ett max i ett tätare intervall istället för derivera. 
max_value = max(L)
my = t[L.index(max_value)]

if my<x[0]:
        if sign(f(my)-p(my)) == sign(f(x[0])-p(x[0])):
            x[0]=my
        else: 
            x[-1]=my
    
if my>x[-1]:
    if sign(f(my)-p(my)) == sign(f(x[-1])-p(x[-1])):
        x[-1]=my
    else:
        x[0]=my
else:
    for i in range(len(x)):
        if (my<x[i]) :
            x_i=i-1
            break
    if sign(f(my)-p(my)) == sign(f(x[x_i])-p(x[x_i])):
        x[x_i] = my
    else:
        x[x_i+1] = my

t=linspace(a,10,int(1e4))
figure(0)
plot(t,p(t)) #Plottar f och p i samma graf som jämförelse
plot(t,f(t))
figure(1)
plot(t,(abs(p(t)-f(t)))) #Plottar error funktionen.
yscale('log')
plot(x,1+0*x,'*')

#print(x)
print(f'h={abs(Pcoeff[-1])}')
print(f'my={my}')
print(max(abs(p(t)-f(t))))

createPolynom(x)
L = [(abs(f(i)-p(i))) for i in t] #Används för att hitta ett max i ett tätare intervall istället för derivera. 
max_value = max(L)
my = t[L.index(max_value)]
figure(0)
plot(t,p(t)) #Plottar f och p i samma graf som jämförelse
plot(t,f(t))
figure(1)
plot(t,abs(p(t)-f(t))) #Plottar error funktionen.
print(f'h={abs(Pcoeff[-1])}')
print(f'my={my}')
print(max(abs(p(t)-f(t))))

#if (x[0]<=my) and (my<=x[-1]):
 #       for i in range(len(x)):
  #          if my<x[i]:
   #             x_i = i-1
    #            if x_i == -1:
     #               x_i = 0
      #  if sign(f(my)-p(my)) == sign(f(x[x_i])-p(x[x_i])):
       #     x[x_i] = my
        #else:
         #   x[x_i + 1] = my

class RemezApprox:
    def __init__(self,f):
        self.func=f
    def remez(self,n,ref_a,ref_b):
        a=ref_a
        b=ref_b
        self.degree=n
        x=array(linspace(a,b,n+2))
        
        M=zeros((n+2,n+2))
        for i in range(n+2): 
            for j in range(n+1):
                M[i][j]=x[i]**j
        for i in range(n+2):
            M[i][-1]=(-1)**i
        F=self.func(x)
        Pcoeff= linalg.inv(X)@(F)
        self.p=Polynomial(Pcoeff[:-1])
    def plot(self,xval):
        figure(0)
        plot(xval,self.p(xval))
        plot(xval,self.func(xval))

abs(Pcoeff[-1]-max(abs(p(t)-f(t))))