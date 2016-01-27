# for all our graphing needs
from pylab import *
#from numpy import *

h = 0.01

def printmenu():
    print("1 - analyze one radius of the Quantum Harmonic Oscillator")
    print("2 - analyze the Harmonic Oscillator with radii of 3, 3.5, 4, ..., and 6")
    print("3 - graph the first three energy levels of the Harmonic Oscillator")
    print("    against radius for l = 0, 1, and 2")
    print("4 - analyze the change in the parameters of the quadratic fit for energy")
    print("    vs. n as R changes")
    print("5 - terminate the program")
    print("m - print this menu")

#takes the numerical approximation of the second order ODE 
def numapprox(r, L, e):
    y = [0.0, h**((L+1)*2)]
    u = h**(L+1)
    ubef = 0.0
    utemp = 0.0
    z = h
    A = 0.0
    while z < r:
        utemp = u
        u = (2.0+(h**2.0)*(L*(L+1)/(z**2.0)+(z**2.0)/4.0-e))*u-ubef
        ubef = utemp
        y.append(u**2)
        z = z + h
    return y

#gets the fit of a function with a parabolic shape
def quadfit(x = [], fx = [], *arg):
    (A, B, C) = polyfit(x, fx, 2)
    print("")
    print("A")
    print(A)
    print("B")
    print(B)
    print("C")
    print(C)
    
    i = 0
    xfit = []
    while(i < x[-1]):
        xfit.append(i)
        i = i + 0.1
    fit = polyval((A, B, C), xfit)
    return (xfit, fit)

#gets the fit of a function with in the form of a power law
def powfit(n, x = [], fx = [], *arg):
    xtem = map(lambda a: log(a - n), fx)
    (A, B) = polyfit(x, xtem, 1)
    print("")
    print("A")
    print(A)
    print("B")
    print(B)

    fit = map(lambda a: A/a, xtem)
    return fit


#makes the x axis of the energy vs. n graph
def xaxisgen(func = [], *arg):
    x = []
    i = 0
    
    size = len(func)
    while i < size:
        x.append(2*i)
        i = i + 1

    return x


#finds solutions
def solfinder(R, L, epmax, zar = [], *arg):
    usqra = []
    usqrb = [1]
    elevel = []
    einc = 0.0
    etemp = 50.0
    n = 0.0
    tem1 = 0.0
    tem1t = 0.0
    tem2 = 0.0
    tem3 = 0.0
    
    
    figure("R = " + str(R) + " & l = " + str(L))
    
    # find solutions
    while einc < (epmax + 0.1):
        usqra = numapprox(R, L, einc)
        tem1t = usqra[-1]
        
        #the energy values which are closest to 0 at R are the solutions 
        if((usqrb[-1] < usqra[-1]) and (usqrb[-1] < tem2)):
           #if the energy value has passed this if statement then it is a solution
           #and should be normalized.
           A = sum(usqrb)*h
           usqrb = map(lambda x: x/A, usqrb)
           n = n+1
           plot(zar, usqrb, label = ("n = 2*" + str(n) + ""))
           elevel.append(einc - 0.1)

        tem2 = tem1
        tem1 = tem1t
        del usqrb
        usqrb = usqra
        del usqra
        einc = einc+0.1

    #add labels and title to graph
    title("R = " + str(R) + " & l = " + str(L))
    xlabel("z")
    ylabel("Normalized Wavefunction")
            
    #add legend
    legend(loc="upper left")
        
    #save figure
    savefig("IMP_Proj_Wavefunc_Graph_l=" + str(L) + "_R=" + str(R) + ".pdf")
        
    return elevel


#finds n = 0, 2, and 4
def ft3n(R, L):
    usqra = []
    usqrb = [1]
    einc = 0.0
    tem1 = 0.0
    tem1t = 0.0
    tem2 = 0.0
    tem3 = 0.0
    soln = 1.0
    sol1 = 0.0
    sol2 = 0.0
    sol3 = 0.0
    
    # find solutions
    while soln < 4:
        usqra = numapprox(R, L, einc)
        tem1t = usqra[-1]
        
        #the energy values which are closest to 0 at R are the solutions 
        if((usqrb[-1] < usqra[-1]) and (usqrb[-1] < tem2)):                
            if(soln == 1):
                sol1 = einc - 0.1
            if(soln == 2):
                sol2 = einc - 0.1
            if(soln == 3):
                sol3 = einc - 0.1
            soln = soln + 1
            
        einc = einc+0.1
        tem2 = tem1
        tem1 = tem1t
        del usqrb
        usqrb = usqra
        del usqra
               
    return (sol1, sol2, sol3)


#solfinder for fit vs. radius option
def fitvRsolfinder(R, L, epmax):
    usqra = []
    usqrb = [1]
    nlist = []
    elevel = []
    einc = 0.0
    etemp = 50.0
    n = 0.0
    tem1 = 0.0
    tem1t = 0.0
    tem2 = 0.0
    tem3 = 0.0
    
    # find solutions
    while einc < (epmax + 0.05):
        usqra = numapprox(R, L, einc)
        tem1t = usqra[-1]
        
        #the energy values which are closest to 0 at R are the solutions 
        if((usqrb[-1] < usqra[-1]) and (usqrb[-1] < tem2)):
           #if the energy value has passed this if statement then it is a solution)
           nlist.append(n)
           elevel.append(einc - 0.1)
           n = n + 2

        tem2 = tem1
        tem1 = tem1t
        del usqrb
        usqrb = usqra
        del usqra
        einc = einc+0.05
        
    return nlist, elevel


def Ranalysis(Ra, maxe):
    print("")
    print("")
    print("")
    print("Beginning analysis of R = " + str(Ra) + ".")
    print("")
    print("")
    print("")
    
    elevl0 = []
    nl0 = []
    
    elevl1 = []
    nl1 = []
    
    elevl2 = []
    nl2 = []
    
    zar = [0, h]
    Z = h
    
    while Z < Ra:
        zar.append(Z)
        Z = Z + h    
    
    elevl0 = solfinder(Ra, 0, maxe, zar)
    
    elevl1 = solfinder(Ra, 1, maxe, zar)
    
    elevl2 = solfinder(Ra, 2, maxe, zar)
    

    figure("Energy Levels vs l for R = " + str(Ra))

    print("Each of the following functions will be fit to a quadratic.")
    print("A, B, and C are the constants in the equation Ax^2+Bx+C for")
    print("each fit.")
    nl0 = xaxisgen(elevl0)
    scatter(nl0, elevl0, label = "l = 0", color = "b")
    print("l = 0")
    (nl0, elevl0) = quadfit(nl0, elevl0)
    plot(nl0, elevl0, label = "Fit for l = 0", color = "b") 
    
    nl1 = xaxisgen(elevl1)
    scatter(nl1, elevl1, label = "l = 1", color = "g")
    print("")
    print("l = 1")
    (nl1, elevl1) = quadfit(nl1, elevl1)
    plot(nl1, elevl1, label = "Fit for l = 1", color = "g")
    
    nl2 = xaxisgen(elevl2)
    scatter(nl2, elevl2, label = "l = 2", color = "r")
    print("")
    print("l = 2")
    (nl2, elevl2) = quadfit(nl2, elevl2)
    plot(nl2, elevl2, label = "Fit for l = 2", color = "r")
    
    
    #add in labels and title
    xlabel('n')
    ylabel('Unitless Energy Number')
    title("Energy Level vs. n for R = " + str(Ra))
    
    #add legend
    legend(loc="lower right")
    
    #save figure
    savefig("IMP_Proj_Energies_of_solutions_R=" + str(Ra) + ".pdf")


def EvRfindsol(Ra):
    (sol1l0, sol2l0, sol3l0) = ft3n(Ra, 0)
    (sol1l1, sol2l1, sol3l1) = ft3n(Ra, 1)
    (sol1l2, sol2l2, sol3l2) = ft3n(Ra, 2)

    return (sol1l0, sol2l0, sol3l0, sol1l1, sol2l1, sol3l1, sol1l2, sol2l2, sol3l2)

def EvR(maxR):
    l0gRs1 = []
    l0gRs2 = []
    l0gRs3 = []
    l1gRs1 = []
    l1gRs2 = []
    l1gRs3 = []
    l2gRs1 = []
    l2gRs2 = []
    l2gRs3 = []
    xl = []
    Rinc = 0.1
    cR = 0.5

    #make lists of 
    while(cR <= maxR):
        (l0s1, l0s2, l0s3, l1s1, l1s2, l1s3,  l2s1, l2s2, l2s3) = EvRfindsol(cR)
        l0gRs1.append(l0s1)
        l0gRs2.append(l0s2)
        l0gRs3.append(l0s3)
        l1gRs1.append(l1s1)
        l1gRs2.append(l1s2)
        l1gRs3.append(l1s3)
        l2gRs1.append(l2s1)
        l2gRs2.append(l2s2)
        l2gRs3.append(l2s3)
        xl.append(cR)
        cR = cR + Rinc
    
    #make figure and graphs of energy vs. R for l = 0
    figure("n = 0, 1, and 2 for l = 0 vs. R")
    scatter(xl, l0gRs1, label = "n = 0", color = "b")
    aa = l0gRs1
    Rfitl0s1 = powfit(3/2+0+0, xl, aa)
    plot(xl, Rfitl0s1, label = "fit n = 0", color = "b")

    scatter(xl, l0gRs2, label = "n = 2", color = "g")
    aa = l0gRs2
    Rfitl0s2 = powfit(3/2+1+0, xl, aa)
    plot(xl, Rfitl0s2, label = "fit n = 2", color = "g")
    
    scatter(xl, l0gRs3, label = "n = 4", color = "r")
    aa  = l0gRs3
    Rfitl0s3 = powfit(3/2+2+0, xl, aa)
    plot(xl, Rfitl0s3, label = "fit n = 4", color = "b")
    
    #add in labels and title
    xlabel('R')
    ylabel('Unitless Energy Number')
    title("n = 0, 2, and 4 for l = 0 vs. R")

    #add legend
    legend(loc="upper right")

    #save figure
    savefig("IMP_Proj_Energies_of_Solutions_vs_R_l=0.pdf")

    #make figure and graphs of energy vs. R for l = 1
    figure("n = 0, 2, and 4 for l = 1 vs. R")
    scatter(xl, l1gRs1, label = "n = 0", color = "b")
    aa  = l1gRs1
    Rfitl1s1 = powfit(3/2+0+1, xl, aa)
    plot(xl, Rfitl1s1, label = "fit n = 0", color = "b")

    scatter(xl, l1gRs2, label = "n = 2", color = "g")
    aa  = l1gRs2
    Rfitl1s2 = powfit(3/2+1+1, xl, aa)
    plot(xl, Rfitl1s2, label = "fit n = 0", color = "b")

    scatter(xl, l1gRs3, label = "n = 4", color = "r")
    aa  = l1gRs3
    Rfitl1s3 = powfit(3/2+2+1, xl, aa)
    plot(xl, Rfitl1s3, label = "fit n = 0", color = "b")

    #add in labels and title
    xlabel('R')
    ylabel('Unitless Energy Number')
    title("n = 0, 2, and 4 for l = 1 vs. R")

    #add legend
    legend(loc="upper right")

    #save figure
    savefig("IMP_Proj_Energies_of_Solutions_vs_R_l=1.pdf")

    #make figure and graphs of energy vs. R for l = 2
    figure("n = 0, 2, and 4 for l = 2 vs. R")
    scatter(xl, l2gRs1, label = "n = 0", color = "b")
    aa  = l2gRs1
    Rfitl2s1 = powfit(3/2+0+2, xl, aa)
    plot(xl, Rfitl2s1, label = "fit n = 0", color = "b")

    scatter(xl, l2gRs2, label = "n = 2", color = "g")
    aa  = l2gRs2
    Rfitl2s2 = powfit(3/2+1+2, xl, aa)
    plot(xl, Rfitl2s2, label = "fit ", color = "g")

    scatter(xl, l2gRs3, label = "n = 4", color = "r")
    aa  = l2gRs3
    Rfitl2s3 = powfit(3/2+2+2, xl, aa)
    plot(xl, Rfitl2s3, label = "fit n = 0", color = "r")
    
    #add in labels and title
    xlabel('R')
    ylabel('Unitless Energy Number')
    title("n = 0, 2, and 4 for l = 2 vs. R")

    #add legend
    legend(loc="upper right")

    #save figure
    savefig("IMP_Proj_Energies_of_Solutions_vs_R_l=2.pdf")


def fitvR(maxR, minR, epmax):
    R = minR
    Al0list = []
    Al1list = []
    Al2list = []
    Bl0list = []
    Bl1list = []
    Bl2list = []
    Cl0list = []
    Cl1list = []
    Cl2list = []
    Rlist = []
    while(R <= maxR):
        (nl0, el0) = fitvRsolfinder(R, 0, epmax)
        (nl1, el1) = fitvRsolfinder(R, 1, epmax)
        (nl2, el2) = fitvRsolfinder(R, 2, epmax)
        (Al0, Bl0, Cl0) = polyfit(nl0, el0, 2)
        (Al1, Bl1, Cl1) = polyfit(nl1, el1, 2)
        (Al2, Bl2, Cl2) = polyfit(nl2, el2, 2)
        Al0list.append(Al0)
        Al1list.append(Al1)
        Al2list.append(Al2)
        Bl0list.append(Bl0)
        Bl1list.append(Bl1)
        Bl2list.append(Bl2)
        Cl0list.append(Cl0)
        Cl1list.append(Cl1)
        Cl2list.append(Cl2)
        Rlist.append(R)
        R = R + 0.1

    #make figure and graphs
    figure("The A coefficient in the quadratic fit of energy vs n for l = 0, 1, and 2 vs. R")
    scatter(Rlist, Al0list, label = "l = 0", color = "b")
    plot(Rlist, Al0list, color = "b")
    scatter(Rlist, Al1list, label = "l = 1", color = "g")
    plot(Rlist, Al1list, color = "g")
    scatter(Rlist, Al2list, label = "l = 2", color = "r")
    plot(Rlist, Al2list, color = "r")

    #add in labels and title
    xlabel('R')
    ylabel('Value of A')
    title("The A coefficient in the quadratic fit of energy vs n for l = 0, 1, and 2 vs. R")

    #add legend
    legend(loc="upper right")

    #save figure
    savefig("IMP_Proj_quad_fit_param_of_energy_vs_n_param_A_vs_R.pdf")

    #make figure and graphs
    figure("The B coefficient in the quadratic fit of energy vs n for l = 0, 1, and 2 vs. R")
    scatter(Rlist, Bl0list, label = "l = 0", color = "b")
    plot(Rlist, Bl0list, color = "b")
    scatter(Rlist, Bl1list, label = "l = 1", color = "g")
    plot(Rlist, Bl1list, color = "g")
    scatter(Rlist, Bl2list, label = "l = 2", color = "r")
    plot(Rlist, Bl2list, color = "r")

    #add in labels and title
    xlabel('R')
    ylabel('Value of B')
    title("The B coefficient in the quadratic fit of energy vs n for l = 0, 1, and 2 vs. R")

    #add legend
    legend(loc="upper right")

    #save figure
    savefig("IMP_Proj_quad_fit_param_of_energy_vs_n_param_B_vs_R.pdf")
    
    #make figure and graphs
    figure("The C coefficient in the quadratic fit of energy vs n for l = 0, 1, and 2 vs. R")
    scatter(Rlist, Cl0list, label = "l = 0", color = "b")
    plot(Rlist, Cl0list, color = "b")
    scatter(Rlist, Cl1list, label = "l = 1", color = "g")
    plot(Rlist, Cl1list, color = "g")
    scatter(Rlist, Cl2list, label = "l = 2", color = "r")
    plot(Rlist, Cl2list, color = "r")

    #add in labels and title
    xlabel('R')
    ylabel('Value of C')
    title("The C coefficient in the quadratic fit of energy vs n for l = 0, 1, and 2 vs. R")

    #add legend
    legend(loc="upper right")

    #save figure
    savefig("IMP_Proj_quad_fit_param_of_energy_vs_n_param_C_vs_R.pdf")

    


def oneanalyse():
    print("Please choose the size of your radius.")
    sz = input("-->")
    print("Please choose the maximum energy value")
    print("that the program will search.")
    maxe = input("-->")
    Ranalysis(sz, maxe)
    show()


def oursixanalyses():
    R = 3
    while(R <= 6):
        Ranalysis(R, 30)
        R = R + 0.5
    show()

          
def EvRanalysis():
    print("Please enter the maximum radius you want the program to graph.")
    maxR = input("-->")
    EvR(maxR)
    show()


def fitanalysis():
    print("Please enter the minimum radius you want the program to graph.")
    minR = input("-->")
    print("Please enter the maximum radius you want the program to graph.")
    maxR = input("-->")
    print("Please enter the maximum energy you want the program to search.")
    epmax = input("-->")
    fitvR(maxR, minR, epmax)
    show()


def terminate():
    exit()

print("Please, make a selection.\n")
printmenu()    
opt = {'1' : oneanalyse,
       '2' : oursixanalyses,
       '3' : EvRanalysis,
       '4' : fitanalysis,
       '5' : terminate,
       'M' : printmenu,
       'm' : printmenu
}

while (1):
    ch = raw_input("-->")
    opt[ch]()
