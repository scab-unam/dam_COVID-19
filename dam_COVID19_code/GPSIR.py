import matplotlib.pyplot as plt
import random

class Persona(object):
    """
    A class used to represent an individual participating in the dynamics 
    of the epidemy 

    ...

    Attributes
    ----------
    contactos: dict
        An empty dictionary in which the individual's contacts are recorded 
        with the contact time as the key and a list of the individuals with 
        whom there was contact as the value.
    
    tao : int
        The number of days that te individual will be infected
        (default 0).
    
    indice : int
        The index with which the person is identified.
    
    estatus : dict
        A dictionary in which the individual's status at each time is recorded 
        with the time as the key and the population to which it belongs.
    
    edad : int
        The age of the individual.
    
    comorbilidades : list
        A list with the comorbidities of the individual 
        (default an empty list []).

    Methods
    -------
    establece_tao():
        Sets the number of days that the infection will last depending
        on the comorbidities of the individual.
    
    agrega_estatus(t, s):
        Adds the status, s, of the individual at time, t, to the dictionary
        estatus.
        
    agrega_contacto(t, person):
        It adds a contact with a person at the time t to the dictionary of 
        contacts.
        
    muestra_contactos():
        Prints the contacts of the individual at each time showing the time of
        the contact and a tuple containing the index of the person with whom 
        there was contact and the population to which it belongs.
    """
    
    def __init__(self, indice, estatus, edad, comorbilidades = []):
        """
        Parameters
        ----------
        indice : int
            The index with which the person is identified.
        
        estatus : dict
            A dictionary in which the individual's status at each time is recorded 
            with the time as the key and the population to which it belongs.
    
        edad : int
            The age of the individual.
    
        comorbilidades : list
            A list with the comorbidities of the individual 
            (default an empty list []).
        """
        
        self.contactos = {}
        self.tao = 0
        self.indice = indice
        self.estatus = estatus
        self.edad = edad
        self.comorbilidades = comorbilidades
        
    def establece_tao(self):
        """Sets the number of days that the infection will last depending
        on the comorbidities of the individual."""
        
        gamma = 1/7
        self.tao = Geometrica(gamma)

    def agrega_estatus(self, t, s):
        """Adds the status, s, of the individual at time, t, to the dictionary
        estatus.
        
        Parameters
        ----------
        t : int
            The time for which the status will be seted.
        s : str
            The population to which it the individual belongs at time t
            (S = susceptible, I = infected, R = removed).
        """
        
        self.estatus[t] = s

    def agrega_contacto(self, t, persona):
        """It adds a contact with a person at the time t to the dictionary of 
        contacts.
        
        This method receives as parameters the time in which the contact was made 
        and the person with whom it was made. If there has not been a contact in 
        that time, it creates a new time key and assigns it a list with the person 
        with whom there was contact. In case there was already a contact at that time, 
        the person is added to the contacts list.
        
        Parameters
        ----------
        t : int
            The time in which the contact happened.
        persona : Persona
            The individual with whom there was contact.
        """

        if t in self.contactos:
            self.contactos[t].append(persona)
        else:
            self.contactos[t] = [persona]

    def muestra_contactos(self):
        """Prints the contacts of the individual at each time showing the time of
        the contact and a tuple containing the index of the person with whom 
        there was contact and the population to which it belongs."""
        
        print("\nContactos para la persona {}:".format(self.indice))
        for t in self.contactos:
            s = str(t)+ " :"
            for c in self.contactos[t]:
                s += " ({},{})".format(c.indice,c.estatus[t])
            print(s)

    def __eq__(self, other):
        """
        We overwrite the method __eq__ in orther to give sense to the boolean
        operator == in the class Persona. In orther to overwirte this method
        we recibe as parameters the object itself and other object to compare.
        The new coparison criterion is the index.
        """
        return (self.indice == other.indice)

def Bernoulli(p):
    """Simulates a random variable with Bernoulli distribution with 
    parameter p.
        
    Parameters
    ----------
    p : float
        The parameter of the Bernoulli distribution.

    Returns
    -------
    int
        0 or 1 depending of the result of the experiment
        
    """
    
    u = random.random()
    if u < p:
        return 1
    else:
        return 0

def Geometrica(p):
    """Simulates a random variable with Gemoetric distribution with 
    parameter p.
        
    Parameters
    ----------
    p : float
        The parameter of the Geometric distribution.

    Returns
    -------
    L : int
        The number L of Bernoulli trials needed to get one success.
        
    """
    L = 0
    aux = 0
    while aux == 0:
        if Bernoulli(p) == 1:
            aux = 1
        L += 1
    return L

def Simulacion_Contagios(S0,I0,R0,T,alpha,beta):
    """Simulates the dynamics of the epidemy
        
    Parameters
    ----------
    S0 : int
        Initial population of susceptible individuals.
    
    I0 : int
        Initial population of infected individuals.
    
    R0 : int
        Initial population of recovered individuals.
    
    T : int
        Number of days of the simulation.
    
    alpha: float
        Probability of formation of an edge.
    
    beta: float
        Probability of contagion given infectious contagion.

    Returns
    -------
    Poblacion : list
        List with the individuals and their modified attributes after
        T days of contacts, infections and recoveries. 
        
    """

    Poblacion = []

    # The susceptible are added to the population
    for i in range(1,S0+1):
        p = Persona(i,{0:'S'},21)
        p.establece_tao()
        Poblacion.append(p)
    
    # The infected are added to the population
    for i in range(S0+1,S0+I0+1):
        p = Persona(i,{0:'I'},21)
        p.establece_tao()
        Poblacion.append(p)

    # We iterate over the days of the epidemy
    for t in range(1,T+1):

        # Iteration per day over each individual of the population
        for p1 in Poblacion:

            # Iteration per individual over the rest of the population in orther
            # to generate contats between individuals.
            for p2 in Poblacion:
                if p1.estatus[t-1] == 'S' and p1 != p2 and Bernoulli(alpha):
                    p1.agrega_contacto(t-1,p2)

            # Here starts the modification of the populations S, I, and R.

            # Modifications over the Suceptible population
            if p1.estatus[t-1] == 'S':
                
                # Case in which there were contacts
                if (t-1) in p1.contactos:
                    temp = 0
                    
                    # for each contact we check if there was efective infection
                    for contacto in p1.contactos[t-1]:
                        if contacto.estatus[t-1] == 'I' and Bernoulli(beta):
                            temp +=1
                    
                    # Case in which there was at least one contagious contact
                    if temp >= 1:
                        p1.agrega_estatus(t,'I')
                    
                    # Case in which there were no contagious contacts
                    else:
                        p1.agrega_estatus(t,'S')
                
                # Case in which there were no contacts
                else:
                    p1.agrega_estatus(t,'S')

            # Modifications over the Infected population
            elif p1.estatus[t-1] == 'I':
                
                if p1.tao > 0:
                    p1.agrega_estatus(t,'I')
                    p1.tao -= 1
                
                else:
                    p1.agrega_estatus(t,'R')

            # Modifications over the Removed population
            else:
                p1.agrega_estatus(t,'R')
    
    return Poblacion

def Grafica_Dinamicas(simulacion, T):
    """Plots the evolution of the populations of susceptible,
    infected and removed individuals over time for a simulation.
        
    Parameters
    ----------
    simulacion : list
        List of Persona objects obtained from the function Simulacion_Contagios. 

    T : int
        Number of days in the simulation
    """
        
    S = []; I = []; R = [];
    for t in range(T+1):
        St = 0; It = 0; Rt = 0;
        for persona in simulacion:
            if persona.estatus[t] == 'S':
                St += 1
            elif persona.estatus[t] == 'I':
                It += 1
            else:
                Rt += 1
        S.append(St)
        I.append(It)
        R.append(Rt)

    plt.plot(range(T+1), S, 'm', label = 'S')
    plt.plot(range(T+1), I, 'b', label = 'I')
    plt.plot(range(T+1), R, 'g', label = 'R')
    plt.legend()
    plt.xlabel('Time (days)')
    plt.ylabel('Population Size')
    plt.title('')
    plt.show()

def main():

	# Results
	T = 365
	Poblacion = Simulacion_Contagios(95,5,0,T,alpha = 0.1,beta = 0.02)

	# Final Number of Susceptible, Infected, and Removed Individuals
	Sf = 0; If = 0; Rf = 0;
	for persona in Poblacion:
		if persona.estatus[T] == 'S':
			Sf += 1
		elif persona.estatus[T] == 'I':
			If += 1
		else:
			Rf += 1
	print("Susceptibles: {}\nInfectados: {}\nRecuperados: {}".format(Sf,If,Rf))

	# Dynamics for Each Individual
	for persona in Poblacion:
		print("\nDinámica para la persona {}:".format(persona.indice))
		print(persona.estatus)

	# Contacts for One Individual
	Poblacion[0].muestra_contactos()

	# Plot of The Simulation
	Grafica_Dinamicas(Poblacion, T)

if __name__ == '__main__':
	main()

