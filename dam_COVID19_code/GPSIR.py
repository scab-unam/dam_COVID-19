import matplotlib.pyplot as plt
import random

class Person(object):
    """
    A class used to represent an individual participating in the dynamics of the epidemy 

    ...

    Attributes
    ----------
    contacts: dict
        An empty dictionary in which the individual's contacts are recorded 
        with the contact time as the key and a list of the individuals with 
        whom there was contact as the value.
    
    tau : int
        The number of days that the individual will be infected
        (default 0).
    
    index : int
        The index with which the person is identified.
    
    status : dict
        A dictionary in which the individual's status at each time is recorded 
        with the time as the key and the population to which it belongs.
    
    age : int
        The age of the individual.
    
    comorbidities : list
        A list with the comorbidities of the individual 
        (default an empty list []).

    Methods
    -------
    set_tau():
        Sets the number of days that the infection will last depending
        on the comorbidities of the individual.
    
    set_status(t, s):
        Adds the status, s, of the individual at time, t, to the dictionary
        estatus.
        
    set_contact(t, person):
        It adds a contact with a person at the time t to the dictionary of 
        contacts.
        
    get_contacts():
        Prints the contacts of the individual at each time showing the time of
        the contact and a tuple containing the index of the person with whom 
        there was contact and the population to which it belongs.
    """
    
    def __init__(self, index, status, age, comorbidities = []):
        """
        Parameters
        ----------
        index : int
            The index with which the person is identified.
        
        status : dict
            A dictionary in which the individual's status at each time is recorded 
            with the time as the key and the population to which it belongs.
    
        age : int
            The age of the individual.
    
        comorbidities : list
            A list with the comorbidities of the individual 
            (default an empty list []).
        """
        
        self.contacts = {}
        self.tau = 0
        self.index = index
        self.status = status
        self.age = age
        self.comorbidities = comorbidities
        
    def set_tau(self, distribution, parameter):
        """Sets the number of days that the infection will last depending
        on the comorbidities of the individual.
        
        Parameters
        ----------
        distribution : function
            probability distribution of the length of the infectious state for an
            individual (tau).
        parameter : float
            parameter of the probability distribution of the length of the infectious 
            state for an individual (tau).
        """
        
        self.tau = distribution(parameter)

    def set_status(self, t, s):
        """Sets the status, s, of the individual at time, t, to the dictionary
        estatus.
        
        Parameters
        ----------
        t : int
            The time for which the status will be seted.
        s : str
            The population to which it the individual belongs at time t
            (S = susceptible, I = infected, R = removed).
        """
        
        self.status[t] = s

    def set_contact(self, t, person):
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
        person : Person
            The individual with whom there was contact.
        """

        if t in self.contacts:
            self.contacts[t].append(person)
        else:
            self.contacts[t] = [person]

    def get_contacts(self):
        """Prints the contacts of the individual at each time showing the time of
        the contact and a tuple containing the index of the person with whom 
        there was contact and the population to which it belongs."""
        
        print("\nContactos para la persona {}:".format(self.index))
        for t in self.contacts:
            s = str(t)+ " :"
            for c in self.contacts[t]:
                s += " ({},{})".format(c.index,c.status[t])
            print(s)

    def __eq__(self, other):
        """
        We overwrite the method __eq__ in orther to give sense to the boolean
        operator == in the class Persona. In orther to overwirte this method
        we recibe as parameters the object itself and other object to compare.
        The new coparison criterion is the index.
        """
        return (self.index == other.index)

def Bernoulli(p):
    """Simulates a random variable with Bernoulli distribution
    with parameter p.
        
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

def Geometric(p):
    """Simulates a random variable with Gemoetric distribution
    with probability of success p.
        
    Parameters
    ----------
    p : float
        The probability of success of a Bernoulli trial.

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

def Constant(c):
    """Simulates a random variable, X, with constant distribution.
        
    Parameters
    ----------
    c : float
        The constant such that P(X = c) = 1.

    Returns
    -------
    c : float
        The constant such that P(X = c) = 1.
        
    """
    return c

def Uniform(limits):
    """Simulates a random variable with Uniform distribution
    over the interval [limits[0],limits[1]].
        
    Parameters
    ----------
    limits : list
        List with the lower bownd of the interval where the random
        variable is defined as first item and the upper bownd
        as second item.

    Returns
    -------
    float
        A number between [limits[0],limits[1]].
    """
    return limits[0] + (limits[1]-limits[0])*random.random()

def GPSIR(S0, I0, R0, T, alpha, beta, tau_distribution = Geometric, 
          tau_parameter = 1/7, inoculation_distribution = Constant,
          inoculation_parameter = 1, inoculation_threshold = 1):
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
    
    tau_distribution : function, optional
        Distribution of the number of days that the individual 
        will be infected (default is Gemoetric).
    
    tau_parameter: float, optional
        Parameter of the distribution of the number of days that
        the individual will be infected (default is 1/7).
        
    inoculation_distribution: function, optional
        Distribution of the viral load to which a susceptible 
        individual is exposed in an infectious contact
        (default is Constant).
        
    inoculation_parameter: float, optional
        Parameter of the distribution of the viral load to 
        which a susceptible is exposed in an infectious contact
        (default is 1).
    
    inoculation_threshold: float, optional
        Viral load threshold for contagion (default is 1).
    

    Returns
    -------
    Population : list
        List with the individuals and their modified attributes after
        T days of contacts, infections and recoveries. 
        
    """

    Population = []

    # The susceptible are added to the population
    for i in range(1,S0+1):
        p = Person(i,{0:'S'},21)
        p.set_tau(tau_distribution,tau_parameter)
        Population.append(p)
    
    # The infected are added to the population
    for i in range(S0+1,S0+I0+1):
        p = Person(i,{0:'I'},21)
        p.set_tau(tau_distribution,tau_parameter)
        Population.append(p)

    # We iterate over the days of the simulation
    for t in range(1,T+1):

        # Iteration per day over each individual of the population
        for p1 in Population:

            # Iteration per individual over the rest of the population in 
            # orther to generate contats between individuals.
            for p2 in Population:
                if p1.status[t-1] == 'S' and p1 != p2 and Bernoulli(alpha):
                    p1.set_contact(t-1,p2)

            # Here starts the modification of the populations S, I, and R.

            # Modifications over the Suceptible population
            if p1.status[t-1] == 'S':
                
                # Case in which there were contacts
                if (t-1) in p1.contacts:
                    viral_load = 0
                    
                    # for each contact we check if there was efective infection
                    for contact in p1.contacts[t-1]:
                        if contact.status[t-1] == 'I' and Bernoulli(beta):
                            viral_load += inoculation_distribution(inoculation_parameter)
                    
                    # Case in which the viral load surpassed the inoculation threshold
                    if viral_load >= inoculation_threshold:
                        p1.set_status(t,'I')
                    
                    # Case in which the viral load did not surpassed the inoculation 
                    # threshold
                    else:
                        p1.set_status(t,'S')
                
                # Case in which there were no contacts
                else:
                    p1.set_status(t,'S')

            # Modifications over the Infected population
            elif p1.status[t-1] == 'I':
                
                # Case in which the days of infection aren´t over
                if p1.tau > 0:
                    p1.set_status(t,'I')
                    p1.tau -= 1
                
                # Case in which the days of infection are over
                else:
                    p1.set_status(t,'R')

            # Modifications over the Removed population
            else:
                p1.set_status(t,'R')
    
    return Population

def final_populations(Population):
    """Counts de number of susceptible, infected, and recovered individuals at the
    end of the simulation.
        
    This function receives as parameter the list of the population after the simulation 
    has been run. It counts the number of susceptible, infected, and recovered individuals 
    at the end of the simulation and prints them.
        
    Parameters
    ----------
    Population : list
        List with the individuals and their modified attributes after
        T days of contacts, infections and recoveries. 

    """
    Sf = 0; If = 0; Rf = 0;
    for person in Population:
        if person.status[365] == 'S':
            Sf += 1
        elif person.status[356] == 'I':
            If += 1
        else:
            Rf += 1
    print("Susceptible: {}\nInfectious: {}\nRemoved: {}".format(Sf,If,Rf))

def Plot_Dynamics(simulation, T):
    """Plots the evolution of the populations of susceptible,
    infected and removed individuals over time for a simulation.
        
    Parameters
    ----------
    simulation : list
        List of Persona objects obtained from the function Simulacion_Contagios. 

    T : int
        Number of days in the simulation
    """
        
    S = []; I = []; R = [];
    for t in range(T+1):
        St = 0; It = 0; Rt = 0;
        for person in simulation:
            if person.status[t] == 'S':
                St += 1
            elif person.status[t] == 'I':
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

    ## 1.1 SIMPLE MODEL
	
    # Results
    T = 365
    Population_11 = GPSIR(95,5,0,T,alpha = 0.1,beta = 0.03)

	# Final Number of Susceptible, Infected, and Removed Individuals
    print("\nFinal Number of Susceptible, Infected, and Removed Individuals"
         + "for the Simple Model\n")
    final_populations(Population_11)

	# Dynamics for Each Individual
    print("\nDynamics for Each Individual in the Simple Model\n")
    for person in Population_11:
        print("\nDynamics for the individual {}:".format(person.index))
        print(person.status)

	# Contacts for One Individual
    print("\nContacts for One Individual\n")
    Population_11[0].get_contacts()

	# Plot of The Simulation of the Simple Model
    Plot_Dynamics(Population_11, T)

    ## 1.2 GENERAL PROBABILITY DISTRIBUTION OF THE LENGHT OF THE
    ## INFECTIOUS STATE OF AN INDIVIDUAL

    # Results
    Population_12 = GPSIR(95, 5, 0, T, alpha = 0.1, beta = 0.03, tau_distribution = Geometric, tau_parameter = 1/10)

    # Final Number of Susceptible, Infected, and Removed Individuals
    print("\nFinal Number of Susceptible, Infected, and Removed Individuals"
         + "for the Model 1.2\n")
    final_populations(Population_12)

    # Plot of The Simulation of the Model 1.2
    Plot_Dynamics(Population_12, T)

    ## 1.4 EFFECTIVE INFECTION AS A THRESHOLD PROCESS

    # Results
    Population_14 = GPSIR(95,5,0,T,alpha = 0.1, beta = 0.03, inoculation_distribution = Uniform,
          inoculation_parameter = [1,3], inoculation_threshold = 2)

    # Final Number of Susceptible, Infected, and Removed Individuals
    print("\nFinal Number of Susceptible, Infected, and Removed Individuals"
         + "for the Model 1.4\n")
    final_populations(Population_14)

    # Plot of The Simulation of the Model 1.4
    Plot_Dynamics(Population_14, T)


if __name__ == '__main__':
	main()

