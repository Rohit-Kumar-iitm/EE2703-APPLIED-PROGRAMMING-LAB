"""
        EE2703: Applied Programming Lab 
        Assignment-2
        Spice Part-2
        Name: Rohit Kumar 
        Roll no.: EE20B111
        Date: 11-02-2022
"""
#Importing required libraries
from re import I
import numpy as np #To solve the matrix
from sys import argv, exit #To pass the arguments(2 in number)
import math #Using math library for sine and cosine function and constant pi

#Making classes for various components to store the data concisely and precisely

#For voltage source
class VoltageSource:
    def __init__(self, ElementName, N1, N2, Val):
        self.ElementName = ElementName
        self.N1 = N1 #N1 == Node 1 
        self.N2 = N2 #N2 == Node 2
        self.Val = Val #Val == Value of the element

#For current source
class CurrentSource:
    def __init__(self, ElementName, N1, N2, Val):
        self.ElementName = ElementName
        self.N1 = N1 #N1 == Node 1 
        self.N2 = N2 #N2 == Node 2
        self.Val = Val #Val == Value of the element

#For resistor
class Resistor:
    def __init__(self, ElementName, N1, N2, Val):
        self.ElementName = ElementName
        self.N1 = N1 #N1 == Node 1 
        self.N2 = N2 #N2 == Node 2
        self.Val = Val #Val == Value of the element               

#For inductor
class Inductor:
    def __init__(self, ElementName, N1, N2, Val):
        self.ElementName = ElementName
        self.N1 = N1 #N1 == Node 1 
        self.N2 = N2 #N2 == Node 2
        self.Val = Val #Val == Value of the element

#For capacitor
class Capacitor:
    def __init__(self, ElementName, N1, N2, Val):
        self.ElementName = ElementName
        self.N1 = N1 #N1 == Node 1 
        self.N2 = N2 #N2 == Node 2
        self.Val = Val #Val == Value of the element       

if len(argv) == 2: #If two files are given(.py and .netlist), then only the loop proceeds
    if argv[1].split('.')[1] != 'netlist': #Checking the type of file extension the second argument(:= argv[1]) is having, it should be a 'netlist' file
        print("Error. Incorrect type of input file is given.\nExpected '.netlist' type of file.")
    else: #Start solving if correct files are obtained and given

        CIRCUIT_START = '.circuit' #To check for beginning of the circuit
        CIRCUIT_END = '.end' #To check for end of the circuit
        CIRCUIT_AC = '.ac' #To check for ac component in the circuit 
        INDEX_AC = 0 #It stores the index of the line containing '.ac' , and is initialised to zero
        INDEX_CIRCUIT = -1 #It stores the index of where the '.circuit' is present in the given file and initialised to -1
        INDEX_END = -2 #It stores the index of where the '.end' is present in the given file and initialised to -2
        TOKENS = [] #It stores the tokens from the netlist file
        NODE_COUNT = 0 #Keeps track of total number of Nodes! And is initialised to zero
        VoltageSource_Count = 0 #Keeps the track of total number of voltage sources, and is initialised to zero 
        VoltageTemp = 0 ####
        
        with open(argv[1]) as netlist: #Opening the netlist file to extract the data in it, naming it as 'netlist' file
            lines = netlist.readlines() #Extracting each and every line from the 'netlist' and storing it in a list by name 'lines'
            #Now each element of the list stores each line(0th element = 0th line of netlist etc.) which was extracted from the 'netlist' file`
            
            #Below "for" loop finds the index of the line where '.circuit', '.end', '.ac' are present and also the frequency of the ac source
            for i in range(0, len(lines)):
                
                if lines[i][:len(CIRCUIT_START)] == CIRCUIT_START: #Comparing the string elements of the line(each line of lines) with CIRCUIT_START(i.e., '.circuit') 
                    INDEX_CIRCUIT = i #Assign the index of the line to the variable where '.circuit' is found
                
                elif lines[i][:len(CIRCUIT_END)] == CIRCUIT_END: #Comparing the string elements of the line(each line of lines) with CIRCUIT_END(i.e., '.end')
                    INDEX_END = i #Assign the index of the line to the variable where '.end' is found
                
                elif lines[i][:len(CIRCUIT_AC)] == CIRCUIT_AC: #Comparing the string elements of the line(each line of lines) with CIRCUIT_AC(i.e., '.ac')
                    INDEX_AC += i #It is incrementing(being added) because there may be multiple ac sources!   
                    TOKENS_AC = lines[i].split() #Extracting the tokens of the line containing '.ac' 
                    OMEGA = 2*math.pi*float(TOKENS_AC[-1]) #OMEGA, Angular frequency = 2*pi*frequency and frequency is the last token 

            if INDEX_AC < INDEX_END: #If '.ac' is within the circuit, then put INDEX_AC = 0  
                INDEX_AC =0      

            if INDEX_CIRCUIT > INDEX_END: #If '.circuit' is present after '.end', tyhen it is invalid circuit definition
                print("ERROR. INVALID CIRCUIT DEFINITION IN THE GIVEN INPUT FILE.")

            else: #Solve the circuit
                
                for i in range(INDEX_CIRCUIT + 1, INDEX_END): #We want tokens after '.circuit' and before '.end' so the range is used like that
                    ElementName, N1, N2, *Val = (lines[i].split('#')[0]).split() #Taking the strings(of each line) till '#'(ignoring '#') and then again splitting to get individual tokens 
                    #Val is a pointer that points to all the stuff thaty is there after node 2!
                    
                    TOKENS.append((lines[i].split('#')[0].split())) #Updating the TOKENS
                    
                    #Using our defined classes to store the values precisely
                    if ElementName[0] == 'V': #Voltage source
                        component = VoltageSource(ElementName, N1, N2, Val)

                    if ElementName[0] == 'I': #Current source
                        component = CurrentSource(ElementName, N1, N2, Val)

                    if ElementName[0] == 'R': #Resistor
                        component = Resistor(ElementName, N1, N2, Val)                     
                    
                    if ElementName[0] == 'L': #Inductor
                        component = Inductor(ElementName, N1, N2, Val)
                    
                    if ElementName[0] == 'C': #Capacitor
                        component = Capacitor(ElementName, N1, N2, Val)

                    ### ~~~~~~~~~ TOTAL NUMBER OF NODES ~~~~~~~~~~~ ###
                    
                    #The below portion of code will find the total number of nodes present in the circuit by using the node tokens given in the netlist file 
                    if component.N1 == 'GND': #Reading the first node(N1) of each component
                        i = 0 #All nodes are numbered wrt "GND", so for GND we take i = 0 
                    else:
                        i = int(component.N1) #Other than GND others are node, "i" stores the first node number(or positon) wrt "GND" 

                    if component.N2 == 'GND': #Reading the second node(N2) of each component
                        j = 0 #Similar to above
                    else: 
                        j = int(component.N2) #Other than GND others are node, "i" stores the second node number(or positon) wrt "GND"

                    if NODE_COUNT <= i: 
                        NODE_COUNT = i  #Assigning the max. value of 'i' to NODE_COUNT based N1(NODE1) values   

                    if NODE_COUNT <= j:
                        NODE_COUNT = j  #Assigning the max. value of 'j' to NODE_COUNT based N2(NODE2) values
                        #Finally the max. value of all N1's and N2's give the NODE_COUNT
                    
                    ### ~~~~~~~~~ TOTAL NUMBER OF VOLTAGE SOURCES ~~~~~~~~~~~ ###
                    
                    if component.ElementName[0] == "V": #If the first token contains "V", then it is a voltage source
                        VoltageSource_Count += 1 #Updating the voltage count

                    ### ~~~~~~~~~ INITIALISING THE MATRIX ~~~~~~~~~~~ ###
                    #Using list comprehension
                    
                    #Conductance matrix 
                    M = [ [ 0 for i in range(NODE_COUNT + VoltageSource_Count) ] for j in range(NODE_COUNT + VoltageSource_Count) ] #Initialising the matrix of required size                         
                    
                    #Independent source Vector(column matrix) 
                    b = [ [ 0 for i in range(NODE_COUNT + VoltageSource_Count) ] for j in range(1) ] #Initialising the matrix of required size

                    
                    ### ~~~~~~~~~ FORMING THE MATRIX, M BY USING THE MNA STAMPS OF ALL THE ELEMENTS ~~~~~~~~~~~ ###

                    for l in range(0, len(TOKENS)): #Running a for loop through all the elements
                        
                        #Adding the MNA Stamp of the Resistor to the conductance matrix  
                        if component.ElementName[0] == "R": #If the element is a resistor
                            XR = float((component.Val[0])) #Reactance due to resistor

                            if component.N1 == "GND": #The MNA stamp of Resistor looks like: |` G -G `| If either N1 = 0 or N2 = 0, then Stamp is [G] 
                                k = 0              #   (2x2) Square matrix                |_-G  G _| (1st column is kth node, 2nd column is nth node)
                            else:
                                k = int(component.N1) #Take Node1 Value

                            if component.N2 == "GND": #If N2 is GND, l = 0  
                                l = 0
                            else:
                                l = int(component.N2) #Else take the Node2 Value   

                            if k == 0:
                                M[l-1][l-1] += 1/XR 

                            if l == 0:
                                M[k-1][k-1] += 1/XR #Also, we are updating the stamp(by using '+=', going on to add) and adding to main matrix M    

                            if k*l != 0: #If none of the two nodes are "GND", then
                                M[k-1][k-1] += 1/XR #11th entry
                                M[k-1][l-1] += -1/XR #12th entry
                                M[l-1][k-1] += -1/XR #21th entry
                                M[l-1][l-1] += 1/XR #22th entry 

                        #Adding the MNA Stamp of the Inductor to the conductance matrix  
                        if component.ElementName[0] == "L" and INDEX_AC != 0: #If the element is a inductor, and also INDEX_AC != 0 indicates that there is an ac source  
                            XL = complex(0, (OMEGA*float(component.Val))) #Reactance due to Inductor, XL = 0 + jwL

                            if component.N1 == "GND": #The MNA stamp of Inductor(only in AC CASE) looks like: |` XL -XL `| If either N1 = 0 or N2 = 0, then Stamp is [XL] 
                                k = 0              #   (2x2) Square matrix                                 |_-XL  XL _| (1st column is kth node, 2nd column is nth node)
                            else:
                                k = int(component.N1) #Take Node1 Value

                            if component.N2 == "GND": #If N2 is GND, l = 0  
                                l = 0
                            else:
                                l = int(component.N2) #Else take the Node2 Value   

                            if k == 0:
                                M[l-1][l-1] += 1/XL 

                            if l == 0:
                                M[k-1][k-1] += 1/XL #Also, we are updating the stamp(by using '+=', going on to add) and adding to main matrix M    

                            if k*l != 0: #If none of the two nodes are "GND", then
                                M[k-1][k-1] += 1/XL #11th entry
                                M[k-1][l-1] += -1/XL #12th entry
                                M[l-1][k-1] += -1/XL #21th entry
                                M[l-1][l-1] += 1/XL #22th entry

                        #Adding the MNA Stamp of the Capacitor to the conductance matrix  
                        if component.ElementName[0] == "C" and INDEX_AC != 0: #If the element is a capacitor, and also INDEX_AC != 0 indicates that there is an ac source  
                            XC = complex(0, -1/(OMEGA*float(component.Val))) #Reactance due to Capacitor, XC = 0 - (j/wc)

                            if component.N1 == "GND": #The MNA stamp of Capacitor(only in AC CASE) looks like: |` XC -XC `| If either N1 = 0 or N2 = 0, then Stamp is [XL] 
                                k = 0              #   (2x2) Square matrix                                 |_-XC  XC _| (1st column is kth node, 2nd column is nth node)
                            else:
                                k = int(component.N1) #Take Node1 Value

                            if component.N2 == "GND": #If N2 is GND, l = 0  
                                l = 0
                            else:
                                l = int(component.N2) #Else take the Node2 Value   

                            if k == 0:
                                M[l-1][l-1] += 1/XC

                            if l == 0:
                                M[k-1][k-1] += 1/XC #Also, we are updating the stamp(by using '+=', going on to add) and adding to main matrix M    

                            if k*l != 0: #If none of the two nodes are "GND", then
                                M[k-1][k-1] += 1/XC #11th entry
                                M[k-1][l-1] += -1/XC #12th entry
                                M[l-1][k-1] += -1/XC #21th entry
                                M[l-1][l-1] += 1/XC #22th entry  

                        ##Adding the MNA Stamp of the  to the conductance matrix                                                                            
                        if component.ElementName[0] == "V":
                            VoltageTemp += 1
                            if component.N1 == "GND":
                                k = 0
                            else:
                                k = int(component.N1)
                                M[int(NODE_COUNT + VoltageTemp - 1)][k-1] = 1
                                M[k-1][int(NODE_COUNT + VoltageTemp - 1)] = -1
                                if component.Val[0] == "ac":								
                                    b[0][int(NODE_COUNT + VoltageTemp - 1)] = complex(float(component.Val[-2])*math.cos(float(component.Val[-1])),float(component.Val[-2])*math.sin(float(component.Val[-1])))
                                else:								
                                    b[0][int(NODE_COUNT + VoltageTemp - 1)] = float(component.Val[-1])
                            if component.N2 == "GND":
                                l = 0
                            else:
                                l = int(component.N2)
                                M[int(NODE_COUNT + VoltageTemp - 1)][l-1] = -1
                                M[l-1][int(NODE_COUNT + VoltageTemp - 1)] = -1
                                if component.Val[0] == "ac":
                                    b[0][NODE_COUNT +  VoltageTemp - 1] = -complex(float(component.Val[-2])*(0.5)*math.cos(float(component.Val[-1])),float(component.Val[-2])*(0.5)*math.sin(float(component.Val[-1])))
                                else:
                                    b[0][NODE_COUNT + VoltageTemp - 1] = -float(component.Val[-1])

                    
                    ### ~~~~~~ SOLVING THE MATRIX USING NUMPY LIBRARY ~~~~~~~ ###
                    Solution = np.linalg.solve(M, b[0]) #Equation is: Mx = b, here x = Solution
                    
                    #Displaying the outpot matrix
                    print("The Required solution matrix is: ")
                    print(Solution)
                    
                    
                    print("The Required values are: ")
                    for i in range(0,NODE_COUNT + VoltageSource_Count):
                        if INDEX_AC == 0:       # for dc sources the node wise voltage and current through the voltage source.
                            if i < NODE_COUNT: 
                                print("Voltage at node %d is"%(i + 1))
                            if i >= NODE_COUNT: 
                                print("Current through V%d source is"%(i - NODE_COUNT + 1))
                            print(float("%.2f" %Solution[i]))
                        else:		# for dc sources the node wise voltage and current through the voltage source.
                            if i < NODE_COUNT: 
                                print("Voltage at node %d is"%(i + 1))
                            if i >= NODE_COUNT: 
                                print("Current through V%d source is"%(i - NODE_COUNT + 1))
                            print(Solution[i]) 

elif len(argv) ==1 : 
    print("Present usage: %s \nError: No Input Found \nExpected: FileName.py FileName.netlist" % argv[0])
else: 
    print("Present usage: %s \nError: Too Many Inputs Found \nExpected: FileName.py FileName.netlist" % argv)                                             		


                     






                







