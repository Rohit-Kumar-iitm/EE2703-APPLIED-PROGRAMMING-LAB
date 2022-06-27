"""
        EE2703 Applied Programming Lab - 2022
            Assignment 1
            NAME: ROHIT KUMAR
            ROLLL NO.: EE20B111
            DATE: 26-01-2022
"""

from sys import argv, exit

#Assigning '.circuit' and '.end' to constant variables "CIRCUIT" and "END"
CIRCUIT_START = '.circuit'
CIRCUIT_END = '.end'

global c #global declaration of variable, c is used to print the line number of the element
c = 1

def WordExtract(line): #Function to extract tokens of each element of a branch, and to print it's connections and value
    global c
    word = line.split('#')[0] #Taking only the strings before the '#', i.e., ignoring the comment
    token = word.split() #Extracting each token(storing it in word)
    
    if len(token) == 4: #For independent sources, number of tokens = 4
        TypeOfElement = token[0]
        FromNode = token[1]
        ToNode = token[2]
        ValueOfElement = token[3]
        print("The Element in Line : ", c)
        print("Type: %s , From Node: %s , To Node: %s , Value of Element: %s" % (TypeOfElement, FromNode, ToNode, ValueOfElement) )
        c += 1 #Incrementing c
        return [TypeOfElement, FromNode, ToNode, ValueOfElement]

    elif len(token) == 5: #For Current controlled sources, number of tokens = 5
        TypeOfElement = token[0]
        FromNode = token[1]
        ToNode = token[2]
        VoltageSource = token[3]
        ValueOfElement = token[4]
        print("The Element in Line : ", c)
        print("Type: %s , From Node: %s , To Node: %s , Voltage Source: %s, Value of Element: %s" % (TypeOfElement, FromNode, ToNode, VoltageSource, ValueOfElement) )
        c += 1 #Incrementing c
        return [TypeOfElement, FromNode, ToNode, VoltageSource, ValueOfElement]
    
    elif len(token) == 6: #For Voltage controlled sources, number of tokens = 6
        TypeOfElement = token[0]
        FromNode = token[1]
        ToNode = token[2]
        VoltageSource1 = token[3]
        VoltageSource2 = token[4]
        ValueOfElement = word[5]
        print("The Element in Line : ", c)
        print("Type: %s , From Node: %s , To Node: %s , Voltage Source1: %s, Voltage Source2: %s, Value of Element: %s" % (TypeOfElement, FromNode, ToNode, VoltageSource1, VoltageSource2, ValueOfElement) )
        c += 1 #Incrementing c
        return [TypeOfElement, FromNode, ToNode, VoltageSource1, VoltageSource2, ValueOfElement]
    else: #If no. of tokens is not equal to either 4 or 5 or 6, it throws an error
        print("Error. The Element in Line : %d is not defined properly" % c)
        c += 1
        return[]
            

#Checking for correct number of arguments are given or not
if len(argv) != 2:
    print('\nUsage: %s <inputfile>' % argv[0])
    print("Make sure that number of arguments are 2")
    exit()

#Checking for the error in the input file, if there's an error an error message will be displayed using the "try"
try:
    with open(argv[1]) as f: 
        lines = f.readlines() 
        f.close() 
        start = -1; end = -2
        NumStart = 0; NumEnd = 0 #NumStart counts how many times '.circuit' was repeated, similarly NumEnd for '.end'
        for line in lines:              #Extracting the circuit definition from the start and end lines
            if CIRCUIT_START == line[:len(CIRCUIT_START)]:
                start = lines.index(line)
                NumStart += 1 
            elif CIRCUIT_END == line[:len(CIRCUIT_END)]:
                end = lines.index(line)
                NumEnd += 1
                break
        
        if start >= end:  #Checking for valid circuit definition, should start with '.circuit' and end with '.end' 
            print("Error. Invalid circuit definition. It should start with '.circuit' and end with '.end'")
            exit(0)
        
        if NumStart != 1 or NumEnd != 1: #There should only one beginning('.circuit') and one end('.end')
            print("Error. The netlist must contain only one '.circuit' and only one '.end'. Please check for it. Your netlist may contain multiple '.circuit' or '.end'")
        
        for line in lines[start+1:end]:
            WordExtract(line) #Calling the function to print the element details
        
       
        print("\nThe required(reversed) output is: \n")
        for line in reversed([' '.join(reversed(line.split('#')[0].split())) 
            for line in lines[start+1:end]]): #This takes care of reversing the lines and then individually reversing the elements of each line
                print(line)     

except IOError:
    print('Error. Invalid file, Make sure that the entered file name is correct.')
    exit()


