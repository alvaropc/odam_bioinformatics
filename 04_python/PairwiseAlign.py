#Python 3.4.3 under Ubuntu 14.04
import PairwiseAlign
import matplotlib.pyplot as plt


def MGenerator(seq1,seq2):
	'''Generates a matrix from 2 sequences 
	All the code works with this sequences Matrix
	If one sequence is shorter than the other, it is saved in the first place,
	in order to create less rows as possible.'''

	lens = [seq1,seq2]
	if len(lens[0])>len(lens[1]):
		maximun=lens[0] 
		minimun=lens[1]
	elif len(lens[0])<len(lens[1]):
		maximun=lens[1] 
		minimun=lens[0]
	else:
		maximun=lens[1]
		minimun=lens[0]

	m1=[]
	for i in minimun:
		m1.append(i)
	m2=[]
	for i in maximun:
		m2.append(i)
	M = [m1,m2]
	return(M)


def M_count(M):	
    '''This function counts the matches, missmatches and gaps points
    of a sequence matrix. It could be the initial sequence Matrix
    or the final alingment Matrix'''
    res = []
    count=0	
    for i in M[0]:
        if i == M[1][count]:
            res.append('*')
        elif M[1][count] == '-' or M[0][count]=='-':
            res.append('-')
        else:
            res.append(' ')
        count +=1
    return(res)

def missmatchesCounter(M):
	'''This function returns the positions of the missmatches inside a sequence Matrix'''
	res = PairwiseAlign.M_count(M)
	miss = []
	c= 0
	for i in res:
		if i == ' ':
			miss.append(c)
		c+=1
	if miss==[]:
		print('no missmatches found')
	return(miss)

def matchesCounter(M): 
	'''Returns the positions of the matches inside a sequence Matrix'''
	res = PairwiseAlign.M_count(M)
	matches = []
	c= 0
	for i in res:
		if i == '*':
			matches.append(c)
		c+=1
	if matches==[]:
		print('No matches found')
	return(matches)

def gapCounter(M):
	'''Returns the positions of gaps inside a sequence Matrix'''
	res = PairwiseAlign.M_count(M)
	gap = []
	c= 0
	for i in res:
		if i == '-':
			gap.append(c)
		c+=1
	if gap==[]:
		print('No gaps found')
	return(gap)




def substitutionM(M):
	'''Creation of the substitution matrix'''
#First we create an empty matrix
	empty=[]
	for row in range(0,len(M[0])):
		empty.append([])
		for column in range(0,len(M[1])):
			empty[row].append('')
#Subtitution matrix generation
	for row in range(0,len(empty)):
		for column in range(0,len(M[1])):
			if M[1][column]==M[0][row]:
				empty[row][column]=5
			elif M[1][column]!=M[0][row]:
				if M[1][column]=='A'or'G' and M[0][row]=='C' or 'T': #Transversion
					empty[row][column]=-4
				if M[1][column]=='C' or 'T' and M[0][row]== 'A' or 'G':
					empty[row][column]=-4
				if M[1][column]=='A' and M[0][row]=='G': #Transition
					empty[row][column]=-2
				if M[1][column]=='G' and M[0][row]== 'A': #transition
					empty[row][column]=-2
				if M[1][column]=='C' and M[0][row]=='T': #Transition
					empty[row][column]=-2
				if M[1][column]=='T' and M[0][row]== 'C': #transition
					empty[row][column]=-2
	subs=empty
	return(subs)

				

def dinamicM(M):
	'''Creation of the dinamic matrix'''
	subs=PairwiseAlign.substitutionM(M)
#We first create an empty dinamic matrix
	dinamicEmpty=[]
	for row in range(0,len(M[0])+1):
		dinamicEmpty.append([])
		for column in range(0,len(M[1])+1):
			dinamicEmpty[row].append('')
#Some parameters are added to the empty dinamic matrix like first row and column with gap penaltys primera fila y columna, gap penalty
	dinamicEmpty[0][0]=0
	gapPenalty=-5
	gapPenalty1=-5
	for i in range(0,len(dinamicEmpty)):
		if i==0:
			for d in range(1, len(dinamicEmpty[i])):
				dinamicEmpty[i][d]=gapPenalty
				gapPenalty-= 5
		else:
			dinamicEmpty[i][0]=gapPenalty1
			gapPenalty1-=5
#Adding score values to the dinamic matrix
	gapPenalty=-5
	for row in range(1,len(dinamicEmpty)):
		for column in range(1,len(dinamicEmpty[0])):
			dinamicEmpty[row][column]=max(dinamicEmpty[row-1][column-1]+subs[row-1][column-1],
										  dinamicEmpty[row-1][column]+gapPenalty,
										  dinamicEmpty[row][column-1]+gapPenalty)
	dinamic=dinamicEmpty
	return(dinamic)

#This function let the user plot the dinamic matrix
def dinamicPlot(dinamic):
	
	plt.imshow(dinamic, interpolation='nearest', cmap=plt.cm.ocean)
	plt.colorbar()
	return(plt.show())
	
def dinamicAnalizer(M):
	'''Dinamic Matrix Analysis.
	As a result you will have the score of the
	alignment, a print of the alignment, and the 
	alignment in an list with the points printed
	in the alignment as second item of the list '''

	dinamic=PairwiseAlign.dinamicM(M)
#First we find the position in the matrix where the score is the highest, to start the analysis
	maximun=[]
	for row in range(1,len(dinamic)):
		for column in range(1,len(dinamic[0])):
			maximun.append(dinamic[row][column])
	for row in range(1,len(dinamic)):
		for column in range(1,len(dinamic[0])):
			if max(maximun)==dinamic[row][column]:
				initialPosition=[[row,column]]
#While we are not in the initial position of the matrix, the code keep runing				
	positions=initialPosition
	counter=0
	sumMaximo=[]
	while positions[len(positions)-1]!=[0,0]:
		maximo = max(dinamic[positions[counter][0]-1][positions[counter][1]-1],
					 dinamic[positions[counter][0]-1][positions[counter][1]],
					 dinamic[positions[counter][0]][positions[counter][1]-1])
		if maximo==dinamic[positions[counter][0]-1][positions[counter][1]-1]:
			positions.append([positions[counter][0]-1,positions[counter][1]-1])
		elif maximo==dinamic[positions[counter][0]-1][positions[counter][1]]:
			positions.append([positions[counter][0]-1,positions[counter][1]])
		else:
			positions.append([positions[counter][0],positions[counter][1]-1])
		sumMaximo.append(maximo)
		counter+=1

	score=sum(sumMaximo)
	print('score:',score)
#Positions are calculated on dinamic matrix, in order to match them with the sequence,
#the positions are corrected and matching with the substituion matrix too
	positionsC=[] #Corrected positions
	for i in range(0,len(positions)-1):
		positionsC.append([positions[i][0]-1,positions[i][1]-1])
	
	positionsR=[]
	for i in reversed(range(0,len(positionsC))):
	    positionsR.append(positionsC[i])

	al0=[] # M[0] equivalent
	al1=[] # M[1] equivalent
	for p in (range(0,len(positionsR))):
		al0.append(M[0][positionsR[p][0]])
		al1.append(M[1][positionsR[p][1]])
	#Resulting alignment
	Al=[al0,al1]
	miss=PairwiseAlign.missmatchesCounter(Al)
	#Addition of gaps
	for i in range(0,len(positionsR)):
		if positionsR[i][0]==positionsR[i-1][0]:
			if i in miss:
				Al[0][i]='-'
			if i-1 in  miss:
				Al[0][i]='-'
		if positionsR[i][1]==positionsR[i-1][1]:
			if i in miss:
				Al[1][i]='-'
			if i-1 in  miss:
				Al[0][i]='-'
	points=PairwiseAlign.M_count(Al)
	print(Al[0],Al[1],points,sep='\n') #This line can give you problems in python 2
#Here I should revise again the alingment and make better the Gap asignment
	return(Al,points)
	
		
		
		
		
		
		
		
		
		
		
		
	
