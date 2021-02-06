#!/usr/bin/python3
import argparse
import sys
import re 
import pandas as pd
import statistics
import numpy as np
import os 

# VERSION AB MI, 3.SEP 2020

#-----------------------------Argument Parser. Input Interface--------------------------------
#					   WORKING
#		add a new argument for turning verbosity on!

parser = argparse.ArgumentParser(description="DESCRIPTION: This program takes .fastq input from the NCBI SRA and gives you general information about all spots (NC contents, means, stdev)", epilog="First specify your input file! Then add additional arguments. A new folder for your project will be created in your wd, based on the inputfile name. To download and run multiple runs of a Project at once, please use sralib_pipeline.sh")

parser.add_argument("-i", "--input", required=True, dest="inputfile", type=argparse.FileType("r"),help="(REQUIRED) specify your input file (fastq)")

parser.add_argument("-o", "--output", dest="outfilename", nargs="?", default=sys.stdout, help="specify outfilename with -o filename . Extension is automatical. If only -o is used, outfile is created with INPUTFILE.txt as default name automatically. If argument is not used, output is printed in stdout (commandline)")

parser.add_argument("-g", "--graphics", dest="graphics", action="store_const", const=True, default=False, help="add if you want graphics to be created. Default is False to reduce runtime.")

parser.add_argument("-s", "--spotlist", dest="spotlist", help="If you want statistics for a subset of spots only, enter a .txt file with line separated spot names.")
#iwas mit spots=[], if spotlist is not False: while open (spotlist, "r") as spotlist: for line in spotlist: spots.append(line).
 


args = parser.parse_args()



# preparing naming convention for folders and graphical outputs.

i= sys.argv[2][:sys.argv[2].rfind(".")]
foldername= "{}_Project".format(i)
folderpath= foldername+"/"

try:
	os.mkdir(foldername)
	print ("\n", "Your folder", foldername, "has been created")
except FileExistsError:
	print ("Directory with INPUTNAME already exists. Please rename existing dir if you want to run the same project again.")
	print (FileExistsError)
	quit()


#-----------------------------------Dictionary creation---------------------------------------
#		create empty dic first WORKING
#		spotlist implementation WORKING	
#		Total A, T, C, G counter for all spots and dict entry counter WORKING
#		format each individual line WORKING BUT DO WE NEED QUALITY HEADER??
#		update dict only when matching the spotlist WORKING	
#		length counter per len(string) instead of int(headerelement) WORKING
#		dict.update with or without myspotlist input WORKING

dict1={}


 
myspots=[]
if args.spotlist is not None: 
	with open (args.spotlist, "r") as spotlist: 
		for line in spotlist: 
			line=line.strip().lstrip("@")
			liste=line.split(" ")
			line=liste[0]
			myspots.append(line) #add (line.rstrip())?? it did print the \n now

#create a list with the input names from the spotlist and trim header same as done for the input file so these two can be matched later on.


c=1
acon=0
tcon=0
ccon=0
gcon=0
# c line counter for the 4 lines. 
# acon, tcon, ccon, gcon are the "a-content" counters for overall A occurence over all entries.
	
for line in args.inputfile:
	if c==1:
		line1=line.strip().lstrip("@")
		liste1=line1.split(" ")
		name=liste1[0]
		c=c+1
	elif c==2:
		line=line.strip()
		seq=line
		counta=line.count("A")
		countt=line.count("T")
		countc=line.count("C")
		countg=line.count("G")
		acon=acon+counta # the just counted A of this line gets added to overall a counter.
		tcon=tcon+countt
		ccon=ccon+countc
		gcon=gcon+countg	
		leng=len(line)
		percentageA= counta/leng*100
		percentageT= countt/leng*100
		percentageC= countc/leng*100
		percentageG= countg/leng*100
		c=c+1
	elif c==3:
		c=c+1
	else:
		line=line.strip()
		qscore=line
# add the condition that dict1 is only expanded if the full line 1 is in myspots!		
		if args.spotlist is not None:
			for el in myspots:
				if el == name:
					dict1.update ({name:{
					"length":leng, 
					"sequence":seq, 
					"qualityScore":qscore, 
					"Acount":counta,
					"A%":percentageA, 
					"Tcount":countt, 
					"T%":percentageT,
					"Ccount":countc, 
					"C%":percentageC,
					"Gcount":countg,
					"G%":percentageG
					}})
		else: 		# if there is NO myspotlist input
			dict1.update ({name:{
			"length":leng, 
			"sequence":seq, 
			"qualityScore":qscore, 
			"Acount":counta,
			"A%":percentageA, 
			"Tcount":countt, 
			"T%":percentageT,
			"Ccount":countc, 
			"C%":percentageC,
			"Gcount":countg,
			"G%":percentageG
			}})
		c=1
# it goes through the lines of the input file (as r),and strips the whitespaces, splits the string from the headers to separate the name from additions like @ and length=x. for line 2 it counts the length, counts all the NCs separately, as well as adds them to the overall count. in else, the dict1 is updated by the remodeled lines, either for all the names in myspotlist, or for all spots from infile. one spot (name key) has the keys length, sequence and qscore as given in the fastq infile, as well as newly added counts of individual NCs. The overall count of NCs is stored in variables and is used in section CLACULATIONS for the generalstats of the input. 
originallen=len(dict1)

print ("----------------------COMMENTS----------------------")
print("Your original library has", originallen, "entries.")

# remove entries that contain NEITHER NC type!!
deletelist=[]
for element in dict1: 
	if (dict1[element]["Acount"] == 0) and (dict1[element]["Tcount"] == 0) and (dict1[element]["Ccount"] == 0) and (dict1[element]["Gcount"] == 0):
		deletelist.append(element)


# Activate for a df of the deleted spots.
#newdict={}
#for name in deletelist: 
#	newdict[name]={}
#	newdict[name]["Acount"]= dict1[name]["Acount"]
#	newdict[name]["Tcount"]= dict1[name]["Tcount"]
#	newdict[name]["Ccount"]= dict1[name]["Ccount"]
#	newdict[name]["Gcount"]= dict1[name]["Gcount"]

#newdf=pd.DataFrame.from_dict(newdict, orient="index")
#print (newdf)



for name in deletelist: 
	del dict1[name]
print (len(deletelist), "spots were deleted, because they were empty (contained neither A, T, C, or G).")	
print ("This makes", len(deletelist)/originallen*100, "% of your library.")
print ("Your new library has", len(dict1), "entries.")


# test if there are empty values in the dict.
#for element in dict1: 
#	if dict1[element]["Acount"] == 0:
#		print (element, "A None")
#	if dict1[element]["Tcount"] == 0:
#		print (element, "T None")
#	if dict1[element]["Ccount"] == 0:
#		print (element, "C None")
#	if dict1[element]["Gcount"] == 0:
#		print (element, "G None")
#	if "N" in dict1[element]["sequence"]:
#		print(element, "contains N")



df=pd.DataFrame.from_dict(dict1, orient="index")
#print (df["Tcount"].dtypes)
# the counts are int64, the sequence is object (so string). should be fine. 


# Activate to print my dict1 as a table with pandas DataFrame function. 
# The axis are reversed for our use, so we have to transpose them. 
# It actually has no index... that one time I wanted one... 	
#print (df)
#df_transposed = df.transpose()   #or df.T   not neccessary bc of parameter orient=
#print (df_transposed)



#---------------------------------DI-NCs-PER-SPOT-CALCULATIONS------------------------------------
#					   WORKING

#					DI-NC PER SPOT


ncliste= ["A", "T", "C", "G"] 		# create a list of poss. Ncs
liste2= [] 				# create empty list to be filled w/ all poss. pairs
for e in ncliste:
	for j in ncliste: 
		liste2.append(e+j)
# to the new list add every element of the first list per element of the first list. 
# so it makes a matrix and pairs the 4 characters with each character.

pairstatsdict= {}
for element in dict1.keys():
	leng= int(dict1[element]["length"]) #just for calc reasons
	possiblepairs= int(leng-1) #the total possible amount of pairs
	pairstatsdict[element]= {}
	for e in liste2:
		count1=int(dict1[element]["sequence"].count(e))
		percentage=count1/possiblepairs*100
		pairstatsdict[element][e]= count1
		pairstatsdict[element][e+"%"]= percentage

pairstatsdf= pd.DataFrame.from_dict(pairstatsdict, orient="index")




#-----------------------------BIG-DICT1 (df)-DATAFRAME FORMATTING---------------------------------
#					   WORKING

spotresult= pd.concat([df, pairstatsdf], axis=1, sort=False)

#print (list(result.columns)) 
# it will give you a list of the column names, so I could rearrange them. left in for maintainance purpose. not very elegant, maybe... 
spotresult= spotresult[['length', 'Acount', 'A%', 'Tcount', 'T%', 'Ccount', 'C%', 'Gcount', 'G%', 'AA', 'AA%', 'AT', 'AT%', 'AC', 'AC%', 'AG', 'AG%', 'TA', 'TA%', 'TT', 'TT%', 'TC', 'TC%', 'TG', 'TG%', 'CA', 'CA%', 'CT', 'CT%', 'CC', 'CC%', 'CG', 'CG%', 'GA', 'GA%', 'GT', 'GT%', 'GC', 'GC%', 'GG', 'GG%', 'sequence', 'qualityScore']]

#print (spotresult)

# basically concated verion of dict1(df) and the pairstatsdf! then columns rearranged


#------------------------------------TOTAL STATISTICS---------------------------------------------

# somehow make the calcs file-printable. maybe rather a df with collumn names, no index and the specific entries. # make separate dfs for general file stats, NC total stats, spot specific stats.
#moved the prints to the verbosity section

# TOTAL NUMBER OF NCs OF ALL SPOTS. Correct length!!! 	WORKING

totalncs = 0
for element in dict1.keys():
	l= int(dict1[element]["length"])
	totalncs= int(totalncs+l)			
# it creates counters of 0, then runs over every dict entry and states that the l is the entry's length value. then it adds the entry's l to the totalncs counter and goes to the next entry.
print ("Total NCs count of all spots: ", totalncs)
print ("The mean length per spot is ", totalncs/len(dict1))

#					Single NC Totals


# SUM OF ALL As, Ts, Cs, Gs with sum()!!! 	WORKING

# Tested against the for element in dict1: allAs=allAs+dict1[element]["Acount"] version and same result!!! 
# hence sum() was chosen, as it is shorter!
#column= list(result["Acount"])
#print (column)		# works in printing the contents of the column as list, but not needed here
Asum= sum(spotresult["Acount"])
Tsum= sum(spotresult["Tcount"])
Csum= sum(spotresult["Ccount"])
Gsum= sum(spotresult["Gcount"])


# PERCENTAGE OF ALL As ON TOTAL NCs LENGTH	WORKING

percenta= Asum/totalncs*100
percentt= Tsum/totalncs*100
percentc= Csum/totalncs*100
percentg= Gsum/totalncs*100


# MEAN OF ALL As... 				WORKING

Amean= statistics.mean(spotresult["Acount"])
Tmean= statistics.mean(spotresult["Tcount"])
Cmean= statistics.mean(spotresult["Ccount"])
Gmean= statistics.mean(spotresult["Gcount"])
# calculates the mean a-content per spot. so total Asum/nr.of.spots. saved in variables.


# SAMPLE STANDARD DEVIATION FROM MEAN (the square root of the sample variance) 	   PROBS WORKING

stdevA= statistics.stdev(spotresult["Acount"])
stdevT= statistics.stdev(spotresult["Tcount"])
stdevC= statistics.stdev(spotresult["Ccount"])
stdevG= statistics.stdev(spotresult["Gcount"])


# MAKE A STANDALONE DF with the SINGLE NC stats. 	WORKING

index1= ["A", "T", "C", "G"] 
# create index list to replace the automatic numbers index.
generalstatsdict= {	"sum":[Asum, Tsum, Csum, Gsum], 
			"%":[percenta, percentt, percentc, percentg], 
			"mean-per-spot":[Amean, Tmean, Cmean, Gmean],
			"stdev":[stdevA, stdevT, stdevC, stdevG]
			}
			
generalstatsdf= pd.DataFrame (generalstatsdict, index=index1)
#print (generalstatsdf)


#					Di-NC Totals


# CALC TOTAL SUM, PERC, MEAN, STDEV FOR Di-NCs			WORKING

newpairstatsdf= pairstatsdf
for col in newpairstatsdf: 
	if "%" in col:
		del newpairstatsdf[col]
#print ("newpairstatsdf= ", newpairstatsdf)

# In order to select all of the rows and some columns, we use single colon [:] 
# to select all of rows and list of some columns which we want to select like this:
#Dataframe.loc[[:, ["column1", "column2", "column3"]]
#or:	row2 = data.iloc [:, [1, 2]
# couldnt find one for "every second col", so I went with excluding "%"-containing cols


dincstatsdict={}
for col in newpairstatsdf: 
	summe= sum(newpairstatsdf[col])
	posspairs= int(totalncs-1)
	percent= summe/posspairs*100
	mean= statistics.mean(newpairstatsdf[col])
	stdev= statistics.stdev(newpairstatsdf[col])
	dincstatsdict[col] = {}
	dincstatsdict[col]["sum"] = summe
	dincstatsdict[col]["%"] = percent
	dincstatsdict[col]["mean-per-spot"] = mean
	dincstatsdict[col]["stdev"] = stdev
	

# MAKE A STANDALONE DF with the Di-NC PAIR stats		WORKING

dincstatsdf= pd.DataFrame.from_dict(dincstatsdict, orient="index")
#print ("dincstatsdf", dincstatsdf)

# CONCAT THE TOTAL SINGLE AND TOTAL Di-NC DFs!			WORKING

generalresult= pd.concat([generalstatsdf, dincstatsdf], sort=False)
#print("generalresult", generalresult)




#---------------------------------OPENING AND WRITING OUTPUT FILES--------------------------------
#						WORKING

# ich wuerds erst spaeter oeffnen, damit es nicht so lange offen sein muss, waehrend man den input einliest und berechent. es koennte schon offen sein und alles was ich printe, geben ich is output. halte ich aber fuer nicht so gut, besser ich speicher das alles in variablen und printe dann den fertigen df hinein! Heisst aber auch ich muss vlt beim add_argument den filetype W raus tun und auch das erst spaeter oeffnen. einfach den eingegbenen namen als string speichern und nachher sagen open stringfile.txt as x.

if args.outfilename==None:
	ii=sys.argv[2][:sys.argv[2].rfind(".")]
	outfile1="{}_TotalStatistics.txt".format(ii)
	generalresult.to_csv(folderpath+outfile1, sep="\t")
	outfile2="{}_SpotStatistics.txt".format(ii)
	spotresult.to_csv(folderpath+outfile2, sep="\t")
# [:x.rfind(".")] finds the first dot from right in i and removes everything from that on. so it cuts the .fastq. # ("{}whatever".format(var)) takes the var after .format and puts it into place where the {} are.

else:
	outfile1="{}_TotalStatistics.txt".format(args.outfilename)
	generalresult.to_csv(folderpath+outfile1, sep="\t")
	outfile2="{}_SpotStatistics.txt".format(args.outfilename)
	spotresult.to_csv(folderpath+outfile2, sep="\t")


# you can also open a file in append mode, so with a instead of w. then it will always be appended and not overwritten. might be problematic when re-running with the same name after an error. alternatively maybe not use .write in the argument but .append?

#________________________________________________________________________________________________
#______________________________STOPS HERE WHEN NO -g FLAG IS GIVEN_______________________________

if args.graphics is False:
	print("END. Find your Outputs in the folder stated above.")
	quit()

print ("Your graphics are being created...")
#-----------------------------------------GRAPHICS----------------------------------------------


#   				     NEW NUMPY ARRAY FOR GRAPHICS!
#						(WORKING)

# THERE WERE EMPTY SPOTS, WHICH IS PROBABLY WHY THE HISTOGRAM DOESNT WORK!!!
# find a way to exclude them. Excluded them and set Count column of gdf to int!!

# creating an empty list and making a second nested list with NC type and count value.
# then extending list1 by list2. 

glist1= []
for element in dict1:
	glist2= [
		["A", dict1[element]["Acount"]], 
		["T", dict1[element]["Tcount"]], 
		["C", dict1[element]["Ccount"]], 
		["G", dict1[element]["Gcount"]]
		]
	glist1.extend(glist2)
#print(glist1)

# convert list to an array, so that we have "columns"
garray1= np.asarray(glist1)
#print (garray1)

# convert array to a pandas df. 
gdf = pd.DataFrame(garray1, columns=["NC", "Count"])
#print (gdf)
#print (gdf["Count"].dtypes)
# Count column was OBJECT!! So converting it to an integer now:

gdf["Count"] = gdf["Count"].astype(int)
#print (gdf["Count"].dtypes)
# NOW IT IS INT64. Perfect :)


# appending an array is working but not the best way. 
# it is supposedly slower than appending a list and converting that to array.
#np1= np.empty((0, 2))
#for element in dict1:
#	np2= np.array(
#		[["A", dict1[element]["Acount"]], 
#		["T", dict1[element]["Tcount"]], 
#		["C", dict1[element]["Ccount"]], 
#		["G", dict1[element]["Gcount"]]], 
#		)
#	np1= np.append (np1, np2, axis=0)                  		
#print(np1)



#					R/RPY2 GGPLOT2

import rpy2.robjects.packages as packages
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as robjects
import rpy2.robjects.lib.grid as grid
R = robjects.r

#-----------------------------------------------------------------------
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


# CONVERTING MY DF TO AN R DF			WORKING!
pandas2ri.activate()
with localconverter(robjects.default_converter + pandas2ri.converter):
	gr_df = robjects.conversion.py2rpy(gdf)
#print (r_df)


# CREATING PLOTS WITH GGPLOT2  		WORKING!	

grp = ggplot2.ggplot(gr_df)


hpp = (grp
     + ggplot2.aes_string(x='Count')
     + ggplot2.geom_histogram(col="darkgrey", fill="lightblue", binwidth=2)
     + ggplot2.facet_grid(robjects.Formula(". ~ NC"))
     + ggplot2.ggtitle("Histogram of NC contents in library")
     + ggplot2.labs(x="NC count", y="nr. spots")
     + ggplot2.theme_light()
     )


hpp.plot()
robjects.r.ggsave(plot=hpp, filename="{}/Histogram.png".format(foldername), width=20, height=3)
# alternative saving:
#R("dev.copy(png,'/home/franziska/Scripts/Rtests/0209FacetHistogram.png')")



bpp = (grp 
	+ ggplot2.aes_string(x="NC", y="Count", fill="NC")
	+ ggplot2.geom_boxplot()
	+ ggplot2.ggtitle("NC contents in library")
	+ ggplot2.theme_light()
	+ ggplot2.theme(legend_position="none")
	+ ggplot2.scale_fill_brewer(palette="BuPu")
	)
	
bpp.plot()
projectname="{}_Boxplot.png".format(i)
#print (projectname)
robjects.r.ggsave(plot=bpp, filename="{}/Boxplot.png".format(foldername))


print("END. Find your Outputs in the folder stated above.")

quit()
# WORKING?
# last plot is never printed.... ? dev.copy problem. it messes up all the orders and stuff. 
# find a way to save them without saving the device copy.
# or do something akin to plt.close(plot), so that it gets deleted from memory before doing the next one. 

with localconverter(robjects.default_converter + pandas2ri.converter):
	r_df = robjects.conversion.py2rpy(spotresult)

gp = ggplot2.ggplot(r_df)

tpp = (gp
     + ggplot2.aes_string(x='Tcount')
     + ggplot2.geom_histogram()
     )

tpp.plot()
R("dev.copy(png,'/home/franziska/Scripts/Rtests/TcountHistogram.png')")

gpp = (gp
     + ggplot2.aes_string(x='Gcount')
     + ggplot2.geom_histogram()
     )

gpp.plot()
R("dev.copy(png,'/home/franziska/Scripts/Rtests/GcountHistogram.png')")

cpp = (gp
     + ggplot2.aes_string(x='Ccount')
     + ggplot2.geom_histogram()
     )

cpp.plot()
R("dev.copy(png,'/home/franziska/Scripts/Rtests/CcountHistogram.png')")


