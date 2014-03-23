
import matplotlib.pyplot as plt
import os


# normalizes bins by dividing the the log change by the number of genes in the bin

def normalizeBins(scatterPlotDict, binSize=None):
    if not binSize == None:
        bins = dict.fromkeys(scatterPlotDict.keys())
        #rename dict keys
        keys  = scatterPlotDict.keys()
        for key in keys:
            scatterPlotDict[key+', bin size = '+str(binSize)] = scatterPlotDict.pop(key)
        #makes a bins dict with all little dicts inside
        for key in scatterPlotDict:
            if not scatterPlotDict[key]['x'] == []:
                chrMax = max(scatterPlotDict[key]['x'])
                binNumber = chrMax/binSize
                bins[key]=dict.fromkeys(range(binNumber+1),0)
        #tallies things within dicts
        for key in scatterPlotDict:
            for base in scatterPlotDict[key]['x']:
                inBin = base/binSize
                bins[key][inBin] +=1
        #alteres scatterPlotDict
        for key in scatterPlotDict:
            for entry in range(len(scatterPlotDict[key]['x'])):
                inBin = scatterPlotDict[key]['x'][entry]/binSize
                scatterPlotDict[key]['y'][entry]/=bins[key][inBin]
            

def scatterPlot(scatterDict):
    try: os.makedirs(F[:-4])
    except: pass
    for key in scatterPlotDict:
        if len(scatterPlotDict[key]['x'])>0:
            print 'plotting ' + key
            plt.figure(1)
            
            plt.scatter(scatterPlotDict[key]['x'],scatterPlotDict[key]['y'])
            
            plt.title('chromosome ' + key)
            plt.xlabel('Distance from telomere')
            plt.ylabel('log2 fold change')
            plt.savefig(F[:-4] + '/scatter chromosome ' + key + '.png')
            plt.show()
            plt.close(1)

def scatterPlotLinear(dataDict, hg):
    
    scatterPlotDict = dict.fromkeys(dataDict.keys())
    for key in dataDict:
        scatterPlotDict[key] = {'x':[],'y':[]}
        for i in dataDict[key]:
            if not i[11] == 'none':
                scatterPlotDict[key]['x'].append(int(i[-2]))
                scatterPlotDict[key]['y'].append(float(i[5]))

    
    try: os.makedirs(F[:-4])
    except: pass
    for key in scatterPlotDict:
        if len(scatterPlotDict[key]['x'])>0:
            print 'plotting ' + key
            plt.figure(1)
            try: cent1,cent2 = hg[key]['centromere']
            except:  cent1,cent2 = 0,0
            telomere1 = hg[key]['telomere']
            #cetntromere line
            plt.plot([cent1, cent2], [1, 1], color='0.75', linestyle='-', linewidth=4)
            #telomere lines
            for tel in telomere1:
                plt.plot([tel, tel], [-1, 1], color='m', linestyle='-', linewidth=4)

            plt.scatter(scatterPlotDict[key]['x'],scatterPlotDict[key]['y'])
            
            plt.title('chromosome ' + key)
            plt.xlabel('Base pair')
            plt.ylabel('log2 fold change')
            plt.savefig(F[:-4] + '/scatter linear chromosome ' + key + '.png')
            plt.close(1)    
    pass
def readHG():
	f = open('mm10.hg','r')
	hg = f.read().split('\n')[2:-2]
	f.close()	
	chromosomes = {}
	for line in hg:
		parts =  line.split('\t')
		CHR = parts[1][3:]
		start = int(parts[2])
		end = int(parts[3])
		type = parts[-2]
		#print(CHR + '\t' + str(start) + '\t' + str(end) + '\t' + type)
		if not(CHR in chromosomes): 
			chromosomes.update({CHR:{'telomere':[],'centromere':[]}})
		chromosomes[CHR][type].append(start)
		chromosomes[CHR][type].append(end)
	for CHR in chromosomes:
		chromosomes[CHR]['telomere']=sorted(chromosomes[CHR]['telomere'])[1:-1]
	return chromosomes

def findClosest(what,CHR,where,hg):
	locations = hg[CHR][what]
	distances = []
	if locations == []:
		return 'none'
	for d in locations:
		distances.append(abs(d-where))
	dist=str(sorted(distances)[0])
	return dist

hg = readHG()


def filterPadj(nestedTable,P=0.05):
    for i in range(len(nestedTable)-1,-1,-1):
        try:
            Padj=float(nestedTable[i][7])
            if Padj>P:
                nestedTable.pop(i)
        except:
            nestedTable.pop(i)
            
F = 'Control_vs_AL_30_months_full_results_genes_named.txt'
f = open(F,'r')
lines = f.read().split('\n')
f.close()
head = lines[0].split('\t')
lines = lines[1:-1]

#adds chromosome, gene mid
head = head + ['gene_mid']
dataDict = {}

for f in lines:
    f=f.split('\t')
    if not f[5] == 'Inf':
        locus = f[-1]
        chromosome = locus.split(':')[0]
        if not chromosome in dataDict: dataDict.update({chromosome:[]})
        
        start = int(locus.split(':')[1].split('-')[0])
        end = int(int(locus.split(':')[1].split('-')[1]))
        midGene = str((start + end)/2)
        f.append(midGene)
        dataDict[chromosome].append(f)
for key in dataDict:
    filterPadj(dataDict[key],0.05)
    

#by this stage I have a table (newLines) with chromosome, and mid gene which low adjusted P values filtered

# now to add telomeric distances
head.append('tel_dist')
for key in dataDict:
    for i in dataDict[key]:
        try:
            i.append(findClosest('telomere',key,int(i[-1]),hg))
        except: i.append('none')
#scaps only important info for scatter plot
scatterPlotDict = dict.fromkeys(dataDict.keys())
for key in dataDict:
    scatterPlotDict[key] = {'x':[],'y':[]}
    for i in dataDict[key]:
        if not i[11] == 'none':
            scatterPlotDict[key]['x'].append(int(i[11]))
            scatterPlotDict[key]['y'].append(float(i[5]))




#scatterPlotLinear(dataDict,hg)

normalizeBins(scatterPlotDict,1000000)
scatterPlot(scatterPlotDict)
