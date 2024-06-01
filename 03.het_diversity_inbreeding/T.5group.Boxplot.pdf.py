
import matplotlib
##
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

plt.switch_backend('agg')



matplotlib.rcParams['font.sans-serif'] = ['Arial']  
matplotlib.font_manager._rebuild()    
#font10 = FontProperties(family='Arial', size=10) 
font8 = matplotlib.font_manager.FontProperties(fname='/usr/share/fonts/arial/Arial.ttf',family='Arial', size=8) 
font6 = matplotlib.font_manager.FontProperties(fname='/usr/share/fonts/arial/Arial.ttf',family='Arial', size=6) 




raw=[]
all_data =[]
colors=[]
pops=[]

popsLabel={}
popsDict={}

#regionOrderDict={"Sumatra":1,"SouthAsian":2,"SoutheastAsian":3,"EastAsian":4,"European":5}
regionOrderDict={"South_Asia":2,"MSEA":3,"East_Asia":4,"Europe":5,"ISEA":6}



filename = open('test.res4edit.txt')
for line in filename:
    lines=line.strip("\n").split("\t")
    
    value=1-float(lines[3]) / float(lines[5]) #Observed Heterozygosity   
    #print value,int(lines[3]),int(lines[5])
    #print lines[22]
    
    #value=float(lines[6])#Inbreeding Coefficient
    
    if lines[22] == "keep":
        print lines[1]
        if lines[18] in popsDict.keys():
            
            popsDict[lines[18]].append(value)
        else:
            popsDict[lines[18]]=[value]
            #print lines[17]
            regionNum=regionOrderDict[lines[18]]
            popsLabel[lines[18]]=[lines[18],lines[18],lines[19],lines[20],value,regionNum] #cut -f18-21 -> China_Bamei     EastAsian       blue    #7719AA
filename.close()





for i in popsDict.keys():
    #print(popsDict[i])
    medianV=np.nanmedian(popsDict[i])
    #print(i,popsLabel[i],medianV,popsDict[i])
    temp=popsLabel[i]
    temp.append(medianV)
    temp.append(popsDict[i])
    raw.append(temp)
    
raw.sort(key=lambda s:(s[5],-s[6]))   #Observed Heterozygosity

#raw.sort(key=lambda s:(s[5],s[6]))   #Inbreeding Coefficient

#print(raw)


for element in raw:
    col=element[3]
    pop=element[0]
    #meanL=element[10]
    medianL=element[7]
    colors.append(col)
    pops.append(pop)
    all_data.append(medianL)   



numBoxes=len(all_data)
print(numBoxes)
 
fig = plt.figure(figsize=(3, 3))

bplot = plt.boxplot(all_data,whis=1.5,
                    showfliers=False,
                    widths=0.7,
                    whiskerprops=dict(linestyle='dashed', color='black', linewidth=1),
                    boxprops=dict(linestyle='solid', color='black', linewidth=0),
                    medianprops=dict(linestyle='solid',  linewidth=1),
                    notch=False,  # notch shape
                    #vert=False,  # vertical box aligmnent
                    patch_artist=True,
                    zorder=20
                    )  # fill with color ,alpha=0.5

for i in range(numBoxes):
    y = all_data[i]
    #x = [i+1]*len(y)
    
    #rand=(np.random.rand()-0.5)/2
    x=(np.random.rand(1,len(y))-0.5)/3 + (i+1)
    x=x.tolist()[0]
    #print(x)
    plt.plot(x, y, 'r.',color='black', alpha=1,ms=1,zorder=30)

for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    
plt.xticks([x + 1 for x in range(len(all_data))], pops,fontproperties=font8) #rotation=90,

plt.ylabel('Observed Heterozygosity', fontproperties=font8) #font={'family':'Arial', 'size':16}

#plt.ylabel('Inbreeding Coefficient', fontproperties=font8) #font={'family':'Arial', 'size':16}


#color=colors,

#/data/xydu/01.work/pig-public-sra/745/S5.sec.745Her.het

fig.savefig("Cjh5groups.Heterozygosity.autosome.Boxplot.3-3.pdf", bbox_inches='tight') 

#inbreeding coefficient
#fig.savefig("Cjh5groups.Inbreeding.Coefficient.autosome.Boxplot.3-3.pdf", bbox_inches='tight') 
 


#blue    #7719AA
#green   #FF4365
#magenta #FAC05E
#red     #3FA7D6  

#s/#green\t#FF4365/gold\t#DAA520/g
    
