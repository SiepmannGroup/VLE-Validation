#You'll need to provide the file name for the input file(ex:isobutane.inp) if it's in the current directory and the full path if it's not.
#You'll also need to provide the molecular weight(kg/mol) and the critical temperature(K).
#It will then run through the input and check for outliers for the Pressure and Z parameters.
#All the data and results will then be put into a csv file in the current directory.
#I attached an example of an input file from the VLE data as well as the generated csv file.
#If you have any questions or get any errors please let me know and I'll get those fixed right away.
#EMAIL: ahsan033@umn.edu

import re
import statistics
import numpy as np
import csv
def mainRun():
    filefound = True
    parsingThrough = True
    zbFound = True
    zmFound = True
    ztFound = True
    uPresent = True
    lPresent = True
    textAll = []
    while filefound:
        try:
            fileName = input("File name for analysis: ")
            with open(fileName, 'r') as dataFile:
                atomicWeight = float(input("Atomic Weight in g/mol: "))
                critTemp = float(input("Critical Temp in Kelvin: "))
                outputName = input("Filename for output CSV file(without extension): ")
                for line in dataFile:
                    textAll.append(line)
                #print(textAll[15])
                atomicWeight = atomicWeight/1000
                filefound = False
                trials = textAll[2].split()
                trials = int(trials[0])
                temps = getTemps(textAll,trials)
                intConvert(temps)
                vaporDensity = getVaporDensity(textAll,atomicWeight,temps,trials)
                pressures = getPressures(textAll,temps,trials)
                intConvert(pressures)
                timesPreassure(pressures)
                #upOut,downOut = findOutliers(pressures,"pressures")
                zFinal = calualteZ(pressures,temps,atomicWeight,vaporDensity)
                #findOutliers(zFinal, "z")
                #calculateZexd(temps,407.8)
                tRR = caluclateTR(temps,critTemp)
                zexd = calculateZexd(temps,critTemp)
                pexd = calculatePexd(temps,pressures)
                lnP = lnp(pressures,pexd)
                lnZ = lnz(zFinal,zexd)
                ranges,low,mid,up = findOutliersz(tRR,lnZ)
                range, lowp,midp,topp = findOutliersPressures(tRR,lnP)
                #lower,mid,upper = sortBoundaries(tRR,lnZ)
                #calculatePexd(temps,pressures)
                writeOut(temps,pressures,vaporDensity,atomicWeight,zFinal,zexd,lnZ,tRR,ranges,low,mid,up,topp,lowp,midp,trials,outputName)
            #    print(q75)
        except FileNotFoundError:
            print("File not found")

def intConvert(nums):
    for i in range(len(nums)):
        nums[i] = float(nums[i])
        #    nums.remove(nums[i])
        #print(nums[i])

def timesPreassure(nums):
    for i in range(len(nums)):
        nums[i] = nums[i]*1000

def getTemps(allText,trials):
    tf = []
    for i in range(4,trials+4):
        ti = allText[i].split()
        #print(ti[0])
        tf.append(ti[0])
    return tf

def getVaporDensity(allText,molarMass,temps,trials):
    vf = []
    for i in range(4,trials+4):
        vi = allText[i].split()
        #print(vi)
        intConvert(vi)
        vi[1] = vi[1]*(1/1000)
        vi[1] = vi[1]*1e6
        vf.append(vi[1])
    return vf


def getPressures(allText,temps,trials):
    p = []
    for i in range(4,trials+4):
        tp = allText[i].split()
        #print(tp[5])
        p.append(tp[5])
    return p
    #print(p)

def getIQR(dataSet):
    median = statistics.median(dataSet)
    q75, q25 = np.percentile(dataSet, [75 ,25], interpolation = 'midpoint')
    iqr = q75 - q25
    iqrx = iqr * 1.5
    return median,q75,q25,iqr,iqrx

def calualteZ(pressures,temp,weight,specificVapor):
    z = []
    # weight = 0.05812#kg/mol
    # specificVapor = 2.01#kh/m^3
    R = 8.3144598#J/(K mol)
    for i in range(len(pressures)):
        top = pressures[i]*weight
        RT = R * temp[i]
        bottom = specificVapor[i]*RT
        zf = top/bottom
        z.append(zf)
    return z

def caluclateTR(temps,critTemp):
    tR = []
    for i in range(len(temps)):
        tRI = temps[i]/critTemp
        tR.append(tRI)
    return tR

def calculateZexd(temps,critTemp):
    tr = caluclateTR(temps,critTemp)
    zed = []
    for i in range( len(temps)):
        full = 0.2732*np.log(1.070-tr[i]**2.921)+1.007
        zed.append(full)
    return zed

def calculatePexd(temps,pressures):
    pexd = []
    for i in range(len(temps)):
        if i < len(temps)-2:
            slopeFinal = findCCSlope(temps,pressures,i,1,2)
            c = getC(slopeFinal,temps,pressures,i)
            final = np.exp(slopeFinal*(1/(1000/temps[i]))+c)
            pexd.append(final)
        else:
            slopeFinal = findCCSlope(temps,pressures,i,-1,-2)
            c = getC(slopeFinal,temps,pressures,i)
            final = np.exp(slopeFinal*(1/(1000/temps[i]))+c)
            pexd.append(final)
    return pexd




def findCCSlope(temps,pressures,i,num1,num2):
    slope1 = np.log(pressures[i+num2]/pressures[i+num1])
    slop2a = 1/(1000/temps[i+num2])
    slop2b = 1/(1000/temps[i+num1])
    slope2 = slop2a - slop2b
    slopeFinal = slope1/slope2
    return slopeFinal

def getC(slope,temps,pressures,i):
    c = np.log(pressures[i]/1000)-slope*(1/(1000/temps[i]))
    return c

def lnp(psim,pcc):
    lnP = []
    for i in range(len(psim)):
        m = np.log((psim[i]/1000)/pcc[i])
        lnP.append(m)
    return lnP

def lnz(zsim,zexd):
    lnF = []
    for i in range(len(zsim)):
        m = np.log(zsim[i]/zexd[i])
        lnF.append(m)
    return lnF

def sortBoundaries(tR,lnz):
    lower = []
    mid = []
    upper = []
    for i in range(len(lnz)):
        if tR[i] < .7:
            lower.append(lnz[i])
        elif tR[i] > .7 and tR[i]<.85:
            mid.append(lnz[i])
        elif tR[i]>.85:
            upper.append(lnz[i])
    return lower,mid,upper

def findOutliersz(tR,lnz):
    global zbFound
    global zmFound
    global ztFound
    ranges = []
    lowerOutliers = []
    midOutliers = []
    topOutliers = []
    lower,mid,upper = sortBoundaries(tR,lnz)
    lowerTRange = 0.05
    lowerBRange = -0.02
    midTRange = 0.05
    midBRange = -.06
    upperTRange = 0.12
    upperBRange = -0.12
    ranges.extend((lowerTRange,lowerBRange,np.mean(lower),np.std(lower),midTRange,midBRange,np.mean(mid),np.std(mid),upperTRange,upperBRange,np.mean(upper),np.std(upper)))
    zbFound = False
    zmFound = False
    ztFound = False
    for i in range(len(tR)):
        if tR[i] < .7:
            if lnz[i]<lowerBRange or lnz[i]>lowerTRange:
                print("Lower outlier Detected")
                print(lnz[i])
                lowerOutliers.append(lnz[i])
                zbFound = True
            else:
                print("No lower outlier Detected")
        elif tR[i] > .7 and tR[i]<.85:
            if lnz[i]<midBRange or lnz[i]>midTRange:
                print("Mid outlier Detected")
                midOutliers.append(lnz[i])
                print(lnz[i])
                zmFound = True
            else:
                print("No mid outlier Detected")
        elif tR[i]>.85:
            if lnz[i]<upperBRange or lnz[i]>upperTRange:
                print("Upper outlier Detected")
                topOutliers.append(lnz[i])
                print(lnz[i])
                ztFound = True
            else:
                print("No upper outlier Detected")
    return ranges,lowerOutliers,midOutliers,topOutliers

def findOutliersPressures(tR,lnz):
    global uPresent
    global mPresent
    global lPresent
    ranges = []
    lowerOutliers = []
    midOutliers = []
    topOutliers = []
    lower,mid,upper = sortBoundaries(tR,lnz)
    lowerTRange = 0.03
    lowerBRange = -0.06
    midTRange = 0.7
    midBRange = -.06
    upperTRange = 0.07
    upperBRange = -0.07
    ranges.extend((lowerTRange,lowerBRange,np.mean(lower),np.std(lower),midTRange,midBRange,np.mean(mid),np.std(mid),upperTRange,upperBRange,np.mean(upper),np.std(upper)))
    uPresent = False
    mPresent = False
    lPresent = False
    for i in range(len(tR)):
        if tR[i] < .7:
            if lnz[i]<lowerBRange or lnz[i]>lowerTRange:
                print("Lower Pressure outlier Detected")
                print(lnz[i])
                lowerOutliers.append(lnz[i])
                lPresent = True
            else:
                print("No lower Pressure outlier Detected")
        elif tR[i] > .7 and tR[i]<.85:
            if lnz[i]<midBRange or lnz[i]>midTRange:
                print("Mid Pressure outlier Detected")
                midOutliers.append(lnz[i])
                print(lnz[i])
                mPresent = True
            else:
                print("No mid Pressure outlier Detected")
        elif tR[i]>.85:
            if lnz[i]<upperBRange or lnz[i]>upperTRange:
                print("Upper Pressure outlier Detected")
                topOutliers.append(lnz[i])
                print(lnz[i])
                uPresent = True
            else:
                print("No upper Pressure outlier Detected")
    return ranges,lowerOutliers,midOutliers,topOutliers


def writeOut(temps,pressures,specificVapor,weight,zsim,zexd,lnz,tR,ranges,lower,mid,upper,upOut,downOut,middleOut,trials,outputName):
    with open(outputName+'.csv', 'w') as outlier_file:
        fullList = list(zip(temps,pressures,specificVapor,zsim,zexd,lnz,tR))
        outlier_writer = csv.writer(outlier_file, delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
        writer = csv.writer(outlier_file, lineterminator='\n')
        outlier_writer.writerow(['weight(kg/mol)',str(weight)])
        outlier_writer.writerow(['Gas Constant(J/(K*mol)',str(8.3144598)])
        outlier_writer.writerow(['Trials',str(trials)])
        writer.writerow('\n')
        outlier_writer.writerow(['Temperature(k)','Pressure(Pa)','Specific Density(kg/m^3)','Zsim','Zexd','ln(zsim/zexd)','Tr'])
        for i in range(len(zsim)):
            writer.writerow(fullList[i])
        writer.writerow('\n')
        outlier_writer.writerow(['','Tr<0.7','','','0.7<Tr<0.85','','','Tr>.85','',''])
        outlier_writer.writerow(['','Mean',str(ranges[2]),'','Mean',str(ranges[6]),'','Mean',str(ranges[10])])
        outlier_writer.writerow(['','standard Deviation',str(ranges[3]),'','standard Deviation',str(ranges[7]),'','standard Deviation',str(ranges[11])])
        outlier_writer.writerow(['','Upper Bounds,Lower Bounds',str(ranges[0]),str(ranges[1]),'Upper Bounds,Lower Bounds',str(ranges[4]),str(ranges[5]),'Upper Bounds,Lower Bounds',str(ranges[8]),str(ranges[9])])
        writer.writerow('\n')
        outlier_writer.writerow(['','Outliers'])
        if zbFound:
            outlier_writer.writerow(['','Z','Lower Bound Outliers Found',lower])
        else:
            outlier_writer.writerow(['','Z','No Lower outliers found'])
        if zmFound:
            outlier_writer.writerow(['','','Middle Outliers Found',mid])
        else:
            outlier_writer.writerow(['','','No Middle Outliers Found'])
        if ztFound:
            outlier_writer.writerow(['','','Upper Outliers Found',upper])
        else:
            outlier_writer.writerow(['','','No Upper Outliers Found'])
        writer.writerow('\n')
        if uPresent:
            outlier_writer.writerow(['','P','Upper Bound Outliers Found',upOut])
        else:
            outlier_writer.writerow(['','P','No Upper outliers found'])
        if mPresent:
            outlier_writer.writerow(['','','Lower Outliers Found',middleOut])
        else:
            outlier_writer.writerow(['','','No Mid Outliers Found'])
        if lPresent:
            outlier_writer.writerow(['','','Lower Outliers Found',downOut])
        else:
            outlier_writer.writerow(['','','No Lower Outliers Found'])
        print("Data Outputted to " + outputName)

if __name__ == '__main__':
    mainRun()
