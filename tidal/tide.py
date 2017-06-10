import sys
sys.path.append('../')
from parsers.parserSELAFIN import SELAFIN
from runSELAFIN import scanSELAFIN
import os, glob
from datetime import datetime, timedelta
import numpy as np
from pyproj import Proj, transform
import matplotlib.pyplot as plt

class VALIDOBJECT():
    '''Reading and plotting, validating Telemac results files'''
    def __init__(self,resfile,**kwargs):
        self.resfile = resfile
        self.useropts = kwargs        
        #first we create a new selfin object
        fullpath = self.useropts['respath'] + self.resfile
        if os.path.isfile(fullpath):
            newscan = scanSELAFIN(fullpath)          
            newscan.printHeader()
            self.slf = SELAFIN(fullpath) 
            self.variables = self.checkVar()            
        else:
            raise IOError('%s does not exist. Check file name and path' % self.resfile)

    def checkVar(self):
        variables = {}
        keys = []
        values = []
        for i,name in enumerate(self.slf.VARNAMES):
            keys.append(name)
            values.append(i)
            variables[name.strip()] = i
            print "%i   %s " %(i, name )        
        print(variables['FREE SURFACE'])
        if self.useropts['variable'] not in variables:
            raise ValueError('%s does not exist in the results file' % self.useropts['variable'])
        return variables


    def plotTS(self):
        if self.useropts['locations'][0]=='all':
            print("Plotting all locations where a model node exists within %.1f metres of the tide gauge." % self.useropts['sampledist'])
            self.plotAll()   
            
    def plotAll(self):
        #Get the number of timesteps in the results file
        nsteps = 0
        for i,time in enumerate(self.slf.tags["times"]):
            nsteps = nsteps + 1        
        #Get a list off all the files in the observations directory for the correct year
        year = self.useropts['startdate'].split("/")[0]
        obsfiles = glob.glob(self.useropts['obspath'] + year + "*.txt")
        xmesh,ymesh = self.getMesh()
        for obsfile in obsfiles:
            print(obsfile)
            odate,tidewl,totwl,res,longauge,latgauge,site = self.readOBS(obsfile)
            nodes = self.nearestNodes(longauge,latgauge,xmesh,ymesh)
            if len(nodes)>0:
                telarray = self.gettelArray(nodes,nsteps) 
                newodate,obsarray = self.getobsArray(odate,tidewl,nsteps)
                telstart = int(self.useropts['spinup'])
                rmsearray = self.getRmsepks(obsarray,telarray[telstart:,:])  
                minrmse = np.argmin(rmsearray)
                telarray = telarray[:,minrmse].reshape(-1,1)            
                xlocs_minor, xlabels_minor = self.getXlabels(newodate)
                fig = plt.figure()
                ax=fig.add_subplot(111)                
                x_pos = np.arange(len(obsarray))
                obs, = plt.plot(x_pos,obsarray,'r-.',linewidth=3.0,label='Observations') 
                for i in range(len(telarray[1,:])):                    
                    omwl = np.mean(obsarray)                    
                    tmwl = np.mean(telarray[telstart:,i])                    
                    dmwl = tmwl-omwl                    
                    tel, = plt.plot(telarray[telstart:,i]-dmwl,'g',linewidth=2.0,label='Model')
                plt.title(site,fontsize=22)
                ax.set_ylabel('Elevation (m)',fontsize=22)
                ax.set_xlabel('Date',fontsize=22)
                ax.set_xticks(xlocs_minor)
                ax.set_xticklabels(xlabels_minor)
                plt.legend(handles=[tel, obs])
                ax.tick_params(labelsize=18)
                ax.set_xlim([0,len(obsarray)])
                #ax.set_xlim([72,144])
                #ax.set_xlim([96,168])
                plt.show()
                plt.close() 
    
    def readOBS(self,obs):
        f = open(obs,'r')
        hdr = f.readlines()[:11]
        f.close()
        obsperiod = 15.*60.
        telperiod = self.useropts['lprintout']*self.useropts['timestep']
        ss = int(telperiod/obsperiod)
        for i in range(0,len(hdr)):
            if "Site" in hdr[i]:
                site = (hdr[i].split())[1]
                latgauge = (hdr[i+1].split())[1]
                longauge = (hdr[i+2].split())[1]
                continue
        latgauge = float(latgauge)
        longauge = float(longauge)
        outProj = Proj(init=self.useropts['meshproj'])
        inProj = Proj(init='epsg:4326')
        longauge,latgauge = transform(inProj,outProj,longauge,latgauge)        
        f = open(obs,'r')
        lines = f.readlines()[11:]
        odate = []    
        elev = []
        res = []
        lcnt = 0
        cnt = 0
        i = 0
        totwl = np.zeros(len(lines)/ss)
        res = np.zeros(len(lines)/ss)    
        for line in lines:
            if cnt==lcnt:
                linedata = line.split()
                odate.append(linedata[1] + " " + linedata[2])            
                totwlstr = linedata[3]
                lchar = (totwlstr[-1:])        
                if lchar.isdigit():
                    totwl[i] = float(totwlstr)
                
                else:
                    totwl[i] = float(totwlstr[:-1])
            
                resstr = linedata[4]
                lchar = (resstr[-1:])        
                if lchar.isdigit():
                    res[i] = float(resstr)
                else:
                    res[i] = float(resstr[:-1])
                i = i + 1
                lcnt = lcnt + ss
            cnt = cnt + 1
        f.close()
        tidewl = totwl-res
        return odate,tidewl,totwl,res,longauge,latgauge,site

    def nearestNodes(self,longauge,latgauge,xmesh,ymesh):
        nodes = []
        dist = np.zeros((len(xmesh),2))
        for i in range(len(dist)):
            dist[i,0] = i
            dist[i,1] = ((xmesh[i]-longauge)**2+(ymesh[i]-latgauge)**2)**0.5
        dist=dist[np.argsort(dist[:,-1])]
        return dist[dist[:,1]<=self.useropts['sampledist'],0].astype(int)        
        
    def getMesh(self):
        xmesh = self.slf.MESHX
        ymesh = self.slf.MESHY
        #inProj = Proj(init=self.useropts['meshproj'])
        #outProj = Proj(init='epsg:4326')
        #xlon,ylon = transform(inProj,outProj,xmesh,ymesh)
        return xmesh,ymesh

    def gettelArray(self,nodes,nsteps):
        telarray = np.zeros((nsteps,len(nodes)))
        valueNumber = self.variables[self.useropts['variable']]
        for j in range(nsteps):
            values = self.slf.getVALUES(j)
            ztri = values[valueNumber]
            telarray[j,:] = ztri[nodes]
        return telarray

    def getobsArray(self,odate,tidewl,nsteps):
        start_time = datetime.strptime(self.useropts['startdate'],'%Y/%d/%m %H:%M:%S')
        start_spin = start_time + timedelta(hours = self.useropts['spinup'])
        end_time = start_time + timedelta(seconds = self.useropts['timestep']*nsteps)
        newodate = []
        obslist = []
        for i in range(len(odate)):
            obsdate = datetime.strptime(odate[i],'%Y/%m/%d %H:%M:%S')
            timediff = (obsdate-start_spin).total_seconds()  
            if timediff==0:            
                for j in range(i,i + nsteps-int(self.useropts['spinup'])):
                    newodate.append(odate[j])
                    obslist.append(tidewl[j])
                continue
        obsarray = np.asarray(obslist,dtype=np.float32)
        return newodate,obsarray

    def getRmsepks(self,obsarray,telarray):
        rmsearray = np.zeros(len(telarray[1,:]))
        obspks = self.getPks(obsarray)
        omwl = np.mean(obsarray) 
        for i in range(len(telarray[1,:])):                               
            tmwl = np.mean(telarray)                    
            dmwl = tmwl-omwl 
            telpks = self.getPks(telarray[:,i]-dmwl)            
            maxi = np.min([len(obspks),len(telpks)])            
            diffsq = np.zeros(maxi)
            for j in range(0, maxi):                
                diffsq[j] = (obspks[j]-telpks[j])**2
            rmsearray[i] = (np.mean(diffsq))**0.5
        return rmsearray        

    def getPks(self,tsarray):
        pks = []
        for i in range(1,len(tsarray)-1):
            if tsarray[i]>=tsarray[i-1] and tsarray[i]>tsarray[i+1] and tsarray[i]>np.mean(tsarray):
                pks.append(tsarray[i])                
        return pks

    def getXlabels(self,newodate):
        xlabels = []
        xlocs = []
        for i in range(len(newodate)):
            splitdate = newodate[i].split()
            time = (splitdate[1].split(":"))[0]
            if time=='00':
                xlabels.append("0")
                xlocs.append(i)
            if time=='12':
                xlocs.append(i)
                xlocs.append(i)
                xlabels.append("12")
                xlabels.append("\n" + splitdate[0])
        return xlocs,xlabels
                
            
        
        
            
            
