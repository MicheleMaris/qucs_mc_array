class qucs_mc_array :
   def __init__(self,nmc,associated) :
      from collections import OrderedDict
      self.clean()
      if nmc==None : return
      self.associated=associated
      if type(nmc) == type('') : 
         self.load(nmc)
      else :
         self.nmc=nmc
   def clean(self) :
      from collections import OrderedDict
      self.__info__=OrderedDict()
      self.associated=None
      self.nmc=None
      self.Freq_GHz=None
      self.Trans=None
      self.Fcent=None
      self.Fmed=None
      self.BW=None
      self.Trans=None
      self.Fcent=None
      self.Fmed=None
      self.BW=None
      self._imc=-1
      self.shape=None
   def copy(self) :
      "creates a copy of itself"
      import copy
      return copy.deepcopy(self)
   def push(self,qo) :
      import numpy as np
      import copy
      if self._imc<0 :
         self.shape=(self.nmc+1,len(qo.Freq_GHz))
         self.Freq_GHz=qo.Freq_GHz
         self.outvolts_V=np.zeros([self.nmc+1,len(qo.Freq_GHz)])
         self.Trans=np.zeros([self.nmc+1,len(qo.Freq_GHz)])
         self.Fcent=np.zeros([self.nmc+1])
         self.Fmed=np.zeros([self.nmc+1])
         self.BW=np.zeros([self.nmc+1])
      if self._imc >= self.nmc : return
      self._imc+=1
      self.outvolts_V[self._imc]=copy.deepcopy(qo.outvolts_V)
      self.Trans[self._imc]=qo.outvolts_V/(qo.outvolts_V.sum()*(qo.Freq_GHz[1]-qo.Freq_GHz[0]))
      self.Trans[self._imc]=self.Trans[self._imc]
      self.Fcent[self._imc]=(self.Trans[self._imc]*self.Freq_GHz).sum()/self.Trans[self._imc].sum()
      xc,bw=self.calc_bw(self.Trans[self._imc])
      self.Fmed[self._imc]=xc
      self.BW[self._imc]=bw
      print self._imc," stored"
   def complete(self) :
      import numpy as np
      from collections import OrderedDict
      dx=(self.Freq_GHz[1]-self.Freq_GHz[0])
      self.TransMean=self.Trans[1:].mean(axis=0)
      self.TransInf=self.Trans[1:].min(axis=0)
      self.TransSup=self.Trans[1:].max(axis=0)
      for k in ['Mean','Inf','Sup'] :
         self.__dict__['Trans'+k]=self.__dict__['Trans'+k]/(self.__dict__['Trans'+k].sum()*(self.Freq_GHz[1]-self.Freq_GHz[0]))
      try :
         self.VMean=self.outvolts_V[1:].mean(axis=0)
      except :
         self.VMean=None
      try :
         self.VInf=self.outvolts_V[1:].min(axis=0)
      except :
         self.VMean=None
      try :
         self.VSup=self.outvolts_V[1:].max(axis=0)
      except :
         self.VMean=None
      self.gamma=np.arange(-3,3.5,.5)
      self.CC=np.zeros([len(self.gamma),self.nmc+1])
      for igamma in range(len(self.gamma)) :
         print self.gamma[igamma],
         print ((self.Freq_GHz/self.Fcent[0])**self.gamma[igamma])[0],
         print ((self.Freq_GHz/self.Fcent[0])**self.gamma[igamma])[-1]
         for imc in range(self.nmc+1) :
            fgamma=(self.Freq_GHz/self.Fcent[imc])**self.gamma[igamma]
            fgamma=fgamma/(fgamma).sum()
            self.CC[igamma,imc] = (fgamma*self.Trans[:,imc]*dx).sum()/(self.Trans[:,imc]*dx).sum()
   def sample_cdf(self,x) :
      import numpy as np
      uniqx, inverse = np.unique(x, return_inverse=True)
      h=np.bincount(inverse)
      return uniqx,h.cumsum()/float(h.sum())
   def calc_bw(self,y) :
      import numpy as np
      cy=y.cumsum()
      cy=cy/cy[-1]
      xc=np.interp(0.5,cy,self.Freq_GHz)
      x1=np.interp(0.25,cy,self.Freq_GHz)
      x2=np.interp(0.75,cy,self.Freq_GHz)
      return xc,x2-x1
   def len(self) :
      return len(self)
   def __len__(self) :
      return self.nmc+1
   def keys(self):
      return self.__dict__.keys()
   def __getitem__(self,*arg) :
      if len(arg) == 0 : return
      if type(arg[0])!=type('') : return
      if len(arg) == 1 : return self.__dict__[arg[0]]
      if type(arg[1])!=type(0) : return
      if len(arg)==2 and len(arg)-1 <= self[arg[0]].ndim : return self[arg[0]][arg[1]]
      if len(arg)==3 and len(arg)-1 <= self[arg[0]].ndim : return self[arg[0]][arg[1],arg[2]]
   def BPDict(self,imc) :
      from collections import OrderedDict
      out=OrderedDict()
      if imc < 0 or imc > self.nmc : return out
      out['Chain']=str(self.__info__['fh'])+str(self.__info__['diode'])
      out['Wavelength']=self.Freq_GHz
      out['FreqGHz']=self.Freq_GHz
      out['Trans']=self.Trans[:,imc]
      return out
   def pickle(self,ofile) :
      import pickle
      pickle.dump(self.__dict__,open(ofile,'w'))
   def load(self,ofile) :
      import pickle
      self.__dict__=pickle.load(open(ofile,'r'))
   def get_from_fits(self,ofile) :
      import pyfits
      self.__info__['f']=pyfits.open(ofile)
      self.__info__['fitsfile']=ofile
      self.__info__['h1']=self.__info__['f'][1].header.copy()
      self.__info__['h2']=self.__info__['f'][2].header.copy()
      self.__info__['fh']=self.__info__['f'][1].header['fh']
      self.__info__['diode']=self.__info__['f'][1].header['diode']
      self.__info__['omt']=self.__info__['f'][1].header['omt']
      self.__info__['fem']=self.__info__['f'][1].header['fem']
      self.__info__['bem']=self.__info__['f'][1].header['bem']
      self.__info__['phaser']=self.__info__['f'][1].header['phaser']
      self.__info__['nmc']=self.__info__['f'][1].header['nmc']
      self.nmc=self.__info__['nmc']
      self.Freq_GHz = self.__info__['f'][1].data['Freq_ghz']
      self.outvolts_V = self.__info__['f'][1].data['outvolts_V']
      self.Trans = self.__info__['f'][1].data['Trans']
      self.Fmed=self.__info__['f'][2].data['Fmed']
      self.Fcent=self.__info__['f'][2].data['Fcent']
      self.BandWidth=self.__info__['f'][2].data['BW']
   def cdf1d(self,arg) :
      """Returns the 1d cdf for the given arg
         the first sample is excluded
         output two columns: values, frequency
      """
      return sample_cdf(self.__dict__[arg][1:])

