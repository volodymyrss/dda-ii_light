import ddosa
import dataanalysis as da
from dataanalysis import hashtools, graphtools
import astropy.io.fits as fits
import pilton
import re
from numpy import *

class TimeBin(ddosa.DataAnalysis):
    tbin=2

    def get_version(self):
        return self.get_signature()+"."+self.version+".%.5lg"%self.tbin
        


class ii_light(ddosa.DataAnalysis):
    input_events=ddosa.ISGRIEvents
    input_scw=ddosa.ScWData
    input_bins=ddosa.ImageBins
    input_spectra=ddosa.ii_spectra_extract
    input_maps=ddosa.BinMapsSpectra
    input_gti=ddosa.ibis_gti
    input_dead=ddosa.ibis_dead

    input_binning=TimeBin

    cached=True

    def get_version(self):
        return self.get_signature()+"."+self.version+".%.5lgs"%self.input_binning.tbin

    def main(self):
        ddosa.construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_events.events.get_path(), \
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
            self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]", \
            self.input_gti.output_gti.get_path(), \
            self.input_dead.output_dead.get_path(), \
        ]) # get separately tc etc

        ddosa.import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND','TSTART','TSTOP','RA_SCX','DEC_SCX','RA_SCZ','DEC_SCZ'])
        ddosa.set_attr({'ISDCLEVL':"BIN_S"})

        lc_fn="lcr.fits"
        lc_tpl="(ISGR-SRC.-LCR-IDX.tpl)"
        ddosa.remove_withtemplate(lc_fn+lc_tpl)

        bin="ii_light"
        ht=ddosa.heatool(bin)
        ht['inSwg'] = "og.fits"
        ht['num_e'] = len(self.input_bins.bins)
        ht['e_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
        ht['e_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
        ht['outLC'] = lc_fn+lc_tpl
        ht['GTIname']="MERGED_ISGRI"
        ht['context']=self.input_scw.revdirpath+"/idx/isgri_context_index.fits[1]"
        ht['idxSwitch']=self.input_scw.revdirpath+"/idx/isgri_pxlswtch_index.fits[1]"
        ht['idxNoise']=self.input_scw.revdirpath+"/idx/isgri_prp_noise_index.fits[1]"
        ht['backDol']=self.input_maps.back.get_path()
        ht['corrDol']=self.input_maps.corr.get_path()
        ht['pifDOL']=self.input_spectra.pifs.get_path()
        ht['source_selectDol']=""
        ht['onlydet']="no"
        ht['chatter']=5
        ht['delta_t']=self.input_binning.tbin

        try:
            ht.run()
        except pilton.HEAToolException as e:
            r=re.search("(terminating with status -25801)",ht.output)
            if r:
                print("this code makes sense",r.groups())
            else:
                raise
        else:
            self.lc=da.DataFile(lc_fn)

            for e in fits.open(lc_fn)[2:]:
                e1=e.header['E_MIN']
                e2=e.header['E_MAX']
                name=e.header['NAME']

                savetxt("lc_%.5lg_%s_%.5lg_%.5lg.txt"%(self.input_binning.tbin,name.replace(" ","_"),e1,e2),e.data)

class PowerSpectrum(da.DataAnalysis):
    input_lc=ii_light

    def main(self):
        lc_fn=self.input_lc.lc.get_path()
        for e in fits.open(lc_fn)[2:]:
            e1=e.header['E_MIN']
            e2=e.header['E_MAX']
            name=e.header['NAME']

            r=e.data['RATE']
            r[isnan(r)]=0
            #r[isnan(r)]=average(r[~isnan(r)])
            ps=abs(fft.fft(r))**2
            freqs=fft.fftfreq(r.size, self.input_lc.tbin)

            savetxt("powerspec_%.5lg_%s_%.5lg_%.5lg.txt"%(self.input_lc.tbin,name.replace(" ","_"),e1,e2),column_stack((freqs,ps)))


class Factorize_ii_light(graphtools.Factorize):
    root='ii_light'
    leaves=["ScWData","Revolution"]

class ScWIILCList(ddosa.DataAnalysis):
    input_scwlist=ddosa.RevScWList
    copy_cached_input=False
    input_iilcsummary=Factorize_ii_light

    allow_alias=True

    version="v0"
    
    maxspec=None

    def main(self):
        self.lcs=[ii_light(assume=scw) for scw in self.input_scwlist.scwlistdata]

        if len(self.lcs)==0:
            raise ddosa.EmptyScWList()
        

class ISGRIIILCSum(ddosa.DataAnalysis):
    input_iilclist=ScWIILCList

    copy_cached_input=False

    cached=True

    version="v3"

    def main(self):
        total={}

        for lc in self.input_iilclist.lcs:
            if not hasattr(lc,'lc'):
                print("empty")
                continue

            print(lc.lc.get_path())
            f=fits.open(lc.lc.get_path())

            for i in range(len(f)-2):
                e=f[2+i]
                k=e.header['NAME'],(e.header['E_MIN'],e.header['E_MAX'])

                if k not in total:
                    total[k]=e
                else:
                    total[k].data=concatenate((total[k].data,e.data))

        for k,e in total.items():
            total_fn="total_ii_light_%s_%.5lg_%.5lg.fits"%(k[0],k[1][0],k[1][1])
            e.writeto(total_fn,overwrite=True)
            setattr(self,total_fn.replace(".fits",""),da.DataFile(total_fn))




#import dataanalysis.callback
#dataanalysis.callback.default_callback_filter.set_callback_accepted_classes([ddosa.mosaic_ii_skyimage, ddosa.ii_skyimage, ddosa.BinEventsImage, ddosa.ibis_gti, ddosa.ibis_dead, ddosa.ISGRIEvents, ddosa.ii_spectra_extract, ddosa.BinEventsSpectra, ddosa.ii_lc_extract, ddosa.BinEventsLC, ISGRISpectraSum])


