import ddosa
import dataanalysis as da
import pyfits
from numpy import *


class ii_light(ddosa.DataAnalysis):
    input_events=ddosa.ISGRIEvents
    input_scw=ddosa.ScWData
    input_bins=ddosa.ImageBins
    input_spectra=ddosa.ii_spectra_extract
    input_maps=ddosa.BinMapsSpectra
    input_gti=ddosa.ibis_gti
    input_dead=ddosa.ibis_dead

    tbin=2

    cached=True

    def get_version(self):
        return self.get_signature()+"."+self.version+".%.5lgs"%self.tbin

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
        ht['delta_t']=self.tbin

        ht.run()

        self.lc=da.DataFile(lc_fn)
        

        for e in pyfits.open(lc_fn)[2:]:
            e1=e.header['E_MIN']
            e2=e.header['E_MAX']
            name=e.header['NAME']

            savetxt("lc_%.5lg_%s_%.5lg_%.5lg.txt"%(self.tbin,name.replace(" ","_"),e1,e2),e.data)

class PowerSpectrum(da.DataAnalysis):
    input_lc=ii_light

    def main(self):
        lc_fn=self.input_lc.lc.get_path()
        for e in pyfits.open(lc_fn)[2:]:
            e1=e.header['E_MIN']
            e2=e.header['E_MAX']
            name=e.header['NAME']

            r=e.data['RATE']
            r[isnan(r)]=0
            #r[isnan(r)]=average(r[~isnan(r)])
            ps=abs(fft.fft(r))**2
            freqs=fft.fftfreq(r.size, self.input_lc.tbin)

            savetxt("powerspec_%.5lg_%s_%.5lg_%.5lg.txt"%(self.input_lc.tbin,name.replace(" ","_"),e1,e2),column_stack((freqs,ps)))

