;+
; :Author: AHACKETT
;-

photonicWaveFunction = READ_ASCII('EPCGP/EPCGP/EPCGP_2D_spin_psicA.dat')
excitonicWaveFunction = READ_ASCII('EPCGP/EPCGP/EPCGP_2D_spin_psixA.dat')

photonicAmp = REFORM(photonicWaveFunction.FIELD001[1,*])
excitonicAmp = REFORM(excitonicWaveFunction.FIELD001[1,*])

;photonicAmp = CONGRID(photonicAmp, 10000)
;excitonicAmp = CONGRID(excitonicAmp, 10000)

consitWaveAmpPlot1 = PLOT(photonicAmp, excitonicAmp,$
   TITLE = 'Plot of Relative Photonic and Excitonic Contributions to Polariton Wavefunction',$
   XTITLE = 'Photon Contribution', YTITLE = 'Exciton Contribution')

END