import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaElectronProducers.defaultLowPtGsfElectronIDExtra_cfi import defaultLowPtGsfElectronIDExtra

lowPtGsfElectronIDExtra = defaultLowPtGsfElectronIDExtra.clone(
    electrons = 'slimmedLowPtElectrons',
    rho = 'fixedGridRhoFastjetAll',
    ModelNames = cms.vstring(['']),
    ModelWeights = cms.vstring([
            'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root',
            ]),
    ModelThresholds = cms.vdouble([-99.])
    )
