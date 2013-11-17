import FWCore.ParameterSet.Config as cms
import sys

from TauAnalysis.SVfit.tools.getInstanceName import getInstanceName
from TauAnalysis.SVfit.tools.composeModuleName import composeModuleName
from TauAnalysis.SVfit.tools.recursiveSetAttr import recursiveSetAttr

#--------------------------------------------------------------------------------
# utility function for generation of a set of modules
# producing combinations of leptonic and hadronic decay products
# (pairs of pat::Electrons, pat::Muons and pat::Taus)
# for estimation of systematic uncertainties
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

class objProdConfigurator(cms._ParameterTypeBase):

    def __init__(self, objProd, systematics = None, pyModuleName = None):
        self.objProd = objProd
        self.systematics = systematics
        self.pyModuleName = pyModuleName,
        self.sequence = cms.Sequence()

    def _addModule(self, objProdItem, sysName, sysAttributes, pyNameSpace = None, process = None):
        # create module
        moduleType = objProdItem.type_()
        module = cms.EDProducer(moduleType)

        # set module attributes
        # to default values
        for objProdAttrName in dir(objProdItem):
            objProdAttr = getattr(objProdItem, objProdAttrName)
            if isinstance(objProdAttr, cms._ParameterTypeBase) and not objProdAttrName in [ "pluginName", "pluginType" ]:
                if isinstance(objProdAttr, cms.PSet):
                    # CV: need to clone configuration parameters of type cms.PSet,...
                    #     in order to avoid that recursiveSetAttr function
                    #     overwrites objProd "template" object passed to objProdConfigurator constructor !!
                    setattr(module, objProdAttrName, objProdAttr.clone())
                else:
                    setattr(module, objProdAttrName, objProdAttr)

        # set names of source collections
        # to objects shifted in energy/transverse momentum, theta, phi...
        for sysAttrName, sysAttrValue in sysAttributes.items():
            recursiveSetAttr(module, sysAttrName, sysAttrValue)
                
        moduleName = composeModuleName([ getInstanceName(objProdItem, pyNameSpace, process), sysName ])
        #print "moduleName = %s" % moduleName
        module.setLabel(moduleName)

        # if process object exists, attach module to process object;
        # else register module in global python name-space
        if process is not None:
            setattr(process, moduleName, module)
        else:
            pyModule = sys.modules[self.pyModuleName[0]]
            if pyModule is None:
                raise ValueError("'pyModuleName' Parameter invalid !!")
            setattr(pyModule, moduleName, module)

        # add module to sequence
        self.sequence += module

    def configure(self, pyNameSpace = None, process = None):
        # configure set of modules
        # producing different combinations of leptonic and hadronic decay products
        # for estimation of systematic uncertainties

        # add original module (for production of central value) to sequence
        self.sequence += self.objProd

        if self.systematics is not None:
            for sysName, sysAttributes in self.systematics.items():
                #print "sysName = %s" % sysName
                self._addModule(self.objProd, sysName = sysName, sysAttributes = sysAttributes,
                                pyNameSpace = pyNameSpace, process = process)

        return self.sequence
