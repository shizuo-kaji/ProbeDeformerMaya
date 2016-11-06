# -*- coding: utf-8 -*-

#  User interface for ProbeDeformer plugin
#  The last object in the selected ones will be the target while the others will serve as "probes."

#  @author      Shizuo KAJI
#  @date        2013/5/13

# for debug
#import debugmaya
#debugmaya.startDebug()

# Import Maya Modules
import maya.cmds as cmds
import pymel.core as pm

#deformerTypes = ["probeDeformer","probeDeformerARAP","probeDeformerPy","probeLocator"]
deformerTypes = ["probeDeformer","probeDeformerARAP","probeLocator"]

for type in deformerTypes:
    try:
        cmds.loadPlugin(type)
    except:
        print("Plugin %s already loaded" %(type))

## prepare interface
class UI_ProbeDeformer:
    uiID = "ProbeDeformer"
    title = "ProbeDeformerPlugin"
    deformers = []
    probes = {}

    ## Constructor
    def __init__(self):
        if pm.window(self.uiID, exists=True):
            pm.deleteUI(self.uiID)
        win = pm.window(self.uiID, title=self.title, menuBar=True)
        with win:
            pm.menu( label='Create', tearOff=True )
            for type in deformerTypes:
                pm.menuItem( label=type, c=pm.Callback( self.initPlugin, type) )
            self._parentLayout = pm.columnLayout( adj=True )
            with self._parentLayout:
                self.createUISet()

    def createUISet(self):
        self._childLayout = pm.columnLayout( adj=True )
        with self._childLayout:
            self.deformers = [pm.ls(type=deformerTypes[i]) for i in range(len(deformerTypes))]
            for i in range(len(deformerTypes)):
                for node in self.deformers[i]:
                    self.probes[node] = pm.listConnections(node.pm)
            # "probeDeformer" specific
            for node in self.deformers[0]:
                frameLayout = pm.frameLayout( label=node.name(), collapsable = True)
                with frameLayout:
                    self.createRamp(node)
                    self.createCommonAttr(node, deformerTypes[0])
                    indices = cmds.getAttr(node+".pm", multiIndices=True)
                    if indices:
                        for j in indices:
                            with pm.rowLayout(numberOfColumns=1) :
                                pm.attrFieldSliderGrp(label=node.prw[j].getAlias(), min=0, max=10.0, attribute=node.prw[j])
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="Frechet sum", attribute= node.fs)
                        pm.attrControlGrp( label="visualisation", attribute= node.vm)
                        pm.attrFieldSliderGrp( label="visualisation multiplier", min=0.001, max=1000, attribute=node.vmp)

            # "probeDeformerARAP" specific
            for node in self.deformers[1]:
                frameLayout = pm.frameLayout( label=node.name(), collapsable = True)
                with frameLayout:
                    self.createRamp(node)
                    self.createCommonAttr(node, deformerTypes[1])
                    indices = cmds.getAttr(node+".pm", multiIndices=True)
                    if indices:
                        for j in indices:
                            with pm.rowLayout(numberOfColumns=2) :
                                pm.attrFieldSliderGrp(label=node.prw[j].getAlias(), min=0, max=1.0, attribute=node.prw[j])
                                pm.attrFieldSliderGrp(label=node.prcr[j].getAlias(), min=0, max=1.0, attribute=node.prcr[j])
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.button( l="Set supervisor", c=pm.Callback( self.setSupervisor, node))
                        pm.attrControlGrp( label="tet mode", attribute= node.tm)
                        pm.attrFieldSliderGrp( label="translation weight", min=0.0, max=1.0, attribute=node.tw)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="constraint mode", attribute= node.ctm)
                        pm.attrFieldSliderGrp( label="constraint weight", min=0.001, max=1000, attribute=node.cw)
                        pm.attrFieldSliderGrp(label="constraint radius", min=0.001, max=10.0, attribute=node.cr)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrFieldSliderGrp( label="iteration", min=1, max=20, attribute=node.it)
                        pm.attrControlGrp( label="visualisation", attribute= node.vm)
                        pm.attrFieldSliderGrp( label="visualisation multiplier", min=0.001, max=1000, attribute=node.vmp)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="stiffness mode", attribute=node.stiffnessMode)

            # "probeDeformerPy" specific
#            for node in self.deformers[2]:
#                frameLayout = pm.frameLayout( label=node.name(), collapsable = True)
#                with frameLayout:
#                    self.createRamp(node)
#                    self.createCommonAttr(node, deformerTypes[2])
#                    with pm.rowLayout(numberOfColumns=1) :
#                        pm.attrControlGrp( label="translation mode", attribute= node.tm)

    # create deformer node and connection
    def initPlugin(self, deformerType):
        if deformerType=="probeLocator":
            cmds.createNode('probeLocator')
            return
        # get transform nodes for the selected objects
        transforms = pm.selected(tr=1)
        if not transforms:
            return
        pm.select( transforms[-1])       # the deformer is attached to the last selected object
        node = pm.ls(cmds.deformer(type=deformerType)[0])[0]
        cmds.makePaintable(deformerType, 'weights', attrType='multiFloat', shapeMode='deformer')
        if len(transforms)>1:
            self.addProbe(node,deformerType,transforms[:-1])
        self.updateUI()

    # add selected transform as a new probe
    def addProbe(self,node,deformerType,newProbes):
        indexes = cmds.getAttr(node+".pm", multiIndices=True)
        if not indexes:
            n=0
        else:
            n=indexes[-1]+1
        # connect pm first to avoid unnecessary arap computations
        for j in range(len(newProbes)):
            cmds.connectAttr(newProbes[j]+".worldMatrix", node+".pm[%s]" %(j+n))
            if deformerType=="probeDeformerARAP" or deformerType=="probeDeformer":
                pm.aliasAttr(newProbes[j].name()+"_weight%s" %(j+n), node.prw[j+n].name())
            if deformerType=="probeDeformerARAP":
                pm.aliasAttr(newProbes[j].name()+"_constraintRadius%s" %(j+n), node.prcr[j+n].name())
        for j in range(len(newProbes)):
            node.ipm[j+n].set(newProbes[j].worldMatrix.get())

    # add selected transform as a new probe
    def addSelectedProbe(self,node,deformerType):
        newProbes = pm.selected(tr=1)
        self.addProbe(node,deformerType,newProbes)
        self.updateUI()

    # delete deformer node
    def deleteNode(self,node):
        cmds.delete(node.name())
        self.updateUI()

    # set selected shapes as supervised mesh
    def setSupervisor(self,node):
        meshes = pm.selected(tr=1)
        if not meshes:
            return
        for i in range(len(meshes)):
            shape=meshes[i].getShapes()[0]
            cmds.connectAttr(shape+".outMesh", node.name()+".supervisedMesh[%s]" %(i), force=True)
        self.updateUI()

    # delete a probe
    def deleteProbe(self,node,j):
        cmds.disconnectAttr(self.probes[node][j]+".worldMatrix", node+".pm[%s]" %(j) )
        self.updateUI()

    # redraw UI
    def updateUI(self):
        pm.deleteUI( self._childLayout )
        pm.setParent(self._parentLayout)
        self.createUISet()

    # create common attributes
    def createRamp(self,node):
        with pm.rowLayout(numberOfColumns=6) :
            pm.text(l='Weight Curve R')
            pm.gradientControl( at='%s.wcr' % node.name() )
            pm.text(l='S')
            pm.gradientControl( at='%s.wcs' % node.name() )
            pm.text(l='L')
            pm.gradientControl( at='%s.wcl' % node.name() )

    def createCommonAttr(self,node,deformerType):
        with pm.rowLayout(numberOfColumns=len(self.probes[node])+2) :
            pm.button( l="Delete deformer", c=pm.Callback( self.deleteNode, node))
            pm.button( l="Add selection to probes", c=pm.Callback( self.addSelectedProbe, node, deformerType) )
            for j in range(len(self.probes[node])):
                pm.button( l=self.probes[node][j].name(), c=pm.Callback( self.deleteProbe, node, j) )
        with pm.rowLayout(numberOfColumns=5) :
            pm.attrControlGrp( label="blend mode", attribute= node.bm)
            pm.attrControlGrp( label="world mode", attribute= node.worldMode)
            pm.attrControlGrp( label="rotation consistency", attribute= node.rc)
            pm.attrControlGrp( label="area weight", attribute= node.aw)
            pm.attrControlGrp( label="neighbour weighting", attribute= node.nghbrw)
        with pm.rowLayout(numberOfColumns=4) :
            pm.attrControlGrp( label="Weight mode", attribute= node.wtm)
            pm.attrFieldSliderGrp(label="effect radius", min=0.001, max=20.0, attribute=node.er)
            pm.attrControlGrp( label="normalise weight", attribute= node.nw)
            pm.attrControlGrp( label="normExponent", attribute=node.ne)

