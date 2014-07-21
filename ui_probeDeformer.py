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

try:
    cmds.loadPlugin("probeDeformer.py")
    cmds.loadPlugin("ProbeDeformer")
    cmds.loadPlugin("ProbeDeformerARAP")
    cmds.loadPlugin("probeLocator.py")
except:
    print("Plugin already loaded")


## prepare interface
class UI_ProbeDeformer:
    uiID = "ProbeDeformer"
    title = "ProbeDeformerPlugin"
    
    deformers = [0,1,2]
    probes = []
    
    ## Constructor
    def __init__(self):
        if pm.window(self.uiID, exists=True):
            pm.deleteUI(self.uiID)
        win = pm.window(self.uiID, title=self.title, menuBar=True)
        with win:
            pm.menu( label='Create', tearOff=True )
            pm.menuItem( label="Probe", c=pm.Callback( self.initPlugin, "probe") )
            pm.menuItem( label="Probe ARAP", c=pm.Callback( self.initPlugin, "probeDeformerARAP") )
            pm.menuItem( label="Probe Python", c=pm.Callback( self.initPlugin, "probePy") )
            pm.menuItem( label="Probe Locator", c=pm.Callback( self.createProbeLocator ) )
            self._parentLayout = pm.columnLayout( adj=True )
            with self._parentLayout:
                self.createUISet()
    
    def createUISet(self):
        self._childLayout = pm.columnLayout( adj=True )
        with self._childLayout:
            self.deformers[0] = pm.ls(type="probe")
            self.probes = [pm.listConnections(self.deformers[0][i].pm) for i in range(len(self.deformers[0]))]
            for i in range(len(self.deformers[0])):
                frameLayout = pm.frameLayout( label=self.deformers[0][i].name(), collapsable = True)
                with frameLayout:
                    with pm.rowLayout(numberOfColumns=6) :
                        pm.text(l='Weight Curve R')
                        pm.gradientControl( at='%s.wcr' % self.deformers[0][i].name() )
                        pm.text(l='S')
                        pm.gradientControl( at='%s.wcs' % self.deformers[0][i].name() )
                        pm.text(l='L')
                        pm.gradientControl( at='%s.wcl' % self.deformers[0][i].name() )
                    with pm.rowLayout(numberOfColumns=len(self.probes[i])+4) :
                        pm.button( l="Del", c=pm.Callback( self.deleteNode, 0, i))
                        pm.attrControlGrp( label="blend mode", attribute= self.deformers[0][i].bm)
                        pm.button( l="Add", c=pm.Callback( self.addProbe, 0, i) )
                        pm.text(l="Press to reset")
                        for j in range(len(self.probes[i])):
                            pm.button( l=self.probes[i][j].name(), c=pm.Callback( self.resetProbe, 0, i, j) )
                    with pm.rowLayout(numberOfColumns=4) :
                        pm.attrControlGrp( label="world mode", attribute= self.deformers[0][i].worldMode)
                        pm.attrControlGrp( label="rotation consistency", attribute= self.deformers[0][i].rc)
                        pm.attrControlGrp( label="Frechet sum", attribute= self.deformers[0][i].fs)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="Weight mode", attribute= self.deformers[0][i].wtm)
                        pm.attrFieldSliderGrp(label="effect radius", min=0.001, max=20.0, attribute=self.deformers[0][i].md)
                        pm.attrControlGrp( label="normExponent", attribute=self.deformers[0][i].ne)
            # probe ARAP
            self.deformers[1] = pm.ls(type="probeDeformerARAP")
            self.probes = [pm.listConnections(self.deformers[1][i].pm) for i in range(len(self.deformers[1]))]
            for i in range(len(self.deformers[1])):
                frameLayout = pm.frameLayout( label=self.deformers[1][i].name(), collapsable = True)
                with frameLayout:
                    with pm.rowLayout(numberOfColumns=6) :
                        pm.text(l='Weight Curve R')
                        pm.gradientControl( at='%s.wcr' % self.deformers[1][i].name() )
                        pm.text(l='S')
                        pm.gradientControl( at='%s.wcs' % self.deformers[1][i].name() )
                        pm.text(l='L')
                        pm.gradientControl( at='%s.wcl' % self.deformers[1][i].name() )
                    with pm.rowLayout(numberOfColumns=len(self.probes[i])+4) :
                        pm.button( l="Del", c=pm.Callback( self.deleteNode, 1, i))
                        pm.attrControlGrp( label="blend mode", attribute= self.deformers[1][i].bm)
                        pm.button( l="Add", c=pm.Callback( self.addProbe, 1, i) )
                        pm.text(l="Click to delete")
                        for j in range(len(self.probes[i])):
                            pm.button( l=self.probes[i][j].name(), c=pm.Callback( self.deleteProbe, 1, i, j) )
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="world mode", attribute= self.deformers[1][i].worldMode)
                        pm.attrControlGrp( label="rotation consistency", attribute= self.deformers[1][i].rc)
                        pm.attrControlGrp( label="normExponent", attribute=self.deformers[1][i].ne)
                    with pm.rowLayout(numberOfColumns=2) :
                        pm.attrControlGrp( label="translationWeight", attribute=self.deformers[1][i].tw)
                        pm.attrControlGrp( label="constraintWeight", attribute=self.deformers[1][i].cw)
                    with pm.rowLayout(numberOfColumns=2) :
                        pm.attrControlGrp( label="Weight mode", attribute= self.deformers[1][i].wtm)
                        pm.attrFieldSliderGrp(label="effect radius", min=0.001, max=20.0, attribute=self.deformers[1][i].md)
            # probe Python
            self.deformers[2] = pm.ls(type="probePy")
            self.probes = [pm.listConnections(self.deformers[2][i].pm) for i in range(len(self.deformers[2]))]
            for i in range(len(self.deformers[2])):
                frameLayout = pm.frameLayout( label=self.deformers[2][i].name(), collapsable = True)
                with frameLayout:
                    with pm.rowLayout(numberOfColumns=6) :
                        pm.text(l='Weight Curve R')
                        pm.gradientControl( at='%s.wcr' % self.deformers[2][i].name() )
                        pm.text(l='S')
                        pm.gradientControl( at='%s.wcs' % self.deformers[2][i].name() )
                        pm.text(l='L')
                        pm.gradientControl( at='%s.wcl' % self.deformers[2][i].name() )
                    with pm.rowLayout(numberOfColumns=len(self.probes[i])+3) :
                        pm.button( l="Del", c=pm.Callback( self.deleteNode, 2, i ))
                        pm.button( l="Add", c=pm.Callback( self.addProbe, 2, i) )
                        pm.text(l="Press to reset")
                        for j in range(len(self.probes[i])):
                            pm.button( l=self.probes[i][j].name(), c=pm.Callback( self.resetProbe, 2, i, j) )
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="blend mode", attribute= self.deformers[2][i].bm)
                        pm.attrControlGrp( label="trans mode", attribute= self.deformers[2][i].tm)
                        pm.attrControlGrp( label="rotation consistency", attribute= self.deformers[2][i].rc)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="Weight mode", attribute= self.deformers[2][i].wtm)
                        pm.attrFieldSliderGrp(label="effect radius", min=0.001, max=20.0, attribute=self.deformers[2][i].md)
                        pm.attrControlGrp( label="normExponent", attribute=self.deformers[2][i].ne)
    
    # create deformer node and connection
    def initPlugin(self, deformerType):
        # get transform nodes for the selected objects
        meshes = pm.selected(tr=1)
        if not meshes:
            return
        pm.select( meshes[-1])       # the deformer is attached to the last selected object
        deformer = pm.ls(cmds.deformer(type=deformerType)[0])[0]
        for i in range(len(meshes)-1):
            cmds.connectAttr(meshes[i]+".worldMatrix", deformer+".pm[%s]" %(i))  # connect current matrices
        for i in range(len(meshes)-1):    # to avoid uncessary ARAP computation, we divide into two for loops
            deformer.ipm[i].set(meshes[i].worldMatrix.get()) # set initial matrices
        cmds.makePaintable(deformerType, 'weights', attrType='multiFloat', shapeMode='deformer')
        self.updateUI()

    # delete deformer node
    def deleteNode(self,type,i):
        cmds.delete(self.deformers[type][i].name())
        self.updateUI()

    # set the selected probes's current matrix as the initial matrix
    def resetProbe(self,type,i,j):
        self.probes = [pm.listConnections(self.deformers[type][i].pm) for i in range(len(self.deformers[type]))]
        self.deformers[type][i].ipm[j].set(self.probes[i][j].worldMatrix.get())
    
    # set the selected probes's current matrix as the initial matrix
    def deleteProbe(self,type,i,j):
        self.probes = [pm.listConnections(self.deformers[type][i].pm) for i in range(len(self.deformers[type]))]
    #TODO

    # add a new probe
    def addProbe(self,type,i):
        self.probes = [pm.listConnections(self.deformers[type][i].pm) for i in range(len(self.deformers[type]))]
        numOldProbes = len(self.probes[i])
        newProbes = pm.selected(tr=1)
        self.probes[i].extend(newProbes)
        for j in range(numOldProbes,len(self.probes[i])):
            self.deformers[type][i].ipm[j].set(self.probes[i][j].worldMatrix.get())
        for j in range(numOldProbes,len(self.probes[i])):
            cmds.connectAttr(self.probes[i][j]+".worldMatrix", self.deformers[type][i]+".pm[%s]" %(j))
        self.updateUI()
    
    # create new probe locator
    def createProbeLocator(self):
        cmds.createNode('probeLocator')
    
    # redraw UI
    def updateUI(self):
        pm.deleteUI( self._childLayout )
        pm.setParent(self._parentLayout)
        self.createUISet()