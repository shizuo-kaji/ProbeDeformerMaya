# -*- coding: utf-8 -*-

# "Probe" locator plugin for Maya

import sys
import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import maya.OpenMayaRender as OpenMayaRender

nodeName    = "probeLocator"
id      = OpenMaya.MTypeId(0x00000111)
gl = OpenMayaRender.MHardwareRenderer.theRenderer().glFunctionTable()

class probeLocator(OpenMayaMPx.MPxLocatorNode):
    type	 	= OpenMaya.MObject()
    def __init__(self):
        OpenMayaMPx.MPxLocatorNode.__init__(self)
    
    def draw(self, view, path, style, status):
        view.beginGL()
        gl.glBegin(OpenMayaRender.MGL_LINES)
        gl.glColor3f(1, 0, 0)
        gl.glVertex3f(-1.0, 0.0, 0.0)
        gl.glVertex3f(1.0, 0.0, 0.0)
        gl.glColor3f(0, 1, 0)
        gl.glVertex3f(0.0, -1.0, 0.0)
        gl.glVertex3f(0.0, 1.0, 0.0)
        gl.glColor3f(0, 0, 1)
        gl.glVertex3f(0.0, 0.0, -1.0)
        gl.glVertex3f(0.0, 0.0, 1.0)
        gl.glEnd()
        view.endGL()

def nodeCreator():
	return OpenMayaMPx.asMPxPtr(probeLocator())

def nodeInitializer():
    a = 0

def initializePlugin(obj):
    plugin = OpenMayaMPx.MFnPlugin(obj)
    plugin.registerNode(nodeName, id, nodeCreator, nodeInitializer, OpenMayaMPx.MPxNode.kLocatorNode)

def uninitializePlugin(obj):
    plugin = OpenMayaMPx.MFnPlugin(obj)
    plugin.deregisterNode(id)
