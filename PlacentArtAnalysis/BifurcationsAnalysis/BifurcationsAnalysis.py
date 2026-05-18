####################################################################################
######################## MASTER 2 INTERNSHIP KNEIB ANTOINE #########################
################# 3D SLICER PLUGIN FOR PLACENTA VASCULAR ANALYSIS ##################
####################################################################################
import logging
import os
import csv 
from typing import Annotated, Optional

import vtk
import qt 
import slicer
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.parameterNodeWrapper import (
    parameterNodeWrapper,
    WithinRange,
)

from slicer import vtkMRMLScalarVolumeNode
try:
    import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry # Import of the vmtk library (import vmtk don't work with slicer) @AKNE
    import vtkvmtkMiscPython as vtkvmtkMisc # Part of VMTK script library @AKNE 
except ImportError: 
    print('Error importing VMTK-related libraries')
import uuid 
try:
    import vtkSegmentationCorePython as vtkSegmentationCore
    import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
    import vtkvmtkMiscPython as vtkvmtkMisc
except ImportError:
    print('Error importing VMTK-related libraries')
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import slicer
from slicer import vtkMRMLMarkupsNode
import time

class RANSACConnection:
    def __init__(self, parent):
        self.parent = parent
        self.startingPointNode = None
        self.isPlacingPoints = False
        self.points = []
        
        # Nous aurons besoin d'importer les modules RANSAC
        # Assurez-vous que ces imports fonctionnent dans votre environnement
        try:
            from ransac_slicer.ransac import run_ransac
            from ransac_slicer.volume import Volume
            self.run_ransac = run_ransac
            self.Volume = Volume
            self.ransac_available = True
        except ImportError:
            print("Modules RANSAC non disponibles")
            self.ransac_available = False
    
    def initialize(self):
        """Initialiser les nœuds nécessaires pour le placement des points"""
        self.startingPointNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
        self.startingPointNode.SetName("RANSACPoints")
        self.startingPointNode.CreateDefaultDisplayNodes()
        displayNode = self.startingPointNode.GetDisplayNode()
        displayNode.SetSelectedColor(1.0, 0.0, 0.0)  # Rouge
        
        # Ajouter un observateur pour détecter l'ajout de points
        self.pointAddedObserverTag = self.startingPointNode.AddObserver(
            slicer.vtkMRMLMarkupsNode.PointPositionDefinedEvent, 
            self.onPointAdded
        )
    
    def startPlacePoint(self):
        """Démarrer le placement de points"""
        if not self.startingPointNode:
            self.initialize()
        
        self.isPlacingPoints = True
        
        # Activer le mode de placement des points
        slicer.app.applicationLogic().GetSelectionNode().SetActivePlaceNodeID(
            self.startingPointNode.GetID()
        )
        slicer.modules.markups.logic().StartPlaceMode(1)
    
    def onPointAdded(self, caller, event):
        """Callback appelé lorsqu'un point est ajouté"""
        if not self.isPlacingPoints:
            return
        
        point_count = self.startingPointNode.GetNumberOfControlPoints()
        if point_count > 0:
            # Récupérer la position du point ajouté
            point_position = [0, 0, 0]
            self.startingPointNode.GetNthControlPointPosition(point_count - 1, point_position)
            self.points.append(point_position)
            
            print(f"Point ajouté à la position: {point_position}")
            
            # Si nous avons exactement 2 points, nous pouvons arrêter le placement
            if point_count == 2:
                self.stopPlacePoint()
    
    def stopPlacePoint(self):
        """Arrêter le placement de points"""
        self.isPlacingPoints = False
        slicer.modules.markups.logic().StartPlaceMode(False)
    
    def connectPoints(self, inputVolume):
        """Connecter les deux points en utilisant l'algorithme RANSAC"""
        if not self.ransac_available:
            slicer.util.errorDisplay("Les modules RANSAC ne sont pas disponibles")
            return False
        
        if len(self.points) != 2:
            slicer.util.errorDisplay("Deux points sont nécessaires pour la connexion")
            return False
        
        try:
            # Vérifier que nous utilisons bien le bon volume
            if inputVolume.GetName() != "vtkMRMLScalarVolumeNode1":
                # Chercher le volume correct dans la scène
                correctVolume = slicer.util.getNode("vtkMRMLScalarVolumeNode1")
                if correctVolume:
                    print(f"Utilisation du volume correct: {correctVolume.GetName()} au lieu de {inputVolume.GetName()}")
                    inputVolume = correctVolume
                else:
                    print(f"Volume utilisé: {inputVolume.GetName()}")
            
            # Préparation des données pour RANSAC
            vol = self.Volume.from_scalar_volume(inputVolume)
            
            starting_point = np.array(self.points[0])
            direction_point = np.array(self.points[1])
            
            # Créer une structure minimale pour GraphBranches
            from ransac_slicer.branch_tree import BranchTree
            
            class MinimalBranchTree:
                def __init__(self):
                    pass
                
                def insertAfterNode(self, nodeId, parentNodeId=None, becomeIntermediaryParent=False):
                    return
            
            class MinimalGraphBranches:
                def __init__(self):
                    self.branch_list = []
                    self.names = []
                    self.nodes = []
                    self.edges = []
                    self.centerline_markups = []
                    self.contour_points_markups = []
                    self.tree_widget = MinimalBranchTree()
                    
                def create_new_branch(self, edge, branch, parent_node=None, isFromSplitBranch=False):
                    self.branch_list.append(branch)
                    self.edges.append(edge)
                    new_name = "b" + str(len(self.edges))
                    self.names.append(new_name)
                    self.centerline_markups.append(None)
                    self.contour_points_markups.append(None)
                    return new_name
                
                def create_new_markup(self, name):
                    return
                
                def update_markup(self, branch_idx):
                    return
                    
            simple_graph_branches = MinimalGraphBranches()
            
            # Paramètres RANSAC ajustés
            params = {
                'vol': vol,
                'starting_point': starting_point,
                'direction_point': direction_point,
                'starting_radius': 1,  # À ajuster selon le diamètre des vaisseaux
                'percent_inlier_points': 80.0,  # Abaissé pour être plus permissif
                'inlier_threshold': 110.0,
                'centerline_resolution': 2,
                'maximum_turn_angle': np.radians(30),  # Augmenté pour permettre des tournants plus serrés
                'min_number_of_attempts': 10000,
                'max_number_of_attempts': 50000,
                'max_number_of_cylinders': 500,  # Augmenté pour des trajets plus longs
                'smart_diameter_selection': False  # Activer l'auto-ajustement du diamètre
            }
            
            from ransac_slicer.popup_utils import CustomStatusDialog
            progress_dialog = CustomStatusDialog(
                windowTitle="Reconnexion en cours...",
                text="Veuillez patienter",
                width=300,
                height=50,
            )
            
            # Exécuter l'algorithme RANSAC
            result = self.run_ransac(**params, graph_branches=simple_graph_branches, progress_dialog=progress_dialog)
            
            # Visualiser les résultats
            if len(simple_graph_branches.branch_list) > 0:
                print(f"RANSAC a généré {len(simple_graph_branches.branch_list)} branches")
                
                # Créer un modèle 3D simplifié pour la ligne centrale
                modelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "RANSAC_Connection_Centerline")
                modelNode.CreateDefaultDisplayNodes()
                modelNode.GetDisplayNode().SetColor(1.0, 0.0, 0.0)  # Rouge
                
                # Créer des points pour la ligne centrale
                points = vtk.vtkPoints()
                lines = vtk.vtkCellArray()
                
                point_index = 0
                for branch_index, branch in enumerate(simple_graph_branches.branch_list):
                    for cyl_index, cylinder in enumerate(branch):
                        points.InsertPoint(point_index, cylinder.center)
                        if point_index > 0:
                            line = vtk.vtkLine()
                            line.GetPointIds().SetId(0, point_index-1)
                            line.GetPointIds().SetId(1, point_index)
                            lines.InsertNextCell(line)
                        point_index += 1
                
                polyData = vtk.vtkPolyData()
                polyData.SetPoints(points)
                polyData.SetLines(lines)
                
                # Créer un tube autour de la ligne centrale
                tubeFilter = vtk.vtkTubeFilter()
                tubeFilter.SetInputData(polyData)
                tubeFilter.SetRadius(0.5)  # Rayon du tube
                tubeFilter.SetNumberOfSides(16)
                tubeFilter.CappingOn()
                tubeFilter.Update()
                
                modelNode.SetAndObservePolyData(tubeFilter.GetOutput())
                
                slicer.util.infoDisplay("Connexion RANSAC terminée. Une ligne centrale a été créée.")
                return True
            else:
                slicer.util.warningDisplay("RANSAC n'a pas pu générer de connexion entre les points. Essayez d'ajuster les paramètres ou les points.")
                return False
            
        except Exception as e:
            slicer.util.errorDisplay(f"Erreur lors de la connexion RANSAC: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
        
    def cleanup(self):
        """Nettoyer les ressources"""
        if self.startingPointNode and hasattr(self, 'pointAddedObserverTag') and self.pointAddedObserverTag:
            self.startingPointNode.RemoveObserver(self.pointAddedObserverTag)
        
        if self.startingPointNode:
            slicer.mrmlScene.RemoveNode(self.startingPointNode)
            self.startingPointNode = None
        if hasattr(self, 'ransacConnection'):
            self.ransacConnection.cleanup()
        
        self.points = []
        self.isPlacingPoints = False
######################################################################################################################
######################################################################################################################
######################################################################################################################

# BifurcationsAnalysis

######################################################################################################################
######################################################################################################################
######################################################################################################################

class BifurcationsAnalysis(ScriptedLoadableModule): # Cette classe permet de configurer des informations sur le module et des informations
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("BifurcationsAnalysis")  # TODO: make this more human readable by adding spaces
        # TODO: set categories (folders where the module shows up in the module selector)
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "PlacentArtAnalysis")]
        self.parent.dependencies = []  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["John Doe (AnyWare Corp.)"]  # TODO: replace with "Firstname Lastname (Organization)"
        # TODO: update with short description of the module and a link to online module documentation
        # _() function marks text as translatable to other languages
        self.parent.helpText = _("""
This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#BifurcationsAnalysis">module documentation</a>.
""")
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = _("""
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""")

        # Additional initialization step after application startup is complete
        slicer.app.connect("startupCompleted()", registerSampleData)
#
# Register sample data sets in Sample Data module
#
def registerSampleData():
    """Add data sets to Sample Data module."""
    # It is always recommended to provide sample data for users to make it easy to try the module,
    # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

    import SampleData

    iconsPath = os.path.join(os.path.dirname(__file__), "Resources/Icons")

    # To ensure that the source code repository remains small (can be downloaded and installed quickly)
    # it is recommended to store data sets that are larger than a few MB in a Github release.

    # BifurcationsAnalysis1
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category="BifurcationsAnalysis",
        sampleName="BifurcationsAnalysis1",
        # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
        # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
        thumbnailFileName=os.path.join(iconsPath, "BifurcationsAnalysis1.png"),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        fileNames="BifurcationsAnalysis1.nrrd",
        # Checksum to ensure file integrity. Can be computed by this command:
        #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
        checksums="SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        # This node name will be used when the data set is loaded
        nodeNames="BifurcationsAnalysis1",
    )

    # BifurcationsAnalysis2
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category="BifurcationsAnalysis",
        sampleName="BifurcationsAnalysis2",
        thumbnailFileName=os.path.join(iconsPath, "BifurcationsAnalysis2.png"),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        fileNames="BifurcationsAnalysis2.nrrd",
        checksums="SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        # This node name will be used when the data set is loaded
        nodeNames="BifurcationsAnalysis2",
    )

#
# BifurcationsAnalysisParameterNode
#
@parameterNodeWrapper
class BifurcationsAnalysisParameterNode:
    """
    The parameters needed by module.

    inputVolume - The volume to threshold.
    imageThreshold - The value at which to threshold the input volume.
    invertThreshold - If true, will invert the threshold.
    thresholdedVolume - The output volume that will contain the thresholded volume.
    invertedVolume - The output volume that will contain the inverted thresholded volume.
    """
    inputVolume: vtkMRMLScalarVolumeNode
    imageThreshold: Annotated[float, WithinRange(-100, 500)] = 100
    invertThreshold: bool = False
    thresholdedVolume: vtkMRMLScalarVolumeNode
    invertedVolume: vtkMRMLScalarVolumeNode
    segmentationMask: vtkMRMLScalarVolumeNode

##############################################################
##############################################################
# BifurcationsAnalysisWidget #################################
##############################################################
##############################################################

class BifurcationsAnalysisWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    #"""Uses ScriptedLoadableModuleWidget base class, available at:
    #https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    #"""
    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self._parameterNode = None
        self._parameterNodeGuiTag = None
        self.placingPointForDeletion = False
        self.pointPlaced = False
        self.pointObserverTag = None 
        self.markupsNode = None
        self.centerline_associations = None
        self.seedPointArterie1 = None
        self.seedPointArterie2 = None
        self.pointPlacementCounter = 0
        self.ransacConnection = RANSACConnection(self)

    def setup(self) -> None:
        """Called when the user opens the module the first time and the widget is initialized."""
        ScriptedLoadableModuleWidget.setup(self)
        # Step 01 : Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/BifurcationsAnalysis.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)
        self.propertiesTable = self.ui.propertiesTable

        # Step 02 : Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)
        # Step 03 : Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = BifurcationsAnalysisLogic()
        # Step 04 : These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
            
        
        self.ui.maskSelector.setMRMLScene(slicer.mrmlScene)
        self.ui.maskSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onMaskSelectorChanged)   # Lien entre la sélection du mask et la fonction 
        self.ui.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onInputSelectorChanged)
        self.ui.pushButton.connect('clicked(bool)', self.onCorrectionButtonClicked) # Calcul des lignes centrales
        self.ui.centerlineSelector.setMRMLScene(slicer.mrmlScene)
        self.ui.centerlineSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onCenterlineSelectorChanged)
        self.ui.goButton.connect('clicked(bool)', self.onGoButtonClicked)
        self.ui.gooButton.connect('clicked(bool)', self.onGooButtonClicked)
        self.ui.extractDataButton.connect('clicked(bool)', self.onExtraDataClicked) # Lien avec le bouton extraction Excel
        
        self.ui.deleteCenterline.connect('clicked(bool)', self.startPlacePointProcedure)
        self.ui.condelCenterline.connect('clicked(bool)', self.onConfirmDeleteClicked)
        self.ui.correctCenterline.connect('clicked(bool)', self.onCorrectClicked)
        
        self.ui.initKissing.connect('clicked(bool)', self.initKissingClicked)
        
        self.ui.kissingButton.connect('clicked(bool)', self.onClipBranchesButtonClicked)
        self.ui.arteriaMerge.connect('clicked(bool)', self.arteriaMerge)
        self.ui.veinKissing.connect('clicked(bool)', self.onMergeLabelMapsButtonClicked)
        self.ui.arteriaKissing.connect('clicked(bool)', self.arteriaSeedPointsClicked)
        self.ui.correctionButton.connect('clicked(bool)', self.applyCorrectionClicked)

        self.ui.newSectionButton.connect('clicked(bool)', self.onPlacingPoint)
        self.ui.connectButton.connect('clicked(bool)', self.onConnectPoints)

        self.initializeParameterNode()  
    
    def initializeParameterNode(self):
        self._parameterNode = self.logic.getOrCreateParameterNode()
        if not self._parameterNode.GetNodeReference("InputVolume"):
            firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
            if firstVolumeNode:
                self._parameterNode.inputVolume = firstVolumeNode
        if not self._parameterNode.GetNodeReference("InputMask"):
            firstSegmentationNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
            if firstSegmentationNode:
                self._parameterNode.segmentationMask = firstSegmentationNode
    
    def onLoadPointsButton(self) -> None:
        csvFilePath = self.ui.csvPathLineEdit.currentPath
        if not os.path.exists(csvFilePath):
            qt.QMessageBox.critical(slicer.util.mainWindow(), 'Error', 'CSV file not found')
            return
        self.logic.loadPointsFromCSV(csvFilePath)
    
    def onMaskSelectorChanged(self, node):
        if not node: 
            print('No node selected')
            return
        segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
        segmentationNode.SetName(node.GetName() + "_Segmentation")
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(node, segmentationNode)
        segmentationNode.CreateDefaultDisplayNodes()
        segmentID = segmentationNode.GetSegmentation().GetNthSegmentID(0)
        if segmentID:
            segmentationNode.GetSegmentation().GetSegment(segmentID).SetName("Mask")
        displayNode = segmentationNode.GetDisplayNode()
        if displayNode: # TODO : A gérer l'affichage 3D du volume 
            displayNode.SetVisibility3D(True)
            displayNode.SetVisibility2DFill(True)
            displayNode.SetVisibility2DOutline(True)
            print("Affichage 3D du segment réussi")
        slicer.app.processEvents()

    def onInputSelectorChanged(self, node):
        return
    
    def onCenterlineSelectorChanged(self, node):
        if node:
            self._parameterNode.SetNodeReferenceID("CenterlineNode", node.GetID())
    
    def cleanup(self) -> None:
        """Called when the application closes and the module widget is destroyed."""
        self.removeObservers()

    def enter(self) -> None:
        """Called each time the user opens this module."""
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self) -> None:
        """Called each time the user opens a different module."""
        # This simply clears the reference to the parameter node if it exists,
        # ensuring no unintended interactions occur after the module is exited.
        if self._parameterNode:
            self._parameterNode = None  # Clear the reference to the parameter node
    pass

    def onSceneStartClose(self, caller, event) -> None:
        """Called just before the scene is closed."""
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event) -> None:
        """Called just after the scene is closed."""
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    # TODO : A Supprimer surement : 
    #def _checkCanApply(self, caller=None, event=None) -> None:
    #    if self._parameterNode and self._parameterNode.inputVolume and self._parameterNode.thresholdedVolume:
    #        self.ui.applyButton.toolTip = _("Compute output volume")
    #        self.ui.applyButton.enabled = True
    #    else:
    #        self.ui.applyButton.toolTip = _("Select input and output volume nodes")
    #        self.ui.applyButton.enabled = False

    ### A supprimer ??? : 
    def onApplyButton(self) -> None:
        """Run processing when user clicks "Apply" button."""
        with slicer.util.tryWithErrorDisplay(_("Failed to compute results."), waitCursor=True):
            # Compute output
            self.logic.process(self.ui.inputSelector.currentNode(), self.ui.outputSelector.currentNode(),
                               self.ui.imageThresholdSliderWidget.value, self.ui.invertOutputCheckBox.checked)
            if self.ui.invertedOutputSelector.currentNode():
                # If additional output volume is selected then result with inverted threshold is written there
                self.logic.process(self.ui.inputSelector.currentNode(), self.ui.invertedOutputSelector.currentNode(),
                                   self.ui.imageThresholdSliderWidget.value, not self.ui.invertOutputCheckBox.checked, showResult=False)
    
    def getPreprocessedPolyData(self, segmentId):
        segmentationNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLSegmentationNode")
        if not segmentationNode:
            raise ValueError("No segmentation node available. Ensure a segmentation node is active")
        segmentation = segmentationNode.GetSegmentation()
        if not segmentation:
            raise ValueError("No segments found in the segmentation node")
        polyData = vtk.vtkPolyData()
        success = slicer.modules.segmentations.logic().GetSegmentClosedSurfaceRepresentation(segmentationNode, segmentId, polyData)
        if not success or polyData.GetNumberOfPoints() == 0:
            raise ValueError("Valid input surface is required for the segment.")
        return polyData
    
    # Après avoir sélectionné Original Data + Mask on appuie sur le bouton et ça affiche les lignes centrales : 
    def onCorrectionButtonClicked(self):
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        inputVolumeNode = self.ui.inputSelector.currentNode()
        maskNode = self.ui.maskSelector.currentNode()
        centerlineNode = self.ui.centerlineSelector.currentNode()

        if not inputVolumeNode or not maskNode:
            qt.QApplication.restoreOverrideCursor()
            qt.QMessageBox.critical(slicer.util.mainWindow(), 'Error', 'Please select both the input volume and the segmentation mask.')
            return
        
        qt.QMessageBox.information(slicer.util.mainWindow(), "Extraction", "Please be patient, it will only take a few moments.")
        print("Starting Preprocessing of polydata")
        segmentationNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLSegmentationNode")
        segmentation = segmentationNode.GetSegmentation()
        try: 
            slicer.util.showStatusMessage(_("Preprocessing..."))
            slicer.app.processEvents()

            if not centerlineNode:
                segmentationNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLSegmentationNode")
                if not segmentationNode:
                    raise ValueError("Segmentation node not found.")
                print("Segmentation node found.")

                segmentation = segmentationNode.GetSegmentation()
                segmentIDs = vtk.vtkStringArray()
                segmentation.GetSegmentIDs(segmentIDs)
                print(f"Found {segmentIDs.GetNumberOfValues()} segments.")

                endPointsMarkupNode = self._parameterNode.GetNodeReference("EndPoints")
                if not endPointsMarkupNode:
                    endPointsMarkupNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", slicer.mrmlScene.GetUniqueNameByString("Centerline endpoints"))
                    endPointsMarkupNode.CreateDefaultDisplayNodes()
                    self._parameterNode.SetNodeReferenceID("EndPoints", endPointsMarkupNode.GetID())
                self.endPointsMarkupNode = endPointsMarkupNode

                startPointPosition = None
                print("Endpoints markup node prepared.")
                for i in range(segmentIDs.GetNumberOfValues()):
                    segmentId = segmentIDs.GetValue(i)
                
                    preprocessedPolyData = self.getPreprocessedPolyData(segmentId)
                    print(f"Processing segment {segmentId}")
                    print("Starting network extraction.")
                    networkPolyData = self.logic.extractNetwork(preprocessedPolyData, endPointsMarkupNode)
                    print("Network extraction completed.")

                    startPositionIndex = self.logic.startPointIndexFromEndPointsMarkupsNode(endPointsMarkupNode)
                    if startPositionIndex >= 0:
                        startPointPosition = [0.0, 0.0, 0.0]
                        endPointsMarkupNode.GetNthControlPointPosition(startPositionIndex, startPointPosition)
                        print(f"Start point position for segment {segmentId} set to: {startPointPosition}")
                    else:
                        print(f"No start point position found for segment {segmentId}")
                    if startPointPosition is None:
                        startPointPosition = [0.0, 0.0, 0.0]
                    endpointPositions = self.logic.getEndPoints(networkPolyData, startPointPosition)
                    if endpointPositions:
                        for position in endpointPositions:
                            endPointsMarkupNode.AddControlPoint(vtk.vtkVector3d(position))
                        print(f"Endpoints added to markups node for segment {segmentId}.")
                    else:
                        print(f"No endpoints detected for segment {segmentId}.")
            
                endPointsMarkupNode.GetDisplayNode().PointLabelsVisibilityOff()
            
                networkModelNode = self._parameterNode.GetNodeReference("NetworkModel")
                if not networkModelNode:
                    networkModelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", slicer.mrmlScene.GetUniqueNameByString("NetworkModel"))
                    networkModelNode.CreateDefaultDisplayNodes()
                    self._parameterNode.SetNodeReferenceID("NetworkModel", networkModelNode.GetID())
                
                #networkCurveNode = self._parameterNode.GetNodeReference("NetworkCurve")
                networkPropertiesTableNode = self._parameterNode.GetNodeReference("NetworkProperties")
                if networkModelNode or networkPropertiesTableNode:
                    slicer.util.showStatusMessage(_("Extract network..."))
                    slicer.app.processEvents()  # force update
                    networkPolyData = self.logic.extractNetwork(preprocessedPolyData, endPointsMarkupNode, computeGeometry=True)
                if networkModelNode:
                    networkModelNode.SetAndObserveMesh(networkPolyData)
                    if not networkModelNode.GetDisplayNode():
                        networkModelNode.CreateDefaultDisplayNodes()
                        networkModelNode.GetDisplayNode().SetColor(0.0, 0.0, 1.0)
                        segmentation.GetDisplayNode().SetOpacity(0.4)
                #if networkCurveNode:
                    #self.logic.addNetworkCurves(networkPolyData, networkCurveNode)
                if networkPropertiesTableNode:
                    self.logic.addNetworkProperties(networkPolyData, networkPropertiesTableNode)
                print("Network model extraction completed.")
            else:
                networkPolyData = centerlineNode.GetPolyData()
                qt.QMessageBox.information(slicer.util.mainWindow(), "Everything selected", "Now that you've selected everything, you can begin propagation...")

            # Afficher un message avant de passer à l'extraction des branches
            #qt.QMessageBox.information(slicer.util.mainWindow(), "Branches Extraction", "Starting branch extraction. Please wait, it's almost finish...")
            #slicer.util.showStatusMessage(_("Pausing before branch extraction..."))             # Ajouter une pause avant de continuer
            #slicer.app.processEvents()  # force update
            #qt.QTimer.singleShot(2000, lambda: self.extractBranchesAndDisplay(networkPolyData))
            #slicer.util.showStatusMessage(_("Branch extraction almost done..."))
            #slicer.app.processEvents()
            #branchPolyData = self.logic.extractBranches(networkPolyData)
            #if branchPolyData:
            #    branchModelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "NetworkBranches")
            #    branchModelNode.SetAndObserveMesh(branchPolyData)
            #    branchModelNode.CreateDefaultDisplayNodes()
            #    branchModelNode.GetDisplayNode().SetColor(0.0, 1.0, 0.0)  # Green for branches
            #    self._parameterNode.SetNodeReferenceID("BranchModel", branchModelNode.GetID())

                networkProperties = self.logic.extractNetworkProperties(networkPolyData)
                self.fillPropertiesTable(networkProperties)
        except Exception as e:
            slicer.util.errorDisplay(str(e))
            import traceback
            traceback.print_exc()
        finally:
            qt.QApplication.restoreOverrideCursor()
            slicer.util.showStatusMessage(_("Automatic endpoint computation complete."), 3000)
    
    def fillPropertiesTable(self, networkProperties):
        """
        Fill the properties table with network properties
        :param networkProperties: A list of dictionaries with network properties
        """
        if not networkProperties:
            return

        # Exclude the 'CellId' column
        headers = [key for key in networkProperties[0].keys() if key != "CellId"]
        self.propertiesTable.setColumnCount(len(headers))
        self.propertiesTable.setHorizontalHeaderLabels(headers)

        self.propertiesTable.setRowCount(len(networkProperties))

        # Fill the table with data
        for rowIndex, properties in enumerate(networkProperties):
            colIndex = 0
            for key, value in properties.items():
                if key == "CellId":
                    continue
                item = qt.QTableWidgetItem(str(value))
                self.propertiesTable.setItem(rowIndex, colIndex, item)
                colIndex += 1

        self.propertiesTable.resizeColumnsToContents()

        for colIndex in range(len(headers)):
            if headers[colIndex] == "CenterlineId":
                self.propertiesTable.setColumnWidth(colIndex, 70)
            elif headers[colIndex] == "Length in mm":
                self.propertiesTable.setColumnWidth(colIndex, 80)
            elif headers[colIndex] == "Tortuosity":
                self.propertiesTable.setColumnWidth(colIndex, 60)
            else:
                self.propertiesTable.setColumnWidth(colIndex, 60)
        self.propertiesTable.setSortingEnabled(True)
    
    def onNetworkLoaded(self):
        networkPolyData = self.logic.loadNetwork()
        networkProperties = self.logic.extractNetworkProperties(networkPolyData)
        self.fillPropertiesTable(networkProperties)
    
    def onGoButtonClicked(self):
        centerlineNode = self.ui.centerlineSelector.currentNode()
        if centerlineNode and isinstance(centerlineNode, slicer.vtkMRMLModelNode):
            centerlinePolyData = centerlineNode.GetPolyData()
            if not centerlinePolyData:
                qt.QMessageBox.critical(slicer.util.mainWindow(), "Error", "No valid polydata in the selected centerline node.")
                return

            qt.QMessageBox.information(slicer.util.mainWindow(), "Centerline selected", "La sélection a fonctionné, l'initialisation va débuter...")

            bifurcation_points_forward = self.logic.extractBifurcationPoints(centerlinePolyData)
            bifurcation_points_reverse = self.logic.extractBifurcationPointsReverse(centerlinePolyData)

            unique_bifurcation_points = self.logic.mergeBifurcationPoints(bifurcation_points_forward, bifurcation_points_reverse)

            end_points, branch_points, centerline_associations = self.logic.classify_bifurcation_points(unique_bifurcation_points, centerlinePolyData)

            self.logic.visualizeBifurcationPoints([p["point"] for p in end_points], "EndPoints")
            self.logic.visualizeBifurcationPoints([p["point"] for p in branch_points], "BranchPoints")

            #for centerline_id, points in centerline_associations.items():
            #    print(f"Centerline ID: {centerline_id}")
            #    print(f"  End Points: {points['end_points']}")
            #    print(f"  Branch Points: {points['branch_points']}")
            ombilical_points,  centerline_radii = self.logic.detectOmbilicalPoints(centerlinePolyData, unique_bifurcation_points, centerline_associations) 
            self.logic.markOmbilicalPoints(ombilical_points, centerline_associations, centerline_radii) # Appel de la fonction de marquage de l'origine des vaisseaux ombilicaux
            
            qt.QMessageBox.information(slicer.util.mainWindow(), "Initialisation réussie", "Veuillez vérifier que les vaisseaux ont été correctement identifié avant d'appuyer sur Go...")
            
            self.centerline_associations = centerline_associations

            bifurcation_count = len(branch_points)
            self.ui.bifurcationCountLabel.setText(str(bifurcation_count))
            endpoint_count = len(end_points)
            self.ui.endpointCountLabel.setText(str(endpoint_count))
            return centerline_associations
    # Etape finale du programme qui permet de sauvegarder les données extraites, corrigées et analysés à la fin du programme : 
    def onExtraDataClicked(self): 
        return
######################################################################################################################
######################################################################################################################
    # ZONE SUPPRESSION D'UNE CENTERLINE ID 
    # Outil de suppression d'une ligne centrale, à utiliser sur de petits segments que l'on veut pas forcément extraire à la fin...
    def startPlacePointProcedure(self):
        # Create a new markup node if it does not exist
        if not hasattr(self, 'pointNode') or self.pointNode is None:
            self.pointNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            self.pointNode.SetName("Centerline a supprimer")
            self.pointNode.CreateDefaultDisplayNodes()
            self.pointNode.GetDisplayNode().SetSelectedColor(1, 0, 0)  # Red color for visibility
        slicer.app.applicationLogic().GetSelectionNode().SetActivePlaceNodeID(
            self.pointNode.GetID()
        )
        slicer.modules.markups.logic().StartPlaceMode(1)
        self.addObserver(
            self.pointNode,
            vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self.pointPlacedCallback
        )
    def pointPlacedCallback(self, caller, event):
        slicer.util.infoDisplay("Point correctement placé")
        slicer.modules.markups.logic().StartPlaceMode(False)
        self.removeObserver(
            self.pointNode,
            vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self.pointPlacedCallback
        )
    # Confirmation pour supprimer une ID 
    def onConfirmDeleteClicked(self):
        centerlineNode = self.ui.centerlineSelector.currentNode()
        if centerlineNode and isinstance(centerlineNode, slicer.vtkMRMLModelNode):
            centerlinePolyData = centerlineNode.GetPolyData()
            if not centerlinePolyData:
                qt.QMessageBox.critical(slicer.util.mainWindow(), "Error", "No valid polydata in the selected centerline node.")
                return
        
            if self.pointNode and self.pointNode.GetNumberOfControlPoints() > 0:
                control_point = [0.0, 0.0, 0.0]
                self.pointNode.GetNthControlPointPosition(0, control_point)
                nearby_centerline_ids = self.find_nearby_centerlines(centerlinePolyData, control_point, search_radius=1)
                slicer.util.infoDisplay(f"Nearby centerline IDs: {nearby_centerline_ids}")

                # Suppression des branches sélectionnées
                newPolyData = self.remove_centerline_branches(centerlinePolyData, set(nearby_centerline_ids))
                # Mise à jour du PolyData sur le nœud
                centerlinePolyData.ShallowCopy(newPolyData)
                centerlineNode.Modified()
                slicer.mrmlScene.Modified()

                self.pointNode.RemoveNthControlPoint(0)

                slicer.app.processEvents()
            
            else:
                slicer.util.errorDisplay("No control point found.")
    # Localisation de l'ID de la ligne centrale la plus proche
    def find_nearby_centerlines(self, centerline_surface, control_point, search_radius=1):
        locator = vtk.vtkPointLocator()
        locator.SetDataSet(centerline_surface)
        locator.BuildLocator()
        points_in_radius = vtk.vtkIdList()
        locator.FindPointsWithinRadius(search_radius, control_point, points_in_radius)
        centerline_ids = set()
        id_array = centerline_surface.GetCellData().GetArray("CenterlineIds")
        for i in range(points_in_radius.GetNumberOfIds()):
            point_id = points_in_radius.GetId(i)
            cell_ids = vtk.vtkIdList()
            centerline_surface.GetPointCells(point_id, cell_ids)
            for j in range(cell_ids.GetNumberOfIds()):
                cell_id = cell_ids.GetId(j)
                centerline_id = id_array.GetTuple1(cell_id)
                centerline_ids.add(centerline_id)
        return list(centerline_ids)
    def remove_centerline_branches(self, centerlinePolyData, centerline_ids_to_remove):
        import vtk
        cell_ids = vtk.vtkIdTypeArray()
        cell_ids.SetNumberOfComponents(1)
        for i in range(centerlinePolyData.GetNumberOfCells()):
            cell_id = centerlinePolyData.GetCellData().GetArray("CenterlineIds").GetTuple1(i)
            if cell_id not in centerline_ids_to_remove:
                cell_ids.InsertNextValue(i)
        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        selectionNode.SetSelectionList(cell_ids)
        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)
        extractSelection = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0, centerlinePolyData)
        extractSelection.SetInputData(1, selection)
        extractSelection.Update()
        newPolyData = vtk.vtkPolyData()
        geometryFilter = vtk.vtkGeometryFilter()
        geometryFilter.SetInputData(extractSelection.GetOutput())
        geometryFilter.Update()
        newPolyData.ShallowCopy(geometryFilter.GetOutput())
        return newPolyData

######################################################################################################################
    # Zone de correction d'une ligne centrale 
    def onCorrectClicked(self):
        return
######################################################################################################################
    # Zone de labeling des lignes centrales (il faut d'abord appuyer sur Initialization) 
    def onGooButtonClicked(self):
        if self.centerline_associations is None:
            qt.QMessageBox.critical(slicer.util.mainWindow(), "Error", "Please run the initialization first by clicking Go.")
            return
        ombilicalPointsNode = slicer.mrmlScene.GetFirstNodeByName("OmbilicalPoints")  # On récup ici les coordonnées des points 
        if not ombilicalPointsNode:
            slicer.util.errorDisplay("Ombilical Points node not found. Please ensure the points are initialized correctly.")
            return
        numberOfPoints = ombilicalPointsNode.GetNumberOfControlPoints()
        if numberOfPoints == 0:
            slicer.util.errorDisplay("No ombilical points found in the fiducial node.")
            return
        centerlineNode = self.ui.centerlineSelector.currentNode()  # On récupère la ligne centrale 
        if not centerlineNode:
            slicer.util.errorDisplay("No centerline node selected.")
            return
        centerlinePolyData = centerlineNode.GetPolyData()  # On charge le polydata de la ligne centrale 
        arteryVeinColors = vtk.vtkUnsignedCharArray()  # On crée un nouvel array pour les couleurs 
        arteryVeinColors.SetNumberOfComponents(3)  # 3 couleurs possibles : Veine, Artère et pas de couleur 
        arteryVeinColors.SetName("ArteryVeinColors")  # Le nom de l'array

        num_cells = centerlinePolyData.GetNumberOfCells()
        arteryVeinMap = {}
        initial_artery_ids = []
        initial_vein_ids = []

        print("Specific Ombilical Points Coordinates:")
        for i in range(ombilicalPointsNode.GetNumberOfControlPoints()):
            label = ombilicalPointsNode.GetNthControlPointLabel(i)
            if label in ["Veine Ombilicale", "Artère Ombilicale"]:
                pointCoord = [0.0, 0.0, 0.0]
                ombilicalPointsNode.GetNthControlPointPosition(i, pointCoord)
                print(f"{label}: {pointCoord}")
                nearby_centerline_ids = self.find_nearby_centerlines(centerlinePolyData, pointCoord, search_radius=1)  # Adjust search_radius as needed
                for cid in nearby_centerline_ids:
                    cid = int(cid)  # Assurez-vous que les ids sont des entiers
                    arteryVeinMap[cid] = (0, 0, 255) if label == "Veine Ombilicale" else (255, 0, 0)
                    if label == "Artère Ombilicale":
                        initial_artery_ids.append(cid)
                    elif label == "Veine Ombilicale":
                        initial_vein_ids.append(cid)
                print(f"Nearby Centerline IDs for {label}: {nearby_centerline_ids}")

        vein_centerline_ids = set(initial_vein_ids)
        artery_centerline_ids = set(initial_artery_ids)
        modelDisplayNode = centerlineNode.GetDisplayNode()

        for i in range(num_cells):
            cellId = int(centerlinePolyData.GetCellData().GetArray("CenterlineIds").GetTuple1(i))
            if cellId in arteryVeinMap:
                arteryVeinColors.InsertNextTuple3(*arteryVeinMap[cellId])
            else:
                arteryVeinColors.InsertNextTuple3(128, 128, 128)  # Grey for undefined
        centerlinePolyData.GetCellData().AddArray(arteryVeinColors)
        centerlinePolyData.GetCellData().SetActiveScalars("ArteryVeinColors")
        modelDisplayNode.SetScalarVisibility(True)
        centerlineNode.Modified()

        def find_branch_points_for_centerline(centerline_id):
            return self.centerline_associations.get(centerline_id, {}).get('branch_points', [])

        def propagate_alternately(vein_ids, artery_ids):
            from collections import deque
            combined_queue = deque()
            visited = set()

            vein_visited_branch_points = set()
            artery_visited_branch_points = set()

            min_length = min(len(vein_ids), len(artery_ids))
            for i in range(min_length):
                combined_queue.append((vein_ids[i], 'vein'))
                combined_queue.append((artery_ids[i], 'artery'))

            remaining_ids = vein_ids[min_length:] if len(vein_ids) > len(artery_ids) else artery_ids[min_length:]
            remaining_type = 'vein' if len(vein_ids) > len(artery_ids) else 'artery'
            for id in remaining_ids:
                combined_queue.append((id, remaining_type))

            while combined_queue:
                current_id, current_type = combined_queue.popleft()
                if current_id in visited:
                    continue
                visited.add(current_id)
                if current_type == 'vein':
                    vein_centerline_ids.add(current_id)
                else:
                    artery_centerline_ids.add(current_id)
                    
                branch_points_current = find_branch_points_for_centerline(current_id)
                for bp in branch_points_current:
                    if current_type == 'vein':
                        vein_visited_branch_points.add(bp)
                    else:
                        artery_visited_branch_points.add(bp)

                    for cid in self.centerline_associations:
                        if bp in self.centerline_associations[cid]['branch_points'] and cid not in visited:
                            combined_queue.append((cid, current_type))

                color = (255, 0, 0) if current_type == 'vein' else (0, 0, 255)
                for i in range(num_cells):
                    cellId = int(centerlinePolyData.GetCellData().GetArray("CenterlineIds").GetTuple1(i))
                    if cellId == current_id:
                        arteryVeinColors.SetTuple3(i, *color)
            centerlinePolyData.GetCellData().Modified()
            centerlineNode.Modified()

            return visited, vein_visited_branch_points, artery_visited_branch_points

        visited, vein_branch_points, artery_branch_points = propagate_alternately(list(vein_centerline_ids), list(artery_centerline_ids))

        intersections = vein_branch_points.intersection(artery_branch_points)
        print(f"Number of branch points containing both a vein ID and an artery ID: {len(intersections)}")
        print(f"Intersections found: {intersections}")

        propertiesTable = self.ui.propertiesTable
        columnCount = propertiesTable.columnCount
        propertiesTable.setColumnCount(columnCount + 1)
        propertiesTable.setHorizontalHeaderItem(columnCount, qt.QTableWidgetItem("Type"))

        for row in range(propertiesTable.rowCount):
            centerline_id = int(propertiesTable.item(row, 0).text())  # Assuming the first column has centerline IDs
            if centerline_id in vein_centerline_ids:
                propertiesTable.setItem(row, columnCount, qt.QTableWidgetItem("Vein"))
            elif centerline_id in artery_centerline_ids:
                propertiesTable.setItem(row, columnCount, qt.QTableWidgetItem("Artery"))
            else:
                propertiesTable.setItem(row, columnCount, qt.QTableWidgetItem("Unknown"))
        
        columnCount = propertiesTable.columnCount
        propertiesTable.setColumnCount(columnCount + 1)
        propertiesTable.setHorizontalHeaderItem(columnCount, qt.QTableWidgetItem("Level"))
        
        column_name = "Type"
        for colIndex in range(propertiesTable.columnCount):
            if propertiesTable.horizontalHeaderItem(colIndex).text() == column_name:
                propertiesTable.setColumnWidth(colIndex, 50)  # Définir la largeur souhaitée ici, par exemple, 100 pixels
        
        column_name = "Level"
        for colIndex in range(propertiesTable.columnCount):
            if propertiesTable.horizontalHeaderItem(colIndex).text() == column_name: 
                propertiesTable.setColumnWidth(colIndex, 50)

        print("Propagation complète.")
        
        #self.project_colors_to_surface(centerlinePolyData)
    
    def project_colors_to_surface(self, centerlinePolyData):
        if not hasattr(self, 'clippedModelNode'):
            slicer.util.errorDisplay("No clipped surface node found. Please run the clipping step first.")
            return

        surfaceNode = self.clippedModelNode
        surfacePolyData = surfaceNode.GetPolyData()

        # Utiliser vtkvmtkSurfaceProjection pour projeter les couleurs sur la surface
        surfaceProjection = vtkvmtkMisc.vtkvmtkSurfaceProjection()
        surfaceProjection.SetInputData(surfacePolyData)
        surfaceProjection.SetReferenceSurface(centerlinePolyData)
        surfaceProjection.Update()

        projectedColors = surfaceProjection.GetOutput().GetPointData().GetArray("ArteryVeinColors")

        surfacePolyData.GetPointData().AddArray(projectedColors)
        surfacePolyData.GetPointData().SetActiveScalars("ArteryVeinColors")

        surfaceNode.Modified()
        slicer.util.showStatusMessage("Projection des couleurs terminée.", 3000)

    ######################################################################################################################
    ######################################################################################################################
    # Zone de code pour la correction de kissing vessels 
    # Etape n°1 : On supprime tout ce qui doit être supprimé : donc surface et centerline
    def initKissingClicked(self):
        if not hasattr(self, 'pointNode') or self.pointNode is None:
            self.pointNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            self.pointNode.SetName("Identification Kissing")
            self.pointNode.CreateDefaultDisplayNodes()
            self.pointNode.GetDisplayNode().SetSelectedColor(1, 0, 0)  # Red color for visibility
        slicer.app.applicationLogic().GetSelectionNode().SetActivePlaceNodeID(
            self.pointNode.GetID()
        )
        slicer.modules.markups.logic().StartPlaceMode(1)
        self.addObserver(
            self.pointNode,
            vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self.onPointPlaced
        )

    def onPointPlaced(self, caller, event):
        slicer.util.infoDisplay("Point correctement placé")
        slicer.modules.markups.logic().StartPlaceMode(False)
        self.removeObserver(
            self.pointNode,
            vtkMRMLMarkupsNode.PointPositionDefinedEvent
        )
        centerlineNode = self.ui.centerlineSelector.currentNode()
        if centerlineNode and isinstance(centerlineNode, slicer.vtkMRMLModelNode):
            centerlinePolyData = centerlineNode.GetPolyData()
            if not centerlinePolyData:
                qt.QMessageBox.critical(slicer.util.mainWindow(), "Error", "No valid polydata in the selected centerline node.")
                return

            if self.pointNode and self.pointNode.GetNumberOfControlPoints() > 0:
                control_point = [0.0, 0.0, 0.0]
                self.pointNode.GetNthControlPointPosition(0, control_point)
                print(f"Control point position: {control_point}")
                nearby_centerline_ids = self.find_nearby_centerlines(centerlinePolyData, control_point, search_radius=1)
                print(f"Nearby centerline IDs: {nearby_centerline_ids}")
                slicer.util.infoDisplay(f"Nearby centerline IDs: {nearby_centerline_ids}")
                self.nearby_centerline_ids = [int(cid) for cid in nearby_centerline_ids]
                self.pointNode.RemoveNthControlPoint(0) # Permet de supprimer le point après l'avoir utilisé 
                #TODO : supprimer le branch_point associé à l'intersection

    def onClipBranchesButtonClicked(self):
        segmentationNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLSegmentationNode")
        if not segmentationNode:
            slicer.util.errorDisplay("Aucun segmentationNode trouvé dans la scène.")
            return

        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        folderItemId = shNode.CreateFolderItem(shNode.GetSceneItemID(), "Exported Models")
        logic = slicer.modules.segmentations.logic()  # Utiliser la logique du module de segmentations pour la conversion
        if not logic.ExportVisibleSegmentsToModels(segmentationNode, folderItemId):  # Convertir tous les segments visibles en modèles et les stocker dans le dossier spécifié
            slicer.util.errorDisplay("L'exportation des segments a échoué.")
            return
        folderItemChildren = vtk.vtkIdList()  # Récupérer le modèle exporté
        shNode.GetItemChildren(folderItemId, folderItemChildren)

        for i in range(folderItemChildren.GetNumberOfIds()):
            modelItemId = folderItemChildren.GetId(i)
            modelNode = shNode.GetItemDataNode(modelItemId)
            if isinstance(modelNode, slicer.vtkMRMLModelNode):
                print("Modèle exporté:", modelNode.GetName())

                centerlineNode = self.ui.centerlineSelector.currentNode()
                if not centerlineNode or not isinstance(centerlineNode, slicer.vtkMRMLModelNode):
                    slicer.util.errorDisplay("No valid centerline node selected.")
                    return

                surfacePolyData = modelNode.GetPolyData()
                centerlinePolyData = centerlineNode.GetPolyData()

                # Remplacer les GroupIds par les CenterlineIds
                cell_data = centerlinePolyData.GetCellData()
                centerline_ids = vtk.util.numpy_support.vtk_to_numpy(cell_data.GetArray("CenterlineIds"))

                if cell_data.HasArray("GroupIds"):
                    cell_data.RemoveArray("GroupIds")
                    print("Ecrasement de l'array effectué avec succès")

                new_group_ids = vtk.util.numpy_support.numpy_to_vtk(centerline_ids)
                new_group_ids.SetName("GroupIds")
                cell_data.AddArray(new_group_ids)

                # Vérifiez si les nearby_centerline_ids sont disponibles
                if not hasattr(self, 'nearby_centerline_ids') or not self.nearby_centerline_ids:
                    slicer.util.errorDisplay("No nearby centerline IDs found. Please place a point first.")
                    return

                for group_id_to_clip in self.nearby_centerline_ids:
                    print(f"Clipping with group_id_to_clip: {group_id_to_clip}")
                    clipped_surface = self.clip_branch(surfacePolyData, centerlinePolyData, [group_id_to_clip], inside_out=True)

                    if clipped_surface.GetNumberOfPoints() == 0:
                        slicer.util.errorDisplay(f"Clipping failed for group ID {group_id_to_clip}.")
                        continue

                    # Mettez à jour surfacePolyData pour le prochain clipping
                    surfacePolyData.DeepCopy(clipped_surface)
                    
                    # Supprimer les cellules associées à l'ID de la ligne centrale clippée
                    centerlinePolyData = self.remove_centerline_branches(centerlinePolyData, set([group_id_to_clip]))

                    slicer.util.infoDisplay(f"Clipping successful for group ID {group_id_to_clip}. Waiting for 5 seconds before next operation.")
                    time.sleep(5)

                # Mettez à jour le nœud de la ligne centrale après toutes les suppressions
                centerlineNode.SetAndObservePolyData(centerlinePolyData)
                centerlineNode.Modified()
                
                slicer.mrmlScene.Modified()
                clipped_model_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Clipped Surface")
                clipped_model_node.SetAndObservePolyData(surfacePolyData)
                clipped_model_node.CreateDefaultDisplayNodes()
                slicer.mrmlScene.AddNode(clipped_model_node)

                self.clippedModelNode = clipped_model_node

                slicer.app.processEvents()
                
    def clip_branch(self, surface, centerlines, group_ids, radius_array_name='Radius', inside_out=True):
        try:
            branchClipper = vtkvmtkComputationalGeometry.vtkvmtkPolyDataCenterlineGroupsClipper()
        except AttributeError:
            slicer.util.errorDisplay("La classe vtkvmtkPolyDataCenterlineGroupsClipper n'a pas été trouvée dans vtkvmtkMiscPython.")
            return

        branchClipper.SetInputData(surface)
        branchClipper.SetCenterlines(centerlines)
        branchClipper.SetCenterlineGroupIdsArrayName("GroupIds")
        branchClipper.SetGroupIdsArrayName("GroupIds")
        branchClipper.SetCenterlineRadiusArrayName(radius_array_name)
        branchClipper.SetBlankingArrayName("Blanking")
        branchClipper.SetCutoffRadiusFactor(1E16)
        branchClipper.SetClipValue(0.0)
        branchClipper.SetUseRadiusInformation(1)

        print(f"Clipping with group_ids: {group_ids} and inside_out: {inside_out}")

        if group_ids:
            groupIdsList = vtk.vtkIdList()
            for group_id in group_ids:
                groupIdsList.InsertNextId(group_id)
            branchClipper.SetCenterlineGroupIds(groupIdsList)
            branchClipper.ClipAllCenterlineGroupIdsOff()
        else:
            branchClipper.ClipAllCenterlineGroupIdsOn()

        if inside_out:
            branchClipper.GenerateClippedOutputOn()
        else:
            branchClipper.GenerateClippedOutputOff()

        branchClipper.Update()
        clipped_surface = branchClipper.GetClippedOutput() if inside_out else branchClipper.GetOutput()

        print(f"Clipped surface has {clipped_surface.GetNumberOfPoints()} points and {clipped_surface.GetNumberOfCells()} cells")
        
        return clipped_surface
    # Etape n°2 : On reconstruit ce qui doit être reconstruit grâce à la méthode colliding fronts 
    def arteriaSeedPointsClicked(self):
        self.seedPointArterie1 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
        self.seedPointArterie1.SetName("SeedPointArterie1")
        self.seedPointArterie1.CreateDefaultDisplayNodes()
        self.seedPointArterie1.GetDisplayNode().SetSelectedColor(1, 0, 0)  # Red color for visibility
        
        self.seedPointArterie2 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
        self.seedPointArterie2.SetName("SeedPointArterie2")
        self.seedPointArterie2.CreateDefaultDisplayNodes()
        self.seedPointArterie2.GetDisplayNode().SetSelectedColor(0, 1, 0)  # Green color for visibility

        slicer.app.applicationLogic().GetSelectionNode().SetActivePlaceNodeID(self.seedPointArterie1.GetID())
        slicer.modules.markups.logic().StartPlaceMode(1)

        self.addObserver(self.seedPointArterie1, vtkMRMLMarkupsNode.PointPositionDefinedEvent, self.onSeedPointPlaced)
        self.addObserver(self.seedPointArterie2, vtkMRMLMarkupsNode.PointPositionDefinedEvent, self.onSeedPointPlaced)

    def onSeedPointPlaced(self, caller, event):
        self.pointPlacementCounter += 1
        if self.pointPlacementCounter == 1:
            slicer.util.infoDisplay("First seed point placed. Now place the second seed point.")
            slicer.app.applicationLogic().GetSelectionNode().SetActivePlaceNodeID(self.seedPointArterie2.GetID())
            slicer.modules.markups.logic().StartPlaceMode(1)
        elif self.pointPlacementCounter == 2:
            slicer.util.infoDisplay("Second seed point placed.")
            slicer.modules.markups.logic().StartPlaceMode(False)
            self.removeObserver(self.seedPointArterie1, vtkMRMLMarkupsNode.PointPositionDefinedEvent)
            self.removeObserver(self.seedPointArterie2, vtkMRMLMarkupsNode.PointPositionDefinedEvent)
            self.pointPlacementCounter = 0

    def addObserver(self, node, event, callback):
        node.AddObserver(event, callback)

    def removeObserver(self, node, event):
        if node.HasObserver(event):
            node.RemoveObserver(event)

    def arteriaMerge(self):
        import vtkvmtkSegmentationPython as vtkvmtkSegmentation

        # Premièrement nous avons besoin des nodes :
        inputVolume = self.ui.inputSelector.currentNode()
        if not inputVolume:
            slicer.util.infoDisplay("Merci de sélectionner un volume valide.")
            return False

        seed1 = slicer.util.getNode('SeedPointArterie1')
        seed2 = slicer.util.getNode('SeedPointArterie2')
        if not seed1 or not seed2:
            slicer.util.infoDisplay('Merci de placer des points valides.')
            return False

        # Création d'un nouveau vtkMRMLModelNode pour la sortie
        newModelNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelNode")
        newModelNode.UnRegister(None)
        newModelNode.SetName(slicer.mrmlScene.GetUniqueNameByString("OutputModel"))
        currentModelNode = slicer.mrmlScene.AddNode(newModelNode)
        currentModelNode.CreateDefaultDisplayNodes()

        # Convertir les fiduciaires en vtkIdList
        seeds1 = self.convertFiducialHierarchyToVtkIdList(seed1, inputVolume)
        seeds2 = self.convertFiducialHierarchyToVtkIdList(seed2, inputVolume)
        stoppers = vtk.vtkIdList()

        # Combiner les deux listes de seeds
        seeds = vtk.vtkIdList()
        for i in range(seeds1.GetNumberOfIds()):
            seeds.InsertNextId(seeds1.GetId(i))
        for i in range(seeds2.GetNumberOfIds()):
            seeds.InsertNextId(seeds2.GetId(i))
        
        inputImage = vtk.vtkImageData()
        inputImage.DeepCopy(inputVolume.GetImageData())
        
        scalarRange = inputImage.GetScalarRange()
        thresholdMin = scalarRange[0]
        thresholdMax = scalarRange[1]

        initImageData = vtk.vtkImageData()
        evolImageData = vtk.vtkImageData()

        initImageData.DeepCopy(self.logic.performInitialization(inputImage,
                                                                thresholdMin,
                                                                thresholdMax,
                                                                seeds,
                                                                stoppers,
                                                                'collidingfronts'))

        if not initImageData.GetPointData().GetScalars():
            # something went wrong, the image is empty
            logging.error("Segmentation failed - the output was empty..")
            return False

        evolImageData.DeepCopy(self.logic.performEvolution(inputImage,
                                                    initImageData,
                                                    10,
                                                    0,
                                                    70,
                                                    50,
                                                    'geodesic'))

        labelMap = vtk.vtkImageData()
        labelMap.DeepCopy(self.logic.buildSimpleLabelMap(evolImageData, 5, 0))

        outputLabelMapNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", "OutputLabelMap")
        outputLabelMapNode.CopyOrientation(inputVolume)
        outputLabelMapNode.SetAndObserveImageData(labelMap)

        slicer.util.setSliceViewerLayers(background=inputVolume, label=outputLabelMapNode, labelOpacity=0.5)

        # Update the 3D model node with the segmentation result
        model = vtk.vtkPolyData()
        ijkToRasMatrix = vtk.vtkMatrix4x4()
        outputLabelMapNode.GetIJKToRASMatrix(ijkToRasMatrix)
        model.DeepCopy(self.logic.marchingCubes(evolImageData, ijkToRasMatrix, 0.0))

        currentModelNode.SetAndObservePolyData(model)
        currentModelNode.CreateDefaultDisplayNodes()
        currentModelDisplayNode = currentModelNode.GetDisplayNode()

        # Configure the display node to show the model
        currentModelDisplayNode.SetColor(1.0, 0.55, 0.4)  # Example color (red)
        currentModelDisplayNode.SetBackfaceCulling(0)
        currentModelDisplayNode.SetSliceIntersectionVisibility(0)
        currentModelDisplayNode.SetVisibility(1)
        currentModelDisplayNode.SetOpacity(1.0)

        slicer.util.infoDisplay(f"Modèle créé avec succès : {currentModelNode.GetName()}")

        # Génération de la ligne centrale qui correspond à ce segment : 
        arteriaCenterline = self.generateCenterlineForSegment(currentModelNode)

        brancharteriaCenterline = self.logic.branchextractor(arteriaCenterline)
        # Mettre à jour le modèle de NetworkModel avec le polydata obtenu
        networkModelNode = self._parameterNode.GetNodeReference("NetworkModel")
        if networkModelNode:
            networkModelNode.SetAndObservePolyData(brancharteriaCenterline)
            if not networkModelNode.GetDisplayNode():
                networkModelNode.CreateDefaultDisplayNodes()
                networkModelNode.GetDisplayNode().SetColor(0.0, 1.0, 0.0)  # Exemple de couleur (vert)
                currentModelNode.GetDisplayNode().SetOpacity(0.4)

            # Trouver les points les plus proches de SeedPointArterie1 et SeedPointArterie2
        seed1_position = [0.0, 0.0, 0.0]
        seed1.GetNthControlPointPosition(0, seed1_position)
        seed2_position = [0.0, 0.0, 0.0]
        seed2.GetNthControlPointPosition(0, seed2_position)

        closest_point_seed1 = self.find_closest_point(brancharteriaCenterline, seed1_position)
        closest_point_seed2 = self.find_closest_point(brancharteriaCenterline, seed2_position)

        # Ajouter les points trouvés comme endpoints
        self.add_point_to_markups_node(closest_point_seed1, "ClosestSeedPoint1")
        self.add_point_to_markups_node(closest_point_seed2, "ClosestSeedPoint2")
        
            # Récupérer les points de bifurcation à partir du nœud fiducial BranchPoints
        branch_points_node = slicer.util.getNode("BranchPoints")
        branch_points_list = []
        for i in range(branch_points_node.GetNumberOfControlPoints()):
            branch_point = [0.0, 0.0, 0.0]
            branch_points_node.GetNthControlPointPosition(i, branch_point)
            branch_points_list.append(branch_point)

        # Trouver les points de bifurcation les plus proches des points ClosestSeedPoint1 et ClosestSeedPoint2
        closest_branch_point_seed1 = self.find_closest_point_in_list(branch_points_list, closest_point_seed1)
        closest_branch_point_seed2 = self.find_closest_point_in_list(branch_points_list, closest_point_seed2)

        # Créer une ligne entre les points trouvés
        self.create_line_between_points_with_attributes(closest_point_seed1, closest_branch_point_seed1, networkModelNode)
        self.create_line_between_points_with_attributes(closest_point_seed2, closest_branch_point_seed2, networkModelNode)

        self.combine_network_and_branch_models()

    def generateCenterlineForSegment(self, segmentModelNode):
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            print("Starting Preprocessing of polydata")

            endPointsMarkupNode = self._parameterNode.GetNodeReference("EndPoints")
            if not endPointsMarkupNode:
                endPointsMarkupNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", slicer.mrmlScene.GetUniqueNameByString("Centerline endpoints"))
                endPointsMarkupNode.CreateDefaultDisplayNodes()
                self._parameterNode.SetNodeReferenceID("EndPoints", endPointsMarkupNode.GetID())
            self.endPointsMarkupNode = endPointsMarkupNode
           
            preprocessedPolyData = segmentModelNode.GetPolyData()

            if preprocessedPolyData is None:
                raise ValueError(f"No polydata found for segment model: {segmentModelNode.GetName()}")
        
            print(f"Processing segment {segmentModelNode.GetName()}")
            print("Starting network extraction.")
            networkPolyData = self.logic.extractNetwork(preprocessedPolyData, endPointsMarkupNode)
            print("Network extraction completed.")

            startPointPosition = [0.0, 0.0, 0.0]
            endpointPositions = self.logic.getEndPoints(networkPolyData, startPointPosition)
            if endpointPositions:
                for position in endpointPositions:
                    endPointsMarkupNode.AddControlPoint(vtk.vtkVector3d(position))
                print(f"Endpoints added to markups node for segment {segmentModelNode.GetName()}.")
            else:
                print(f"No endpoints detected for segment {segmentModelNode.GetName()}.")
            
            endPointsMarkupNode.GetDisplayNode().PointLabelsVisibilityOff()

            networkModelNode = self._parameterNode.GetNodeReference("NetworkModel")
            if not networkModelNode:
                networkModelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", slicer.mrmlScene.GetUniqueNameByString("NetworkModel"))
                networkModelNode.CreateDefaultDisplayNodes()
                self._parameterNode.SetNodeReferenceID("NetworkModel", networkModelNode.GetID())
            
            networkPropertiesTableNode = self._parameterNode.GetNodeReference("NetworkProperties")
            if networkModelNode or networkPropertiesTableNode:
                slicer.util.showStatusMessage(_("Extract network..."))
                slicer.app.processEvents()
                networkPolyData = self.logic.extractNetwork(preprocessedPolyData, endPointsMarkupNode, computeGeometry=True)
            if networkModelNode:
                networkModelNode.SetAndObserveMesh(networkPolyData)
                if not networkModelNode.GetDisplayNode():
                    networkModelNode.CreateDefaultDisplayNodes()
                    networkModelNode.GetDisplayNode().SetColor(0.0, 0.0, 1.0)
                    segmentModelNode.GetDisplayNode().SetOpacity(0.4)
            if networkPropertiesTableNode:
                self.logic.addNetworkProperties(networkPolyData, networkPropertiesTableNode)
            print("Network model extraction completed.")
            
            return networkModelNode.GetPolyData()
        
        except Exception as e:
            qt.QMessageBox.critical(slicer.util.mainWindow(), 'Error', str(e))
        finally:
            qt.QApplication.restoreOverrideCursor()

    def convertFiducialHierarchyToVtkIdList(self, hierarchyNode, volumeNode):
        outputIds = vtk.vtkIdList()

        if not hierarchyNode or not volumeNode:
            return outputIds
        if isinstance(hierarchyNode, slicer.vtkMRMLMarkupsFiducialNode) and isinstance(volumeNode, slicer.vtkMRMLScalarVolumeNode):
            image = volumeNode.GetImageData()
            # Boucle sur les fiduciaires
            for n in range(hierarchyNode.GetNumberOfFiducials()):
                currentCoordinatesRAS = [0, 0, 0]

                # Récupérer les coordonnées actuelles
                hierarchyNode.GetNthFiducialPosition(n, currentCoordinatesRAS)

                # Convertir les coordonnées RAS en IJK
                currentCoordinatesIJK = self.ConvertRAStoIJK(volumeNode, currentCoordinatesRAS)

                # Convertir les coordonnées IJK en identifiant de point
                currentCoordinatesIJKlist = (int(currentCoordinatesIJK[0]), int(currentCoordinatesIJK[1]), int(currentCoordinatesIJK[2]))
                outputIds.InsertNextId(int(image.ComputePointId(currentCoordinatesIJKlist)))

        # La liste d'ID a été créée, la retourner même si elle est vide
        return outputIds   
    @staticmethod
    def ConvertRAStoIJK(volumeNode, rasCoordinates):
        rasToIjkMatrix = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(rasToIjkMatrix)

        # Les coordonnées RAS doivent être de 4 éléments
        if len(rasCoordinates) < 4:
            rasCoordinates.append(1)

        ijkCoordinates = rasToIjkMatrix.MultiplyPoint(rasCoordinates)

        return ijkCoordinates
    
    #TODO : MAINTENANT QU'ON A LE MORCEAU, on calcule la ligne centrale, on merge et on rajoute ET C'EST TERMINE !!!!!!!!!!!
    def find_closest_point(self, polydata, target_point):
        point_locator = vtk.vtkPointLocator()
        point_locator.SetDataSet(polydata)
        point_locator.BuildLocator()
        closest_point_id = point_locator.FindClosestPoint(target_point)
        closest_point = polydata.GetPoint(closest_point_id)
        return closest_point

    def add_point_to_markups_node(self, point, node_name):
        markups_node = slicer.mrmlScene.GetFirstNodeByName(node_name)
        if not markups_node:
            markups_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", node_name)
        markups_node.AddControlPoint(vtk.vtkVector3d(point))
    
    def find_closest_point_in_list(self, points_list, target_point):
        import numpy as np
        closest_point = None
        min_distance = float('inf')
        for point in points_list:
            distance = np.linalg.norm(np.array(point) - np.array(target_point))
            if distance < min_distance:
                closest_point = point
                min_distance = distance
        return closest_point
    
    def create_line_between_points_with_attributes(self, point1, point2, networkModelNode):
        line_source = vtk.vtkLineSource()
        line_source.SetPoint1(point1)
        line_source.SetPoint2(point2)
        line_source.Update()

        line_polydata = line_source.GetOutput()

        # Copier les attributs géométriques de l'ancienne ligne centrale
        input_polydata = networkModelNode.GetPolyData()
        point_data = input_polydata.GetPointData()
        cell_data = input_polydata.GetCellData()

        new_point_data = line_polydata.GetPointData()
        new_cell_data = line_polydata.GetCellData()

        for i in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(i)
            new_array = array.NewInstance()
            new_array.DeepCopy(array)
            new_point_data.AddArray(new_array)

        for i in range(cell_data.GetNumberOfArrays()):
            array = cell_data.GetArray(i)
            new_array = array.NewInstance()
            new_array.DeepCopy(array)
            new_cell_data.AddArray(new_array)

        # Combiner la ligne existante et la nouvelle ligne avec les attributs géométriques
        append_filter = vtk.vtkAppendPolyData()
        append_filter.AddInputData(input_polydata)
        append_filter.AddInputData(line_polydata)
        append_filter.Update()

        combined_polydata = append_filter.GetOutput()
        networkModelNode.SetAndObservePolyData(combined_polydata)
        networkModelNode.Modified()

    def combine_network_and_branch_models(self):
        network_model_node = self._parameterNode.GetNodeReference("NetworkModel")
        branch_model_node = slicer.util.getNode("Placenta_77_Branchmodel")

        if not network_model_node or not branch_model_node:
            slicer.util.infoDisplay("Veuillez vérifier que les nœuds 'NetworkModel' et 'Placenta_77_Branchmodel' existent.")
            return

        network_polydata = network_model_node.GetPolyData()
        branch_polydata = branch_model_node.GetPolyData()

        if not network_polydata or not branch_polydata:
            slicer.util.infoDisplay("Les données polydata sont manquantes dans l'un des nœuds.")
            return

        combined_polydata = self.combine_polydata_with_attributes(network_polydata, branch_polydata)

        # Mettre à jour le Placenta_77_Branchmodel avec le polydata combiné
        branch_model_node.SetAndObservePolyData(combined_polydata)
        branch_model_node.Modified()
        
        slicer.util.infoDisplay("Les modèles 'NetworkModel' et 'Placenta_77_Branchmodel' ont été combinés avec succès.")

    def combine_polydata_with_attributes(self, polydata1, polydata2):
        append_filter = vtk.vtkAppendPolyData()
        append_filter.AddInputData(polydata1)
        append_filter.AddInputData(polydata2)
        append_filter.Update()

        combined_polydata = append_filter.GetOutput()

        point_data1 = polydata1.GetPointData()
        cell_data1 = polydata1.GetCellData()
        point_data2 = polydata2.GetPointData()
        cell_data2 = polydata2.GetCellData()

        combined_point_data = combined_polydata.GetPointData()
        combined_cell_data = combined_polydata.GetCellData()

        # Copier les attributs géométriques et les données des points et des cellules de polydata1
        for i in range(point_data1.GetNumberOfArrays()):
            array = point_data1.GetArray(i)
            new_array = array.NewInstance()
            new_array.DeepCopy(array)
            combined_point_data.AddArray(new_array)

        for i in range(cell_data1.GetNumberOfArrays()):
            array = cell_data1.GetArray(i)
            new_array = array.NewInstance()
            new_array.DeepCopy(array)
            combined_cell_data.AddArray(new_array)

        # Copier les attributs géométriques et les données des points et des cellules de polydata2
        for i in range(point_data2.GetNumberOfArrays()):
            array = point_data2.GetArray(i)
            new_array = array.NewInstance()
            new_array.DeepCopy(array)
            combined_point_data.AddArray(new_array)

        for i in range(cell_data2.GetNumberOfArrays()):
            array = cell_data2.GetArray(i)
            new_array = array.NewInstance()
            new_array.DeepCopy(array)
            combined_cell_data.AddArray(new_array)

        combined_point_data.Modified()
        combined_cell_data.Modified()
        combined_polydata.Modified()

        return combined_polydata
    # Zone d'application des changements :
    def applyCorrectionClicked(self):
        network_model_node = self._parameterNode.GetNodeReference("Placenta_77_Branchmodel")

        if not network_model_node:
            slicer.util.infoDisplay("Veuillez vérifier que le nœud 'Placenta_77_Branchmodel' existe.")
            return

        centerlinePolyData = network_model_node.GetPolyData()
        if not centerlinePolyData:
            slicer.util.infoDisplay("Aucune donnée polydata valide trouvée dans le nœud 'Placenta_77_Branchmodel'.")
            return

        # Recalculer les points de bifurcation et les endpoints
        bifurcation_points_forward = self.logic.extractBifurcationPoints(centerlinePolyData)
        bifurcation_points_reverse = self.logic.extractBifurcationPointsReverse(centerlinePolyData)

        unique_bifurcation_points = self.logic.mergeBifurcationPoints(bifurcation_points_forward, bifurcation_points_reverse)

        end_points, branch_points, centerline_associations = self.logic.classify_bifurcation_points(unique_bifurcation_points, centerlinePolyData)

        # Supprimer les anciens `Endpoints` et `Branchpoints` nodes
        existing_endpoints_node = slicer.mrmlScene.GetFirstNodeByName("EndPoints")
        if existing_endpoints_node:
            slicer.mrmlScene.RemoveNode(existing_endpoints_node)

        existing_branchpoints_node = slicer.mrmlScene.GetFirstNodeByName("BranchPoints")
        if existing_branchpoints_node:
            slicer.mrmlScene.RemoveNode(existing_branchpoints_node)

        # Visualiser les nouveaux points
        self.logic.visualizeBifurcationPoints([p["point"] for p in end_points], "EndPoints")
        self.logic.visualizeBifurcationPoints([p["point"] for p in branch_points], "BranchPoints")

        self.centerline_associations = centerline_associations

        slicer.util.infoDisplay("Points de bifurcation et endpoints corrigés et visualisés avec succès.")
    
    def onMergeLabelMapsButtonClicked(self):
        labelmapNames = ["LabelMap1", "LabelMap2", "LabelMap3"]
        segmentationName = "CombinedSegmentation"
        self.logic.combineMultipleLabelMapsAndCreateSegmentation(labelmapNames, segmentationName)
        return
    
    def onPlacingPoint(self):
        self.ransacConnection.startPlacePoint()

    def onConnectPoints(self):
        volumeNode = self.ui.inputSelector.currentNode()
        if not volumeNode:
            slicer.util.errorDisplay("Veuillez sélectionner un volume d'entrée")
            return
        self.ransacConnection.connectPoints(volumeNode)
######################################################################################################################
######################################################################################################################
######################################################################################################################

# BifurcationsAnalysisLogic

######################################################################################################################
######################################################################################################################
######################################################################################################################
class BifurcationsAnalysisLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self) -> None:
        """Called when the logic class is instantiated. Can be used for initializing member variables."""
        ScriptedLoadableModuleLogic.__init__(self)
        self.radiusArrayName = "Radius"
        self.topologyArrayName = "Topology"
        self.marksArrayName = "Marks"
        self.lengthArrayName = "Length"
        self.curvatureArrayName = "Curvature"
        self.torsionArrayName = "Torsion"
        self.tortuosityArrayName = "Tortuosity"
        self.frenetTangentArrayName = "FrenetTangent"
        self.frenetNormalArrayName = "FrenetNormal"
        self.frenetBinormalArrayName = "FrenetBinormal"
        
        self.blankingArrayName = 'Blanking'
        self.groupIdsArrayName = 'GroupIds'
        self.centerlineIdsArrayName = 'CenterlineIds'
        self.tractIdsArrayName = 'TractIds'

    def getParameterNode(self):
        return BifurcationsAnalysisParameterNode(super().getParameterNode())
    
    def getOrCreateParameterNode(self):
        """Ensure a parameter node exists, create if necessary, and return it."""
        node = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScriptedModuleNode")
        if not node:
            node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode", "VSegPluginParameters")
            self.setDefaultParameters(node)
        return node
    
    def setSegmentationMask(self, node):
        self.getParameterNode().segmentationMask = node

    def correctionFunction(self):
        print("Correction function executed.")
        # Ajoutez ici la logique de correction des vaisseaux kiss
        pass
    
    def extractNetwork(self, surfacePolyData, endPointsMarkupsNode, computeGeometry=False):
        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        import vtkvmtkMiscPython as vtkvmtkMisc
        

        decimate = False
        if decimate:
            decimationFilter = vtk.vtkDecimatePro()
            decimationFilter.SetInputData(surfacePolyData)
            decimationFilter.SetTargetReduction(0.99)
            decimationFilter.SetBoundaryVertexDeletion(0)
            decimationFilter.PreserveTopologyOn()
            decimationFilter.Update()

        # @AKNE Clean and triangulate
        cleaner = vtk.vtkCleanPolyData()
        if decimate:
            cleaner.SetInputData(decimationFilter.GetOutput())
        else:
            cleaner.SetInputData(surfacePolyData)
        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(cleaner.GetOutputPort())
        triangleFilter.Update()
        simplifiedPolyData = triangleFilter.GetOutput()

        print(f"Simplified polydata has {simplifiedPolyData.GetNumberOfPoints()} points")

        # Cut hole at start position
        if endPointsMarkupsNode and endPointsMarkupsNode.GetNumberOfControlPoints() > 0:
            startPosition = [0, 0, 0]
            endPointsMarkupsNode.GetNthControlPointPosition(
                self.startPointIndexFromEndPointsMarkupsNode(endPointsMarkupsNode), startPosition)
        else:
            # If no endpoints are specific then use the closest point to a corner
            bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            simplifiedPolyData.GetBounds(bounds)
            startPosition = [bounds[0], bounds[2], bounds[4]]
        self.openSurfaceAtPoint(simplifiedPolyData, startPosition)

        # @AKNE : Extract network
        networkExtraction = vtkvmtkMisc.vtkvmtkPolyDataNetworkExtraction()
        networkExtraction.SetInputData(simplifiedPolyData)
        networkExtraction.SetAdvancementRatio(1.05)
        networkExtraction.SetRadiusArrayName(self.radiusArrayName)
        networkExtraction.SetTopologyArrayName(self.topologyArrayName)
        networkExtraction.SetMarksArrayName(self.marksArrayName)
        networkExtraction.Update()
        networkOutput = networkExtraction.GetOutput()
        print(f"Extracted network has {networkOutput.GetNumberOfPoints()} points")

        if computeGeometry:
            centerlineGeometry = vtkvmtkComputationalGeometry.vtkvmtkCenterlineGeometry()
            centerlineGeometry.SetInputData(networkExtraction.GetOutput())
            centerlineGeometry.SetLengthArrayName(self.lengthArrayName)
            centerlineGeometry.SetCurvatureArrayName(self.curvatureArrayName)
            centerlineGeometry.SetTorsionArrayName(self.torsionArrayName)
            centerlineGeometry.SetTortuosityArrayName(self.tortuosityArrayName)
            centerlineGeometry.SetFrenetTangentArrayName(self.frenetTangentArrayName)
            centerlineGeometry.SetFrenetNormalArrayName(self.frenetNormalArrayName)
            centerlineGeometry.SetFrenetBinormalArrayName(self.frenetBinormalArrayName)
            # centerlineGeometry.SetLineSmoothing(0)
            # centerlineGeometry.SetOutputSmoothedLines(0)
            # centerlineGeometry.SetNumberOfSmoothingIterations(100)
            # centerlineGeometry.SetSmoothingFactor(0.1)
            centerlineGeometry.Update()
            return centerlineGeometry.GetOutput()
        else:
            return networkOutput
    
    def startPointIndexFromEndPointsMarkupsNode(self, endPointsMarkupsNode):
        startPointPosition = [0.0, 0.0, 0.0]
        startPointIndex = -1
        minDistanceSquared = float('inf')
        for i in range(endPointsMarkupsNode.GetNumberOfControlPoints()):
            position = [0.0, 0.0, 0.0]
            endPointsMarkupsNode.GetNthControlPointPosition(i, position)
            distanceSquared = vtk.vtkMath.Distance2BetweenPoints(startPointPosition, position)
            if distanceSquared < minDistanceSquared:
                minDistanceSquared = distanceSquared
                startPointIndex = i
        return startPointIndex
    
    def getEndPoints(self, inputNetworkPolyData, startPointPosition):
        try: 
            import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
            import vtkvmtkMiscPython as vtkvmtkMisc
        except ImportError: 
            print('Erreur import des bibliothèques liées à VMTK')
        
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(inputNetworkPolyData)
        cleaner.Update()
        network = cleaner.GetOutput()
        network.BuildCells()
        network.BuildLinks(0)

        networkPoints = network.GetPoints()
        radiusArray = network.GetPointData().GetArray(self.radiusArrayName)

        startPointId = -1
        maxRadius = 0
        minDistance2 = 0

        endpointIds = vtk.vtkIdList()
        for i in range(network.GetNumberOfCells()):
            numberOfCellPoints = network.GetCell(i).GetNumberOfPoints()
            if numberOfCellPoints < 2:
                continue
            for pointIndex in [0, numberOfCellPoints - 1]:
                pointId = network.GetCell(i).GetPointId(pointIndex)
                pointCells = vtk.vtkIdList()
                network.GetPointCells(pointId, pointCells)
                if pointCells.GetNumberOfIds() == 1:
                    endpointIds.InsertUniqueId(pointId)
                    if startPointPosition is not None: 
                        position = networkPoints.GetPoint(pointId)
                        distance2 = vtk.vtkMath.Distance2BetweenPoints(position, startPointPosition)
                        if startPointId < 0 or distance2 < minDistance2:
                            minDistance2 = distance2
                            startPointId = pointId
                    else: 
                        radius = radiusArray.GetValue(pointId)
                        if startPointId < 0 or radius > maxRadius:
                            maxRadius = radius
                            startPointId = pointId

        endpointPositions = []
        numberOfEndpointIds = endpointIds.GetNumberOfIds()
        print(f"Number of endpoint IDs detected: {numberOfEndpointIds}")
        if numberOfEndpointIds == 0:
            return endpointPositions
        endpointPositions.append(networkPoints.GetPoint(startPointId))
        for pointIdIndex in range(numberOfEndpointIds):
            pointId = endpointIds.GetId(pointIdIndex)
            if pointId == startPointId:
                continue
            endpointPositions.append(networkPoints.GetPoint(pointId))
        print(f"Number of endpoint positions detected: {len(endpointPositions)}")
        return endpointPositions 
    
    def openSurfaceAtPoint(self, polyData, holePosition=None, holePointIndex=None):
        '''
        Modifies the polyData by cutting a hole at the given position.
        '''
        if holePointIndex is None:
            pointLocator = vtk.vtkPointLocator()
            pointLocator.SetDataSet(polyData)
            pointLocator.BuildLocator()
            # find the closest point to the desired hole position
            holePointIndex = pointLocator.FindClosestPoint(holePosition)

        if holePointIndex < 0:
            # Calling GetPoint(-1) would crash the application
            raise ValueError(_("openSurfaceAtPoint failed: empty input polydata"))

        # Tell the polydata to build 'upward' links from points to cells
        polyData.BuildLinks()
        # Mark cells as deleted
        cellIds = vtk.vtkIdList()
        polyData.GetPointCells(holePointIndex, cellIds)
        removeFirstCell = True
        if removeFirstCell:
            # remove first cell only (smaller hole)
            if cellIds.GetNumberOfIds() > 0:
                polyData.DeleteCell(cellIds.GetId(0))
                polyData.RemoveDeletedCells()
        else:
            # remove all cells
            for cellIdIndex in range(cellIds.GetNumberOfIds()):
                polyData.DeleteCell(cellIds.GetId(cellIdIndex))
            polyData.RemoveDeletedCells()   
    
    # Fonction de test qui permet d'afficher les points fiduciaux
    # TODO : A SUPPRIMER
    def loadPointsFromCSV(self, filePath):
        import csv
        import vtk
        nodes = {
            'Branchement': slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode', 'Nœud Branchement'),
            'Terminaison': slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode', 'Nœud Terminaison')
        }
        colors = {
            'Branchement': [0.0, 1.0, 0.0],  # Vert
            'Terminaison': [1.0, 0.5, 0.0]   # Orange
        }
        for node_type, node in nodes.items():
            displayNode = node.GetDisplayNode()
            displayNode.SetVisibility(True)
            displayNode.SetTextScale(0)
            displayNode.SetColor(colors[node_type]) 
            displayNode.SetSelectedColor(colors[node_type]) 
        with open(filePath, 'r', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter=';')
            headers = reader.fieldnames
            if not {'Label', 'X', 'Y', 'Z', 'Type'}.issubset(headers):
                raise ValueError("En-têtes de colonnes manquantes dans le fichier CSV.")
            
            for row in reader:
                x, y, z = -float(row['X']), -float(row['Y']), float(row['Z'])
                label, point_type = row['Label'], row['Type']
                node = nodes[point_type]
                index = node.AddControlPoint(vtk.vtkVector3d(x, y, z))
                node.SetNthControlPointLabel(index, label)
        slicer.app.processEvents()        

    def process(self,
                inputVolume: vtkMRMLScalarVolumeNode,
                outputVolume: vtkMRMLScalarVolumeNode,
                imageThreshold: float,
                invert: bool = False,
                showResult: bool = True) -> None:
        """
        Run the processing algorithm.
        Can be used without GUI widget.
        :param inputVolume: volume to be thresholded
        :param outputVolume: thresholding result
        :param imageThreshold: values above/below this threshold will be set to 0
        :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
        :param showResult: show output volume in slice viewers
        """

        if not inputVolume or not outputVolume:
            raise ValueError("Input or output volume is invalid")

        import time

        startTime = time.time()
        logging.info("Processing started")

        # Compute the thresholded output volume using the "Threshold Scalar Volume" CLI module
        cliParams = {
            "InputVolume": inputVolume.GetID(),
            "OutputVolume": outputVolume.GetID(),
            "ThresholdValue": imageThreshold,
            "ThresholdType": "Above" if invert else "Below",
        }
        cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True, update_display=showResult)
        # We don't need the CLI module node anymore, remove it to not clutter the scene with it
        slicer.mrmlScene.RemoveNode(cliNode)

        stopTime = time.time()
        logging.info(f"Processing completed in {stopTime-startTime:.2f} seconds")
    
    def addNetworkProperties(self, networkPolyData, networkPropertiesTableNode):
        networkPropertiesTableNode.RemoveAllColumns()

        # Cell index column
        numberOfCells = networkPolyData.GetNumberOfCells()
        cellIndexArray = vtk.vtkIntArray()
        cellIndexArray.SetName("CellId")
        cellIndexArray.SetNumberOfValues(numberOfCells)
        for cellIndex in range(numberOfCells):
            cellIndexArray.SetValue(cellIndex, cellIndex)
        networkPropertiesTableNode.GetTable().AddColumn(cellIndexArray)

        # Add length
        lengthArray = networkPolyData.GetCellData().GetArray(self.lengthArrayName)
        if not lengthArray:
            raise ValueError(_("Network polydata does not contain length cell array"))
        networkPropertiesTableNode.GetTable().AddColumn(lengthArray)

        # Add average radius, curvature, torsion values
        for columnName in [self.radiusArrayName, self.curvatureArrayName, self.torsionArrayName]:
            pointDataToCellData = vtk.vtkPointDataToCellData()
            pointDataToCellData.SetInputData(networkPolyData)
            pointDataToCellData.ProcessAllArraysOff()
            pointDataToCellData.AddPointDataArray(columnName)
            pointDataToCellData.Update()
            averageArray = pointDataToCellData.GetOutput().GetCellData().GetArray(columnName)
            if not averageArray:
                raise ValueError(_("Failed to compute array ") + columnName)
            networkPropertiesTableNode.GetTable().AddColumn(averageArray)

        # Add tortuosity
        tortuosityArray = networkPolyData.GetCellData().GetArray(self.tortuosityArrayName)
        if not tortuosityArray:
            raise ValueError(_("Network polydata does not contain length cell array"))
        networkPropertiesTableNode.GetTable().AddColumn(tortuosityArray)

        # Add branch start and end positions
        startPointPositions = vtk.vtkDoubleArray()
        startPointPositions.SetName("StartPointPosition")
        endPointPositions = vtk.vtkDoubleArray()
        endPointPositions.SetName("EndPointPosition")
        for positions in [startPointPositions, endPointPositions]:
            positions.SetNumberOfComponents(3)
            positions.SetComponentName(0, "R")
            positions.SetComponentName(1, "A")
            positions.SetComponentName(2, "S")
            positions.SetNumberOfTuples(numberOfCells)
        for cellIndex in range(numberOfCells):
            pointIds = networkPolyData.GetCell(cellIndex).GetPointIds()
            startPointPosition = [0, 0, 0]
            if pointIds.GetNumberOfIds() > 0:
                networkPolyData.GetPoint(pointIds.GetId(0), startPointPosition)
            if pointIds.GetNumberOfIds() > 1:
                endPointPosition = [0, 0, 0]
                networkPolyData.GetPoint(pointIds.GetId(pointIds.GetNumberOfIds()-1), endPointPosition)
            else:
                endPointPosition = startPointPosition
            startPointPositions.SetTuple3(cellIndex, *startPointPosition)
            endPointPositions.SetTuple3(cellIndex, *endPointPosition)
        networkPropertiesTableNode.GetTable().AddColumn(startPointPositions)
        networkPropertiesTableNode.GetTable().AddColumn(endPointPositions)

        networkPropertiesTableNode.GetTable().Modified()
    
    def addNetworkCurves(self, networkPolyData, centerlineCurveNode, baseName=None):
        shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
        parentItem = shNode.GetItemByDataNode(centerlineCurveNode)

        # remove old children
        shNode.RemoveItemChildren(parentItem)

        if baseName is None:
            baseName = centerlineCurveNode.GetName()

        colorNode = slicer.mrmlScene.GetNodeByID("vtkMRMLColorTableNodeRandom")
        numberOfCells = networkPolyData.GetNumberOfCells()
        slicer.app.pauseRender()
        try:
            radiusArray = networkPolyData.GetPointData().GetArray('Radius')
            for cellId in range(numberOfCells):
                # Create curve node
                curveNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsCurveNode", "{0} ({1})".format(baseName, cellId))
                curveNode.CreateDefaultDisplayNodes()
                color = [0.5, 0.5, 0.5, 1.0]
                colorNode.GetColor(cellId, color)
                curveNode.GetDisplayNode().SetSelectedColor(color[0:3])
                curveNode.SetNumberOfPointsPerInterpolatingSegment(1)
                # Add to subject hierarchy
                curveItem = shNode.GetItemByDataNode(curveNode)
                shNode.SetItemParent(curveItem, parentItem)

                # Add point positions and radius array
                radiusMeasurementArray = vtk.vtkDoubleArray()
                radiusMeasurementArray.SetName('Radius')
                curveNode.SetAttribute("CellId", str(cellId))
                cellPoints = networkPolyData.GetCell(cellId).GetPointIds()
                numberOfCellCurvePoints = cellPoints.GetNumberOfIds()
                for cellPointIdIndex in range(numberOfCellCurvePoints):
                    pointId = cellPoints.GetId(cellPointIdIndex)
                    curveNode.AddControlPointWorld(vtk.vtkVector3d(networkPolyData.GetPoint(pointId)))
                    radiusMeasurementArray.InsertNextValue(radiusArray.GetValue(pointId))

                self._addCurveMeasurementArray(curveNode, radiusMeasurementArray)

                slicer.modules.markups.logic().SetAllControlPointsVisibility(curveNode, False)
        finally:
            slicer.app.resumeRender()
    
    def extractNetworkProperties(self, networkPolyData):
        properties = []
        numberOfCells = networkPolyData.GetNumberOfCells()
        
        centerlineIdsArray = networkPolyData.GetCellData().GetArray(self.centerlineIdsArrayName)
        if not centerlineIdsArray:
            print(_("Network polydata does not contain CenterlineIds cell array"))
        pointDataToCellData = vtk.vtkPointDataToCellData()
        pointDataToCellData.SetInputData(networkPolyData)
        pointDataToCellData.PassPointDataOff()
        pointDataToCellData.Update()
        cellData = pointDataToCellData.GetOutput().GetCellData()

        uniqueCenterlineIds = set()
        for cellIndex in range(numberOfCells):
            centerlineId = centerlineIdsArray.GetValue(cellIndex)
            if centerlineId in uniqueCenterlineIds:
                continue  # Skip if this CenterlineId has already been processed
            uniqueCenterlineIds.add(centerlineId)
            cellProperties = {}
            cellProperties["CenterlineId"] = centerlineId
            cellProperties["CellId"] = cellIndex
            lengthArray = networkPolyData.GetCellData().GetArray(self.lengthArrayName)
            if lengthArray:
                cellProperties["Length in mm"] = f"{lengthArray.GetValue(cellIndex):.2f}"

            radiusArray = cellData.GetArray(self.radiusArrayName)
            if radiusArray:
                cellProperties["Radius"] = f"{radiusArray.GetValue(cellIndex):.2f}"

            curvatureArray = networkPolyData.GetCellData().GetArray(self.curvatureArrayName)
            if curvatureArray:
                cellProperties["Curvature"] = curvatureArray.GetValue(cellIndex)

            torsionArray = networkPolyData.GetCellData().GetArray(self.torsionArrayName)
            if torsionArray:
                cellProperties["Torsion"] = torsionArray.GetValue(cellIndex)

            tortuosityArray = networkPolyData.GetCellData().GetArray(self.tortuosityArrayName)
            if tortuosityArray:
                tortuosity_value = tortuosityArray.GetValue(cellIndex) + 1
                cellProperties["Tortuosity"] = f"{tortuosity_value:.2f}"

            properties.append(cellProperties)
        return properties

    def extractBranches(self, networkPolyData):
        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry
        
        branchExtractor = vtkvmtkComputationalGeometry.vtkvmtkCenterlineBranchExtractor()
        branchExtractor.SetInputData(networkPolyData)
        branchExtractor.SetBlankingArrayName(self.blankingArrayName)
        branchExtractor.SetRadiusArrayName(self.radiusArrayName)
        branchExtractor.SetGroupIdsArrayName(self.groupIdsArrayName)
        branchExtractor.SetCenterlineIdsArrayName(self.centerlineIdsArrayName)
        branchExtractor.SetTractIdsArrayName(self.tractIdsArrayName)
        branchExtractor.Update()

        return branchExtractor.GetOutput()
    
    def extractBifurcationPoints(self, polyData):
        centerlineIds = polyData.GetCellData().GetArray("CenterlineIds")
        bifurcation_points = []
        radii = vtk.vtkFloatArray()  

        for i in range(1, polyData.GetNumberOfCells()):
            currentCellId = centerlineIds.GetTuple1(i)
            previousCellId = centerlineIds.GetTuple1(i - 1)
            if currentCellId != previousCellId:
                cell = polyData.GetCell(i)
                firstPointId = cell.GetPointId(0)
                pointCoords = polyData.GetPoint(firstPointId)
                radius = polyData.GetPointData().GetArray('Radius').GetTuple1(firstPointId)
                bifurcation_points.append(pointCoords)
                radii.InsertNextValue(radius)
        
        return bifurcation_points

    def extractBifurcationPointsReverse(self, polyData):
        centerlineIds = polyData.GetCellData().GetArray("CenterlineIds")
        bifurcation_points = []

        for i in range(1, polyData.GetNumberOfCells()):
            currentCellId = centerlineIds.GetTuple1(i)
            previousCellId = centerlineIds.GetTuple1(i - 1)
            if currentCellId != previousCellId:
                cell = polyData.GetCell(i - 1)
                lastPointId = cell.GetPointId(cell.GetNumberOfPoints() - 1)
                pointCoords = polyData.GetPoint(lastPointId)
                bifurcation_points.append(pointCoords)

        return bifurcation_points

    def mergeBifurcationPoints(self, points1, points2):
        unique_points = set(map(tuple, points1)) | set(map(tuple, points2))
        return unique_points
    
    # Fonction de locator :  
    def find_nearby_centerlines(self, centerline_surface, control_point, search_radius=5):
        locator = vtk.vtkPointLocator()
        locator.SetDataSet(centerline_surface)
        locator.BuildLocator()
        points_in_radius = vtk.vtkIdList()
        locator.FindPointsWithinRadius(search_radius, control_point, points_in_radius)
        centerline_ids = set()
        id_array = centerline_surface.GetCellData().GetArray("CenterlineIds")
        for i in range(points_in_radius.GetNumberOfIds()):
            point_id = points_in_radius.GetId(i)
            cell_ids = vtk.vtkIdList()
            centerline_surface.GetPointCells(point_id, cell_ids)
            for j in range(cell_ids.GetNumberOfIds()):
                cell_id = cell_ids.GetId(j)
                centerline_id = id_array.GetTuple1(cell_id)
                centerline_ids.add(centerline_id)
        return list(centerline_ids)
    
    def classify_bifurcation_points(self, unique_points, centerlinePolyData):
        unique_points = list(unique_points)  # Convertir l'ensemble en liste
        bifurcation_points = vtk.vtkPoints()
        for point in unique_points:
            bifurcation_points.InsertNextPoint(point)

        end_points = []
        branch_points = []
        association_dict = {}

        for i in range(bifurcation_points.GetNumberOfPoints()):
            bifurcation_point = bifurcation_points.GetPoint(i)
            nearby_centerline_ids = self.find_nearby_centerlines(centerlinePolyData, bifurcation_point, search_radius=0.1)

            if len(nearby_centerline_ids) <= 1:  # Si le point est connecté à zéro ou une seule centerline ID
                point_info = {
                    "point": bifurcation_point,
                    "centerline_ids": nearby_centerline_ids
                }
                end_points.append(point_info)
                for centerline_id in nearby_centerline_ids:
                    if centerline_id not in association_dict:
                        association_dict[centerline_id] = {
                            "end_points": [],
                            "branch_points": []
                        }
                    association_dict[centerline_id]["end_points"].append(bifurcation_point)
            else:  # Si le point est connecté à plusieurs centerline IDs
                point_info = {
                    "point": bifurcation_point,
                    "centerline_ids": nearby_centerline_ids
                }
                branch_points.append(point_info)
                for centerline_id in nearby_centerline_ids:
                    if centerline_id not in association_dict:
                        association_dict[centerline_id] = {
                            "end_points": [],
                            "branch_points": []
                        }
                    association_dict[centerline_id]["branch_points"].append(bifurcation_point)

        return end_points, branch_points, association_dict
    
    def visualizeBifurcationPoints(self, points, nodeName="BifurcationPoints"):
        # Create a new fiducial node
        bifurcationFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", nodeName)
    
        # Add points to the fiducial node
        for point in points:
            bifurcationFiducialNode.AddControlPoint(vtk.vtkVector3d(point))
        displayNode = bifurcationFiducialNode.GetDisplayNode()
        if nodeName == "EndPoints":
            displayNode.SetSelectedColor(1, 0, 0)
        else:
            displayNode.SetSelectedColor(0, 1, 0)
        displayNode.SetGlyphScale(1.0)  # Adjust the size as needed
        displayNode.SetTextScale(0.0)  # Hide the text labels if not needed
        displayNode.SetVisibility(True)
        displayNode.SetVisibility3D(True)
        displayNode.SetVisibility2D(True)
    
    def detectOmbilicalPoints(self, centerlinePolyData, unique_bifurcation_points, centerline_associations):
        from vtk.util.numpy_support import vtk_to_numpy
        import numpy as np

        centerline_ids_array = centerlinePolyData.GetCellData().GetArray('CenterlineIds')
        radius_array = centerlinePolyData.GetPointData().GetArray('Radius')
        if centerline_ids_array is None:
            raise ValueError("No centerline Ids were found in the polydata")
        if radius_array is None:
            raise ValueError("No Radius array was found")
    
        centerline_ids = vtk_to_numpy(centerline_ids_array)
        z_min_per_id = {}
        points_per_id = {}
        radius_per_id = {}

        for i in range(centerlinePolyData.GetNumberOfCells()):
            cell = centerlinePolyData.GetCell(i)
            points = cell.GetPoints()
            point_ids = [cell.GetPointId(j) for j in range(points.GetNumberOfPoints())]
            z_values = [points.GetPoint(j)[2] for j in range(points.GetNumberOfPoints())]
            radii = [radius_array.GetTuple1(point_ids[j]) for j in range(points.GetNumberOfPoints())]
            z_min = min(z_values)
            avg_radius = np.mean(radii)
            centerline_id = centerline_ids[i]
            if centerline_id not in z_min_per_id or z_min < z_min_per_id[centerline_id]:
                z_min_per_id[centerline_id] = z_min
                points_per_id[centerline_id] = [points.GetPoint(j) for j in range(points.GetNumberOfPoints())]
                radius_per_id[centerline_id] = avg_radius

        sorted_ids = sorted(z_min_per_id.items(), key=lambda x: x[1])
        ombilical_ids = []
        ombilical_points = {}

        for centerline_id, z_min in sorted_ids:
            if len(ombilical_ids) >= 3:
                break
            is_common_bifurcation = False
            for point in points_per_id[centerline_id]:
                nearby_centerline_ids = self.find_nearby_centerlines(centerlinePolyData, point, search_radius=0.1)
                if any(id in ombilical_ids for id in nearby_centerline_ids):
                    is_common_bifurcation = True
                    break
            if not is_common_bifurcation:
                ombilical_ids.append(centerline_id)
                ombilical_points[centerline_id] = points_per_id[centerline_id][0]
                print(f"Centerline ID gardé: {centerline_id}, Point: {ombilical_points[centerline_id]}")

                # Afficher les endpoints associés
                associated_endpoints = centerline_associations.get(centerline_id, {}).get("end_points", [])
                for endpoint in associated_endpoints:
                    print(f"  Endpoint associé: {endpoint}")

        # Calculer la moyenne des rayons pour chaque CenterlineId sélectionné
        average_radius_per_id = {id: radius_per_id[id] for id in ombilical_ids}

        for ombilical_id in ombilical_ids:
            print(f"Centerline ID: {ombilical_id}, Average Radius: {average_radius_per_id[ombilical_id]}")

        return ombilical_points, average_radius_per_id
    
    # Fonction qui permet la création de points pour les vaisseaux ombilicaux : 
    def markOmbilicalPoints(self, ombilical_points, centerline_associations, centerline_radii):
        ombilical_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "OmbilicalPoints")
        ombilical_node.CreateDefaultDisplayNodes()

        # Trouver le centerline ID avec le plus grand rayon
        max_radius_centerline_id = max(centerline_radii, key=centerline_radii.get)
        used_ids = set()

        for ombilical_id in ombilical_points.keys():
            if ombilical_id == max_radius_centerline_id:
                label = "Veine Ombilicale"
            else:
                label = "Artère Ombilicale"
            associated_endpoints = centerline_associations.get(ombilical_id, {}).get("end_points", [])
            if associated_endpoints:
                min_z_endpoint = min(associated_endpoints, key=lambda pt: pt[2])
                ombilical_node.AddFiducial(min_z_endpoint[0], min_z_endpoint[1], min_z_endpoint[2])
            else:
                point = ombilical_points[ombilical_id]
                ombilical_node.AddFiducial(point[0], point[1], point[2])
            used_ids.add(ombilical_id)
            ombilical_node.SetNthControlPointLabel(ombilical_node.GetNumberOfFiducials() - 1, label)

        # Ajoute les artères restantes s'il y en a moins de 3 points affichés
        if len(used_ids) < 3:
            for ombilical_id in ombilical_points.keys():
                if ombilical_id not in used_ids:
                    label = "Artère Ombilicale"
                    associated_endpoints = centerline_associations.get(ombilical_id, {}).get("end_points", [])
                    if associated_endpoints:
                        min_z_endpoint = min(associated_endpoints, key=lambda pt: pt[2])
                        ombilical_node.AddFiducial(min_z_endpoint[0], min_z_endpoint[1], min_z_endpoint[2])
                    else:
                        point = ombilical_points[ombilical_id]
                        ombilical_node.AddFiducial(point[0], point[1], point[2])
                    ombilical_node.SetNthControlPointLabel(ombilical_node.GetNumberOfFiducials() - 1, label)
                    used_ids.add(ombilical_id)
                    if len(used_ids) >= 3:
                        break
    
    def performInitialization(self, image, lowerThreshold, upperThreshold, sourceSeedIds, targetSeedIds, method="collidingfronts"):
        '''
        '''
        # import the vmtk libraries
        try:
            import vtkvmtkSegmentationPython as vtkvmtkSegmentation
        except ImportError:
            logging.error("Unable to import the SlicerVmtk libraries")

        cast = vtk.vtkImageCast()
        cast.SetInputData(image)
        cast.SetOutputScalarTypeToFloat()
        cast.Update()
        image = cast.GetOutput()

        scalarRange = image.GetScalarRange()

        imageDimensions = image.GetDimensions()
        maxImageDimensions = max(imageDimensions)

        threshold = vtk.vtkImageThreshold()
        threshold.SetInputData(image)
        threshold.ThresholdBetween(lowerThreshold, upperThreshold)
        threshold.ReplaceInOff()
        threshold.ReplaceOutOn()
        threshold.SetOutValue(scalarRange[0] - scalarRange[1])
        threshold.Update()

        thresholdedImage = threshold.GetOutput()

        scalarRange = thresholdedImage.GetScalarRange()

        shiftScale = vtk.vtkImageShiftScale()
        shiftScale.SetInputData(thresholdedImage)
        shiftScale.SetShift(-scalarRange[0])
        shiftScale.SetScale(1.0 / (scalarRange[1] - scalarRange[0]))
        shiftScale.Update()

        speedImage = shiftScale.GetOutput()

        if method == "collidingfronts":
            # ignore sidebranches, use colliding fronts
            logging.debug("Using colliding fronts algorithm")
            logging.debug("number of vtk ids: " + str(sourceSeedIds.GetNumberOfIds()))
            logging.debug("SourceSeedIds:")
            logging.debug(sourceSeedIds)
            collidingFronts = vtkvmtkSegmentation.vtkvmtkCollidingFrontsImageFilter()
            collidingFronts.SetInputData(speedImage)
            sourceSeedId1 = vtk.vtkIdList()
            sourceSeedId1.InsertNextId(sourceSeedIds.GetId(0))
            collidingFronts.SetSeeds1(sourceSeedId1)
            sourceSeedId2 = vtk.vtkIdList()
            sourceSeedId2.InsertNextId(sourceSeedIds.GetId(sourceSeedIds.GetNumberOfIds()-1))
            collidingFronts.SetSeeds2(sourceSeedId2)
            collidingFronts.ApplyConnectivityOn()
            collidingFronts.StopOnTargetsOn()
            collidingFronts.Update()

            subtract = vtk.vtkImageMathematics()
            subtract.SetInputData(collidingFronts.GetOutput())
            subtract.SetOperationToAddConstant()
            subtract.SetConstantC(-10 * collidingFronts.GetNegativeEpsilon())
            subtract.Update()
        
        outImageData = vtk.vtkImageData()
        outImageData.DeepCopy(subtract.GetOutput())

        return outImageData
    
    def buildSimpleLabelMap(self, image, inValue, outValue):

        threshold = vtk.vtkImageThreshold()
        threshold.SetInputData(image)
        threshold.ThresholdByLower(0)
        threshold.ReplaceInOn()
        threshold.ReplaceOutOn()
        threshold.SetOutValue(outValue)
        threshold.SetInValue(inValue)
        threshold.Update()

        outVolumeData = vtk.vtkImageData()
        outVolumeData.DeepCopy(threshold.GetOutput())

        return outVolumeData

    def marchingCubes(self, image, ijkToRasMatrix, threshold):
        
        transformIJKtoRAS = vtk.vtkTransform()
        transformIJKtoRAS.SetMatrix(ijkToRasMatrix)
        marchingCubes = vtk.vtkMarchingCubes()
        marchingCubes.SetInputData(image)
        marchingCubes.SetValue(0, threshold)
        marchingCubes.ComputeScalarsOn()
        marchingCubes.ComputeGradientsOn()
        marchingCubes.ComputeNormalsOn()
        marchingCubes.ReleaseDataFlagOn()
        marchingCubes.Update()

        if transformIJKtoRAS.GetMatrix().Determinant() < 0:
            reverser = vtk.vtkReverseSense()
            reverser.SetInputData(marchingCubes.GetOutput())
            reverser.ReverseNormalsOn()
            reverser.ReleaseDataFlagOn()
            reverser.Update()
            correctedOutput = reverser.GetOutput()
        else:
            correctedOutput = marchingCubes.GetOutput()

        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetInputData(correctedOutput)
        transformer.SetTransform(transformIJKtoRAS)
        transformer.ReleaseDataFlagOn()
        transformer.Update()

        normals = vtk.vtkPolyDataNormals()
        normals.ComputePointNormalsOn()
        normals.SetInputData(transformer.GetOutput())
        normals.SetFeatureAngle(60)
        normals.SetSplitting(1)
        normals.ReleaseDataFlagOn()
        normals.Update()

        stripper = vtk.vtkStripper()
        stripper.SetInputData(normals.GetOutput())
        stripper.ReleaseDataFlagOff()
        stripper.Update()
        stripper.GetOutput()

        result = vtk.vtkPolyData()
        result.DeepCopy(stripper.GetOutput())

        return result
        
    def performEvolution(self, originalImage, segmentationImage, numberOfIterations, inflation, curvature, attraction, levelSetsType='geodesic'):
        '''

        '''
        # import the vmtk libraries
        try:
            import vtkvmtkSegmentationPython as vtkvmtkSegmentation
        except ImportError:
            logging.error("Unable to import the SlicerVmtk libraries")

        featureDerivativeSigma = 0.0
        maximumRMSError = 1E-20
        isoSurfaceValue = 0.0

        logging.debug("NumberOfIterations: " + str(numberOfIterations))
        logging.debug("inflation: " + str(inflation))
        logging.debug("curvature: " + str(curvature))
        logging.debug("attraction: " + str(attraction))

        if levelSetsType == 'geodesic':
            logging.debug("using vtkvmtkGeodesicActiveContourLevelSetImageFilter")
            levelSets = vtkvmtkSegmentation.vtkvmtkGeodesicActiveContourLevelSetImageFilter()
            levelSets.SetFeatureImage(self.buildGradientBasedFeatureImage(originalImage))
            levelSets.SetDerivativeSigma(featureDerivativeSigma)
            levelSets.SetAutoGenerateSpeedAdvection(1)
            levelSets.SetPropagationScaling(inflation * (-1))
            levelSets.SetCurvatureScaling(curvature)
            levelSets.SetAdvectionScaling(attraction * (-1))
        elif levelSetsType == 'curves':
            levelSets = vtkvmtkSegmentation.vtkvmtkCurvesLevelSetImageFilter()
            levelSets.SetFeatureImage(self.buildGradientBasedFeatureImage(originalImage))
            levelSets.SetDerivativeSigma(featureDerivativeSigma)
            levelSets.SetAutoGenerateSpeedAdvection(1)
            levelSets.SetPropagationScaling(inflation * (-1))
            levelSets.SetCurvatureScaling(curvature)
            levelSets.SetAdvectionScaling(attraction * (-1))
        elif levelSetsType == 'threshold':
            raise NotImplementedError()
        elif levelSetsType == 'laplacian':
            raise NotImplementedError()
        else:
            raise NameError('Unsupported LevelSetsType')

        levelSets.SetInputData(segmentationImage)
        levelSets.SetNumberOfIterations(numberOfIterations)
        levelSets.SetIsoSurfaceValue(isoSurfaceValue)
        levelSets.SetMaximumRMSError(maximumRMSError)
        levelSets.SetInterpolateSurfaceLocation(1)
        levelSets.SetUseImageSpacing(1)
        levelSets.Update()

        outImageData = vtk.vtkImageData()
        outImageData.DeepCopy(levelSets.GetOutput())

        return outImageData


    def buildGradientBasedFeatureImage(self, imageData):
        '''
        '''
        # import the vmtk libraries
        try:
            import vtkvmtkSegmentationPython as vtkvmtkSegmentation
        except ImportError:
            logging.error("Unable to import the SlicerVmtk libraries")

        derivativeSigma = 0.0
        sigmoidRemapping = 1

        cast = vtk.vtkImageCast()
        cast.SetInputData(imageData)
        cast.SetOutputScalarTypeToFloat()
        cast.Update()

        if (derivativeSigma > 0.0):
            gradientMagnitude = vtkvmtkSegmentation.vtkvmtkGradientMagnitudeRecursiveGaussianImageFilter()
            gradientMagnitude.SetInputData(cast.GetOutput())
            gradientMagnitude.SetSigma(derivativeSigma)
            gradientMagnitude.SetNormalizeAcrossScale(0)
            gradientMagnitude.Update()
        else:
            gradientMagnitude = vtkvmtkSegmentation.vtkvmtkGradientMagnitudeImageFilter()
            gradientMagnitude.SetInputData(cast.GetOutput())
            gradientMagnitude.Update()

        featureImage = None
        if sigmoidRemapping == 1:
            scalarRange = gradientMagnitude.GetOutput().GetPointData().GetScalars().GetRange()
            inputMinimum = scalarRange[0]
            inputMaximum = scalarRange[1]
            alpha = -(inputMaximum - inputMinimum) / 6.0
            beta = (inputMaximum + inputMinimum) / 2.0

            sigmoid = vtkvmtkSegmentation.vtkvmtkSigmoidImageFilter()
            sigmoid.SetInputData(gradientMagnitude.GetOutput())
            sigmoid.SetAlpha(alpha)
            sigmoid.SetBeta(beta)
            sigmoid.SetOutputMinimum(0.0)
            sigmoid.SetOutputMaximum(1.0)
            sigmoid.Update()
            featureImage = sigmoid.GetOutput()
        else:
            boundedReciprocal = vtkvmtkSegmentation.vtkvmtkBoundedReciprocalImageFilter()
            boundedReciprocal.SetInputData(gradientMagnitude.GetOutput())
            boundedReciprocal.Update()
            featureImage = boundedReciprocal.GetOutput()

        outImageData = vtk.vtkImageData()
        outImageData.DeepCopy(featureImage)

        return outImageData

    def branchextractor(self, arteriabranchCenterline):
        import vtkvmtkComputationalGeometryPython as vtkvmtkComputationalGeometry

        branchExtractor = vtkvmtkComputationalGeometry.vtkvmtkCenterlineBranchExtractor()
        branchExtractor.SetInputData(arteriabranchCenterline)
        branchExtractor.SetBlankingArrayName(self.blankingArrayName)
        branchExtractor.SetRadiusArrayName(self.radiusArrayName)
        branchExtractor.SetGroupIdsArrayName(self.groupIdsArrayName)
        branchExtractor.SetCenterlineIdsArrayName(self.centerlineIdsArrayName)
        branchExtractor.SetTractIdsArrayName(self.tractIdsArrayName)
        branchExtractor.Update()
        centerlines = branchExtractor.GetOutput()

        return centerlines
    
    def combineMultipleLabelMapsAndCreateSegmentation(self, labelMapNames, segmentationNodeName):
        import vtk
        
        labelMapNodes = [slicer.util.getNode(name) for name in labelMapNames if slicer.util.getNode(name)]
        if not labelMapNodes or len(labelMapNodes) < 2:
            slicer.util.errorDisplay("Please provide at least two valid labelmaps to combine.")
            return None

        # Initialiser avec le premier labelmap
        combinedImageData = vtk.vtkImageData()
        combinedImageData.DeepCopy(labelMapNodes[0].GetImageData())

        # Additionner les autres labelmaps
        for labelMapNode in labelMapNodes[1:]:
            addFilter = vtk.vtkImageMathematics()
            addFilter.SetInput1Data(combinedImageData)
            addFilter.SetInput2Data(labelMapNode.GetImageData())
            addFilter.SetOperationToAdd()
            addFilter.Update()
            combinedImageData.DeepCopy(addFilter.GetOutput())

        # Créer un nouveau labelmap à partir du résultat de l'addition
        combinedLabelMap = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", "CombinedLabelMap")
        combinedLabelMap.SetAndObserveImageData(combinedImageData)
        combinedLabelMap.SetSpacing(labelMapNodes[0].GetSpacing())
        combinedLabelMap.SetOrigin(labelMapNodes[0].GetOrigin())
        combinedLabelMap.SetIJKToRASDirections(labelMapNodes[0].GetIJKToRASDirections())

        # Créer un nouveau SegmentationNode à partir du labelmap combiné
        segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", segmentationNodeName)
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(combinedLabelMap, segmentationNode)
        
        # Afficher le SegmentationNode
        segmentationNode.CreateDefaultDisplayNodes()
        segmentationNode.GetDisplayNode().SetVisibility3D(True)
        segmentationNode.GetDisplayNode().SetVisibility2DFill(True)
        segmentationNode.GetDisplayNode().SetVisibility2DOutline(True)
        
        return segmentationNode

######################################################################################################################
######################################################################################################################
######################################################################################################################

# BifurcationsAnalysisTest

######################################################################################################################
######################################################################################################################
######################################################################################################################
class BifurcationsAnalysisTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """Do whatever is needed to reset the state - typically a scene clear will be enough."""
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here."""
        self.setUp()
        self.test_BifurcationsAnalysis1()

    def test_BifurcationsAnalysis1(self):
        """Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")

        # Get/create input data

        import SampleData

        registerSampleData()
        inputVolume = SampleData.downloadSample("BifurcationsAnalysis1")
        self.delayDisplay("Loaded test data set")

        inputScalarRange = inputVolume.GetImageData().GetScalarRange()
        self.assertEqual(inputScalarRange[0], 0)
        self.assertEqual(inputScalarRange[1], 695)

        outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
        threshold = 100

        # Test the module logic

        logic = BifurcationsAnalysisLogic()

        # Test algorithm with non-inverted threshold
        logic.process(inputVolume, outputVolume, threshold, True)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], threshold)

        # Test algorithm with inverted threshold
        logic.process(inputVolume, outputVolume, threshold, False)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], inputScalarRange[1])

        self.delayDisplay("Test passed")