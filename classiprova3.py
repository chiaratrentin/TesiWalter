#!/usr/bin/env python
import sys,math,os, vtk, MC_methods
#import clCenterline, clBranch, clSection
import Tkinter, tkFileDialog
from vmtk import vtkvmtk
from vmtk import pypes
from vmtk import vmtkrenderer
from vmtk import vmtkscripts
#from pypes import pypescript
import vtk_methods
import time

# ------ VARIABILI GLOBALI
warnings = []
errori = []
printings= []
colors=[(0, 0, 0), (1, 0, 0), (0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1), (0.5, 0, 0), (0, 0.5, 0)]
#             black,        red,        blue,       green,    yellow,   purple,     cian,       brown,      dark green

# CODICE INSERITO PER DEBUGGING DA PULIRE:
# vessel.__capping, togliere -ofile, inutile
#vessel.__adjacentBranches togliere il print commentato

# DA SISTEMARE:

# vesse.branching, merge di piu centerline 

class vessel():
    def __init__(self, file):                
        
        if not os.path.exists(file):
            errori.append('(vessel,init): file non esistente')
            pass
        else:
            c = os.path.split(file)
            c = c[1]
            
            # Superfici
            self.Surface = None #surface interactive
            self.Centerline = None #centerline [basic, resampled, merged]
            self.Centerlinetmp = [None, None, None]
            self.Branches = [] # lista degli oggetti branch
            self.GroupId = 0  # qui ho gli ID dei vari branch
            self.BifGroupId= [] #vettore dei GroupId delle biforcazioni
            self.GroupIdsInCenterlineIds = [] # vettore di: vettori contenenti i GroupIds che appartengono al CenterlineId (relativo all'indice)
            self.adjacencesMatrix = None #Tabella delle adiacente, righe centerlineids, colonne groupid. 1 se presente, 0 se non presente
            self.tabellaBiforcazioni = [[], []]
            
            # Filename
            self.dir = os.path.dirname(c)+c[:len(c)-4]+'/'
            if os.path.exists(self.dir)==False:
                os.mkdir(self.dir)
            
            self.fn_base = c[:len(c)-4]
            self.fn = self.dir+self.fn_base+'.vtp' # superficie, una volta scritto come .vtp
            self.fn_bc = self.dir+self.fn_base+'_bc.vtp' # surface clippata
            self.fn_bc_cl = self.dir+self.fn_base+'_bc_cl.vtp' # centerline con calcolo di divisione in branch            
            self.fn_splitted = self.dir+self.fn_base+'_splitted.vtp' # surface intera con calcolo per i branch
            
            # ELABORAZIONE
            self.Reader(file)                                           # Leggo la superficie
            self.Centerlinetmp[0] = centerline(self.fn)   # creo l'oggetto centerline
            self.__branching()                                          # divido il vaso in branch
            self.__adjacentBranches()                             # calcolo delle adiacenze
            self.adjacentMidPoint()                                 # calcolo il punto medio delle Boundaries dei vasi adiacenti
     
    def Reader(self, file=None, surf=None): # legge la superficie dal file
        b = os.path.splitext(file)
        if file != None:
            if b[1]=='.stl':
                self.Surface = vtk_methods.ReadSTL(file)
            elif b[1]=='.vtp':
                self.Surface = vtk_methods.ReadPolyData(file)
            else:
                errori.append('(vessel,reader): Estensione non supportata')
                return None
            if self.Surface == None:
                errori.append('(vessel.reader): La superficie non e stata letta corretamente')
                return None
        elif surf != None:
            self.Surface = surf
        else:
            errori.append('(vessel.reader): parametri mancanti')            
            return None                
        # Scrivo il file
        self.fn = vtk_methods.WritePolyData(self.Surface, self.fn)
    
    def __branching(self):
        # WORKFLOW:
        # 1) Estraggo i branch dalla centerline
        # 2) Divido la superficie
        # 3) Clippo la superficie secondo i GroupIds
        # 4) Threshold delle centerline
        # 5) Se la centerline e divisa -> merge
        
        # 1)
#        if os.path.exists(self.fn_bc_cl):
#           self.Centerline_bc = centerline(None, self.fn_bc_cl)
#           printings.append('CL_bc loaded!')
#        else:
        str0='vmtkbranchextractor -ifile '+self.Centerlinetmp[0].fn
        str1=' -radiusarray MaximumInscribedSphereRadius'
        str2=' -ofile '+ self.fn_bc_cl        
        myArguments = str0+str1+str2
        myPype = pypes.PypeRun(myArguments)
        self.Centerline_bc=myPype.GetScriptObject('vmtkbranchextractor','0').Centerlines
    
        # 2)     
        if os.path.exists(self.fn_bc):
            self.Surface_bc = vtk_methods.ReadPolyData(self.fn_bc)
            printings.append('Surface clippata caricata!')
        else:
            str0='vmtkbranchclipper -ifile '+self.fn+' -centerlinesfile '+ self.fn_bc_cl
            str1=' -groupidsarray GroupIds -blankingarray Blanking -radiusarray MaximumInscribedSphereRadius'        
            str2=' -ofile '+self.fn_bc    
            myArguments1=str0+str1#+str2
            myPype1 = pypes.PypeRun(myArguments1)
            self.Surface_bc = myPype1.GetScriptObject('vmtkbranchclipper','0').Surface
                
        self.Surface_bc.GetPointData().SetActiveScalars('GroupIds')
        tmp=self.Surface_bc.GetPointData().GetScalars().GetRange(0)           
        self.GroupId = range(int(tmp[0]), int(tmp[1])+1)      
        tmp=range(int(tmp[0]), int(tmp[1])+1)
        
        for ids in tmp:
               
            self.Branches.append(branch(ids))
            # creo la cartella contenente tutti i file di ogni branch
            dir = self.dir+'branch_'+str(ids)+'/'
            fn = dir+self.fn_base+'_bc_'+str(ids)
            if os.path.exists(dir)==False:
                os.mkdir(dir)
            # 3)
            self.Branches[ids].fn = fn+'.vtp' #surface
            self.Branches[ids].fn_base = fn
            if os.path.exists(self.Branches[ids].fn):
                self.Branches[ids].Surface = vtk_methods.ReadPolyData(self.Branches[ids].fn)
                printings.append('BC_'+str(ids)+' surface loaded!')
            else:
                str0='vmtkbranchclipper -ifile '+self.fn
                str1=' -centerlinesfile '+self.fn_bc_cl
                str2=' -groupidsarray GroupIds -groupids '+str(ids)+' -blankingarray Blanking -radiusarray MaximumInscribedSphereRadius'
                str3=' -ofile '+self.Branches[ids].fn
                myArguments=str0+str1+str2+str3
                myPype=pypes.PypeRun(myArguments)                 
                self.Branches[ids].Surface = myPype.GetScriptObject('vmtkbranchclipper', '0').Surface
            self.Branches[ids].GroupIds = ids
            # 4) Threshold : Blanking -> GroupIds   
            temp1 = self.__threshold(self.__threshold(self.Centerline_bc, 1, 1, 'Blanking'), ids, ids)
            #Check che siamo nel vaso e non nella biforcazione
            if temp1.GetNumberOfPoints()==0:
                #siamo nel vaso
                print 'Siamo nel vaso'
                temp1 = self.__threshold(self.__threshold(self.Centerline_bc, 0, 0, 'Blanking'), ids, ids, 'GroupIds')
                if os.path.exists(fn+'_cl_basic.vtp'):
                    self.Branches[ids].Centerlinetmp[0] = centerline(None, fn+'_cl_basic.vtp')
                else:
                    self.Branches[ids].Centerlinetmp[0] = centerline(None, 
                                                                     vtk_methods.WritePolyData(temp1, fn+'_cl_basic.vtp') )

                try:                
                    if temp1.GetCellData().GetArray('CenterlineIds').GetNumberOfTuples() > 1 and self.Branches[ids].Blanking == False:
                    #5) controllo se nel vaso ci sono piu centerline, se si faccio il merge
#                    if self.__threshold(temp1, 0, 0, 'CenterlineIds').GetNumberOfPoints()!=0 and self.__threshold(temp1, 1, 1, 'CenterlineIds').GetNumberOfPoints()!=0 and self.Branches[ids].Blanking == False:
                        warnings.append( '(vessel.__branching): Branch ID '+str(ids)+' Ci sono 2 centerline, procedo al merge')
                        if os.path.exists(fn+'_cl_merged.vtp'):
                            self.Branches[ids].Centerlinetmp[2] = centerline(None,fn+'_cl_merged.vtp')
                            printings.append('BC_'+str(ids)+': CL merged loaded!')
                        else:
                            str0 = 'vmtkcenterlinemerge -ifile '+self.Branches[ids].Centerlinetmp[0].fn
                            str1 = ' -radiusarray MaximumInscribedSphereRadius -tractidsarray TractIds -groupidsarray GroupIds -blankingarray Blanking -centerlineidsarray CenterlineIds -length 0.0'#+str(length)
                            str3 = ' -ofile '+fn+'_cl_merged.vtp'
                            myArguments = str0+str1
                            myPype = pypes.PypeRun(myArguments)
                            self.Branches[ids].Centerlinetmp[2] = centerline(None, vtk_methods.WritePolyData(myPype.GetScriptObject('vmtkcenterlinemerge', '0').Centerlines, 
                                                                                                             fn+'_cl_merged.vtp'))                    
                        # ricampiono la centerline mergiata con lo stesso numero di punti che aveva prima
                        try:
                            if temp1.GetNumberOfPoints() != self.Branches[ids].Centerlinetmp[2].Centerline.GetNumberOfPoints():
                                if os.path.exists(fn+'_cl_merged.vtp'):
                                    self.Branches[ids].Centerlinetmp[1] = centerline(None,  fn+'_cl_merged.vtp'  )
                                else:
                                    tempCL = self.Branches[ids].Centerlinetmp[2].resampling(self.Branches[ids].Centerlinetmp[2].Centerline.GetCellData().GetArray('Length').GetTuple(0)[0]/int(temp1.GetNumberOfPoints()))
                                    tempCL_fn = vtk_methods.WritePolyData(tempCL, fn+'_cl_merged.vtp')
                                    warnings.append('(vessel.__branching): Branch ID '+str(ids)+' ricampiono la merged!')
                                    self.Branches[ids].Centerlinetmp[1] = centerline(None,  tempCL_fn  )
                        except:   
                            self.Branches[ids].Centerlinetmp[1]=centerline(None, vtk_methods.WritePolyData(myPype.GetScriptObject('vmtkcenterlinemerge', '0').Centerlines, 
                                                                                                            fn+'_cl_merged.vtp'))                    
                            warnings.append('(vessel.__branching): Branch ID '+str(ids)+' non ho ricampionato!!!')
                except:                
                    warnings.append('(vessel.__branching): Branch ID '+str(ids)+' threshold')
                    self.Branches[ids].Centerlinetmp[1] = centerline(None, vtk_methods.WritePolyData(temp1,
                                                                                                      fn+'_cl_resampled.vtp'))
            else:
                print 'Siamo nella biforcazione'
                self.Branches[ids].Blanking = True #siamo nella biforcazione
                if os.path.exists(fn+'_cl_basic.vtp'):
                    self.Branches[ids].Centerlinetmp[0] = centerline(None,fn+'_cl_basic.vtp', 1)
                    
                else:
                    self.Branches[ids].Centerlinetmp[0] = centerline(None,vtk_methods.WritePolyData(temp1, fn+'_cl_basic.vtp'), 1)
                    
                self.BifGroupId.append(ids)
            
            self.Branches[ids].Execute()
            
            if self.Branches[ids].Centerlinetmp[0] != None:
                temp = self.Branches[ids].Centerlinetmp[0]
            else:
                temp = self.Branches[ids].Centerlinetmp[1]
            
            
                
            temp2 = []
            for j in range(temp.CenterlineIds.GetNumberOfTuples()):
                temp2.append(int(temp.CenterlineIds.GetTuple(j)[0]))
            for j in temp2:    
                try:
                    self.GroupIdsInCenterlineIds[j].append(ids)
                except:
                    self.GroupIdsInCenterlineIds.append([])
                    self.GroupIdsInCenterlineIds[j].append(ids)         
    
    def adjustBifurcationCenterline(self):
        ac=[]
        #creo i 4 punti della centerline
        points = vtk.vtkPoints()
        points1 = vtk.vtkPoints()
        points.SetNumberOfPoints(3)
        points1.SetNumberOfPoints(2)
        #setto i punti
        #punto del branch[0]
        points.SetPoint(0, self.Branches[0].Centerlinetmp[1].Centerline.GetPoints().GetPoint(self.Branches[0].sectioningID[1]))            
        # Centro di massa dei branch adiacenti
        points.SetPoint(1, self.Branches[1].midPoint)
        
        #punto del branch[2]
        points.SetPoint(2, self.Branches[2].Centerlinetmp[0].Centerline.GetPoints().GetPoint(self.Branches[2].sectioningID[0]))
        #punto del branch[3]
        points1.SetPoint(0, self.Branches[1].midPoint)
        points1.SetPoint(1, self.Branches[3].Centerlinetmp[0].Centerline.GetPoints().GetPoint(self.Branches[3].sectioningID[0]))
        #prendo i punti relativi alla prima sezione di tutti i branch vicini
        #creo una linea con questi punti
        for i in range(0, 3):
            ac.append(vtk_methods.CreateSphere(points.GetPoint(i), 0.15, [0, 1, 0]))
            if i < 3:
                ac.append(vtk_methods.CreateSphere(points1.GetPoint(i), 0.15, [0, 0, 1]))
#            for j in range(0, len(self.Branches[1].Centerlinetmp[0].createActor())):
#                ac.append(self.Branches[1].Centerlinetmp[0].createActor()[j])
                
        Lines = vtk.vtkCellArray()
        Lines.InsertNextCell(3)  
        Lines.InsertCellPoint(0)
        Lines.InsertCellPoint(1)
        Lines.InsertCellPoint(2)
        Lines1 = vtk.vtkCellArray()
        Lines1.InsertNextCell(2)  
        Lines1.InsertCellPoint(0)
        Lines1.InsertCellPoint(1)
        Polygon = vtk.vtkPolyData()
        Polygon.SetPolys(Lines)

        ac.append(vtk_methods.CreateActor(Polygon))
        self.viewActors(ac)
        print points
    
    def prova(self):
        ids = self.Branches[1].Centerlinetmp[0].Centerline.GetCellData().GetArray('CenterlineIds').GetRange(0)
        #controllo che la centerline sia divisa in piu parti
        if ids[0]!=ids[1]:
            temp1 = self.__threshold(self.Branches[1].Centerlinetmp[0].Centerline, ids[0], ids[0], 'CenterlineIds')
            temp2 = self.__threshold(self.Branches[1].Centerlinetmp[0].Centerline, ids[1], ids[1], 'CenterlineIds')
            points = vtk.vtkPoints()
                
            points.SetNumberOfPoints(4)

            points.SetPoint(0, self.Branches[0].Centerlinetmp[0].Centerline.GetPoints().GetPoint(self.Branches[0].Centerlinetmp[0].Centerline.GetNumberOfPoints()-1))
            points.SetPoint(1, self.Branches[1].midPoint)
            points.SetPoint(2, temp1.GetPoints().GetPoint(temp1.GetNumberOfPoints()-1))
            points.SetPoint(3, temp2.GetPoints().GetPoint(temp2.GetNumberOfPoints()-1))
             
            polyLine0 = vtk.vtkPolyLine()
            polyLine0.GetPointIds().SetNumberOfIds(2)
            polyLine0.GetPointIds().SetId(0,0)
            polyLine0.GetPointIds().SetId(1,1)
                
            polyLine1 = vtk.vtkPolyLine()
            polyLine1.GetPointIds().SetNumberOfIds(2)
            polyLine1.GetPointIds().SetId(0,1)
            polyLine1.GetPointIds().SetId(1,2)
            
            polyLine2 = vtk.vtkPolyLine()
            polyLine2.GetPointIds().SetNumberOfIds(2)
            polyLine2.GetPointIds().SetId(0,1)
            polyLine2.GetPointIds().SetId(1,3)
            
            cells0 = vtk.vtkCellArray()
            cells0.InsertNextCell(polyLine0)
            cells0.InsertNextCell(polyLine1)
            cells0.InsertNextCell(polyLine2)
            
            polyData = vtk.vtkPolyData()
            polyData.SetPoints(points)
            polyData.SetLines(cells0)
            
            ac=[]
            
            for i in range(0, 4):
                ac.append(vtk_methods.CreateSphere(points.GetPoint(i), 0.15, [1, 1, 0]))
                for j in range(0, len(self.Branches[1].Centerlinetmp[0].createActor())):
                    ac.append(self.Branches[1].Centerlinetmp[0].createActor()[j])
                        
            ac.append(vtk_methods.CreateActor(polyData))
            print polyData
            self.viewActors(ac)

    def adjacentMidPoint(self):#assegna ad ogni branch il punto medio calcolato su tutte le boundaries dei branch a contatto, solo quelle in uscita!
        ac=[]
        for ids in self.GroupId:
            midPoint=[0, 0, 0] #x,y,z
            neighboursID=[]
            neighboursBoundaries = []
            
            #siamo nel vaso o nella biforcazione? se vaso: controllo adjacentIDs che abbiano solo 1 ID!, se biforcazione: pass
            if self.Branches[ids].Blanking == False and self.Branches[ids].adjacentIDs[1] !=[]: # siamo nel vaso, e il vaso ha uscite                
                neighboursID=self.Branches[self.Branches[ids].adjacentIDs[1][0]].adjacentIDs[1] # allora i vicini sono quelli che escono dalla biforcazione
                divider = len(self.Branches[ids].BoundariesTmp[1])
                for i in range(len(self.Branches[ids].BoundariesTmp[1])): # punti delle boundaries branch di partenza ( che sono quelli della zona di uscita) del vaso
                    midPoint[0]= midPoint[0]+self.Branches[ids].BoundariesTmp[1][i][0]
                    midPoint[1]= midPoint[1]+self.Branches[ids].BoundariesTmp[1][i][1]
                    midPoint[2]= midPoint[2]+self.Branches[ids].BoundariesTmp[1][i][2]
                
                for j in neighboursID:
                    divider = divider + len(self.Branches[j].BoundariesTmp[0])
                    for i in range(len(self.Branches[j].BoundariesTmp[0])): #punti dei branch vicini
                        midPoint[0]=midPoint[0]+ self.Branches[j].BoundariesTmp[0][i][0]
                        midPoint[1]=midPoint[1]+ self.Branches[j].BoundariesTmp[0][i][1]
                        midPoint[2]=midPoint[2]+ self.Branches[j].BoundariesTmp[0][i][2]
                midPoint[0]=midPoint[0]/divider
                midPoint[1]=midPoint[1]/divider
                midPoint[2]=midPoint[2]/divider
                self.Branches[self.Branches[ids].adjacentIDs[1][0]].midPoint = midPoint
        
    def __adjacentBranches(self): 
        """assegna ai branch gli id dei branch con i quali confinano seguendo le centerline
        
        """
        for k in range(len(self.Branches)):
            for i in range(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetNumberOfTuples()):
                idx = self.GroupIdsInCenterlineIds[int(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetTuple(i)[0])].index(self.Branches[k].GroupId)
                if self.GroupIdsInCenterlineIds[int(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetTuple(i)[0])][idx] == self.Branches[k].GroupId:                        
                    try:
                        if (self.GroupIdsInCenterlineIds[int(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetTuple(i)[0])][idx+1] not in self.Branches[k].adjacentIDs[1]):
                            self.Branches[k].adjacentIDs[1].append(self.GroupIdsInCenterlineIds[int(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetTuple(i)[0])][idx+1])#quelli in "uscita"
                    except:
                        pass
                    try:
                        if (self.GroupIdsInCenterlineIds[int(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetTuple(i)[0])][idx-1] not in self.Branches[k].adjacentIDs[0]) and k>0:
                            self.Branches[k].adjacentIDs[0].append(self.GroupIdsInCenterlineIds[int(self.Branches[k].Centerlinetmp[0].CenterlineIds.GetTuple(i)[0])][idx-1])#quelli in "entrata"
                    except:
                        pass    
                else:
                    pass                                
            self.Branches[k].adjacentIDs[0].sort()#entranti
            self.Branches[k].adjacentIDs[1].sort()#uscenti
            
    def __threshold(self,line,  high=999, low=0, ID='GroupIds'):
        try:
            thresh = vtk.vtkThreshold()
            thresh.ThresholdBetween(low, high)
            thresh.SetInputArrayToProcess(0, 0, 0, 1, ID)        
            thresh.SetInput(line)
            thresh.Update()
            extract = vtk.vtkGeometryFilter()
            extract.SetInput(thresh.GetOutput())
            extract.Update()
            return extract.GetOutput()
        except:
            errori.append('(vessel.__threshold): errore generico')
        
    def __capping(self, surf):        
        self.__temp=surf[0:len(surf)-4]+'_cap.vtp' #["simple","centerpoint","smooth","annular"]
        str0='vmtksurfacecapper -ifile '+surf
        str1=' -method centerpoint -entityidsarray CellEntityIds -interactive 0'
        str2=' -ofile '+self.__temp
        myArguments = str0+str1+str2
        myPype = pypes.PypeRun(myArguments)        
        return  myPype.GetScriptObject('vmtksurfacecapper','0').Surface
        
    def Viewer(self, surf=None):        
        if surf!=None:
            temp=surf
        else:            
            if os.path.exists(self.fn):                
                temp=vtk_methods.ReadPolyData(self.fn)
            else:
                temp=self.Surface
        actor_surface = vtk_methods.CreateActor(temp)
        self.Surface=temp
        
        prop_surface = actor_surface.GetProperty()
        prop_surface.SetColor(1, 1, 1)
#        prop_surface.SetOpacity(0.5)
        
        renderer = vtk.vtkRenderer()       
        renderer.AddActor(actor_surface)
        
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)
        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)
        
        renderWindow.SetInteractor(renderWindowInteractor)
        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()
        renderWindow.Render()
        renderWindowInteractor.Start()
        
    def ViewerBoth(self, surf=None, cl=None):        
        if self.Centerline!=None:
            tempCL=self.Centerlinetmp[1].Centerline
        elif self.Centerlinetmp[0]!= None:
            tempCL = self.Centerlinetmp[0].Centerline
        else:
            tempCL=self.Centerline
        if surf!=None and cl!=None:            
            tempS=surf
            tempCL=cl
        else:
            tempS=self.Surface
        actor_surface = vtk_methods.CreateActor(tempS)
        actor_cl = vtk_methods.CreateActor(tempCL)      
        prop_surface = actor_surface.GetProperty()
        prop_cl = actor_cl.GetProperty()        
        prop_surface.SetColor(1, 1, 1)
        prop_cl.SetColor(1, 1, 1)
        prop_surface.SetOpacity(0.5)        
        renderer = vtk.vtkRenderer()       
        renderer.AddActor(actor_surface)
        renderer.AddActor(actor_cl)        
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()        
        renderWindow.Render()
        renderWindowInteractor.Start()
        
    def viewActors(self, actor=None):
        """ Visualizza una list di attori
        actor = lista di attori"""
        renderer = vtk.vtkRenderer()
        if actor==None:
            for i in range(0, len(self.actor)):
                renderer.AddActor(self.actor[i])
        else:
            for i in range(0, len(actor)):
                renderer.AddActor(actor[i])            
        renderWindow = vtk.vtkRenderWindow()
        renderer.SetBackground(.1, .2, .3)
        renderWindow.AddRenderer(renderer)        
        renderWindow.SetWindowName(self.fn)
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()        
        renderWindow.Render()
        renderWindowInteractor.Start()
            
    def createActor(self,  input = None, ID= None): #crea tutti gli attori da visualizzare per questa classe
        self.actor = []
        if input==None:
            input = raw_input("Vuoi visualizzare la superficie divisa in branch?(y/n) ")
        if input == 'y':
            self.actor.append(vtk_methods.CreateActor(self.Surface_bc))
            if ID == None:
                ID = raw_input("Secondo quali ID vuoi visualizzarli? (0:CenterlineIds; 1:GroupIds, 2: TractIds,3: Blanking): ")
                if ID == '0':
                    ID = 'CenterlineIds'                    
                elif ID == '1':
                    ID = 'GroupIds'                    
                elif ID == '2':
                    ID = 'TractIds'                    
                elif ID == '3':
                    ID = 'Blanking'
                else:
                    print 'WARNING(vessel.createActor): ID inserito non correttamente, visualizzo ID di default'
                    ID = 'GroupIds'
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInput(self.Surface_bc)
                if ID=='GroupIds':
                    self.Surface_bc.GetPointData().SetActiveScalars(ID)
                    mapper.SetScalarRange(self.Surface_bc.GetPointData().GetScalars().GetRange(0))
                else:
                    self.Surface_bc.GetCellData().SetActiveScalars(ID)            
                    mapper.SetScalarModeToUseCellData()
                mapper.SetScalarRange(self.Surface_bc.GetPointData().GetScalars().GetRange(0))
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                ScalarBarActor = vtk.vtkScalarBarActor()
                ScalarBarActor.SetLookupTable(actor.GetMapper().GetLookupTable())
                ScalarBarActor.GetLabelTextProperty().ItalicOff()
                ScalarBarActor.GetLabelTextProperty().BoldOff()
                ScalarBarActor.GetLabelTextProperty().ShadowOff()
                ScalarBarActor.SetLabelFormat('%1.0f')
                ScalarBarActor.SetTitle(ID)
                self.actor.append(actor) #posizione 0
                self.actor.append(ScalarBarActor) #posizione 1
                prop = []
                tmp=self.Surface_bc.GetPointData().GetScalars().GetRange(0)
                for ids in range(int(tmp[0]), int(tmp[1])):
                    prop.append(self.actor[ids].GetProperty())
                    prop[ids].SetOpacity(0.5)                
        elif input == 'n':
            self.actor.append(vtk_methods.CreateActor(self.Surface)) #posizione 0
            prop = self.actor[0].GetProperty()
            prop.SetOpacity(0.5)
        else:
            print 'ERRORE(vessel.createActor): input inserito NON CORRETTO!!!'
            pass
                
        for i in range(0, len(self.Centerlinetmp[1].createActor())):
            self.actor.append(self.Centerlinetmp[1].createActor()[i])        
        return self.actor
    def ViewerSingleBC(self):
        input = raw_input("Scegli gli ID da visualizzare:"+ "".join(str(self.GroupId)))
        input = int(input)
        self.ViewerBoth(self.Branches[input].Surface, self.Branches[input].Centerline)        
    def ViewerBC(self, ID='GroupIds'):#visualizza la superficie colorata secondo ID: CenterlineIds, GroupIds, TractIds, Blanking
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(self.Surface_bc)
        if ID=='GroupIds':
            self.Surface_bc.GetPointData().SetActiveScalars(ID)
            mapper.SetScalarRange(self.Surface_bc.GetPointData().GetScalars().GetRange(0))
        else:
            self.Surface_bc.GetCellData().SetActiveScalars(ID)            
            mapper.SetScalarModeToUseCellData()
        mapper.SetScalarRange(self.Surface_bc.GetPointData().GetScalars().GetRange(0))
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)  
        
        ScalarBarActor = vtk.vtkScalarBarActor()
        ScalarBarActor.SetLookupTable(actor.GetMapper().GetLookupTable())
        ScalarBarActor.GetLabelTextProperty().ItalicOff()
        ScalarBarActor.GetLabelTextProperty().BoldOff()
        ScalarBarActor.GetLabelTextProperty().ShadowOff()
        ScalarBarActor.SetLabelFormat('%1.0f')
        ScalarBarActor.SetTitle('GroupIds')
        
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.AddActor(ScalarBarActor)
        
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()        
        renderWindow.Render()
        renderWindowInteractor.Start()
    def viewBranchSections(self, GrpId = None, sections=[None, None]):
        """ Visualizza le sezioni nel range sections del branch con GrpId 
        """
        
        ac = []
        prop =[]
        ac.append(vtk_methods.CreateActor(self.Branches[GrpId].Surface))
        prop.append(ac[0].GetProperty())
        prop[0].SetOpacity(0.1)
        for i in range(len(self.Branches[GrpId].Sections_resampled)):
            if self.Branches[GrpId].Sections_resampled[i] != None:                    
                if self.Branches[GrpId].Sections_resampled[i].Surface_resampled != None:
                    if self.Branches[GrpId].Sections_resampled[i] != None:
                        ac1=self.Branches[GrpId].Sections_resampled[i].createActor()
                        ac += ac1
                        m, tempCL = self.Branches[GrpId].selectCenterline()
                        for j in range(tempCL.Centerline.GetNumberOfPoints()):
                            ac.append(vtk_methods.CreateSphere(tempCL.Centerline.GetPoints().GetPoint(j), 0.1, [1, 1, 1]))
                        ac.append(vtk_methods.CreateActor(vtk_methods.ReadPolyData(tempCL.fn)))
        self.viewActors(ac)
        
    def viewAllSections(self, input=None, input1=None):
        """ Visualizza tutte le sezioni della classe Vaso
        """
        ac = []
        prop =[]
        
        ac.append(vtk_methods.CreateActor(self.Surface))
        prop.append(ac[0].GetProperty())
        prop[0].SetOpacity(0.1)
        for i in range(len(self.Branches)):
            if self.Branches[i].Blanking == False:
                prop.append(ac[i].GetProperty())
#                self.Branches[i].sectionizeResampling(input, input1, 1)
                for j in range(len(self.Branches[i].Sections_resampled)):
                    if self.Branches[i].Sections_resampled[j] != None:
#                        self.Branches[i].Sections_resampled[j].resample(input1)
                        prop.append(ac[i].GetProperty())                       
                        if self.Branches[i].Sections_resampled[j].Surface_resampled != None:
                            
                            
                            print 'BC:'+str(i)+' Sections_resampled:'+str(j), self.Branches[i].Sections_resampled[j]
                            if self.Branches[i].Sections_resampled[j] != None:
                                ac1=self.Branches[i].Sections_resampled[j].createActor()
                                ac += ac1
                                
#                                for k in range(self.Branches[i].Sections_resampled[j].Surface_resampled.GetNumberOfPoints()):
#                                if k==0:
#                                    ac.append(vtk_methods.CreateSphere(self.Branches[i].Sections_resampled[j].Surface_resampled.GetPoint(k), 0.2, [1, 1, 0])) #giallo
#                                elif k == 2:
#                                    ac.append(vtk_methods.CreateSphere(self.Branches[i].Sections_resampled[j].Surface_resampled.GetPoint(k), 0.1, [1, 0.7, 0]))#arancione
#                                else:
#                                    ac.append(vtk_methods.CreateSphere(self.Branches[i].Sections_resampled[j].Surface_resampled.GetPoint(k), 0.1, [0, 1, 1]))
                            #-------------------------
#                            ac2=vtk_methods.CreateVersor(
#                                                       self.Branches[i].Sections_resampled[j].clPoint,
#                                                       vtk_methods.vettore(self.Branches[i].Sections_resampled[j].clPoint, 
#                                                                                                 self.Branches[i].Sections_resampled[j].Surface_resampled.GetPoint(k)), 
#                                                        [0, 1, 0]
#                                                       )
#                            for m in range(len(ac2)):
#                                ac.append(ac2[m])
                            
                        #COORDINATE
#                        ac1 = vtk_methods.CreateVersor(
#                                                       self.Branches[i].Sections_resampled[j].clPoint, 
#                                                      self.Branches[i].Sections_resampled[j].NormalsArray)
#                        for m in range(len(ac1)):
#                            ac.append(ac1[m])
#                            prop.append(ac[i].GetProperty())<----
                
#                
#                if self.Branches[i].Centerlinetmp[2] != None:
#                    tempCL = self.Branches[i].Centerlinetmp[2]
#                    print i,'uso CL ricampionata'
#                elif self.Branches[i].Centerlinetmp[1] != None:
#                    tempCL = self.Branches[i].Centerlinetmp[1]
#                    print i,'uso CL'
#                elif self.Branches[i].Centerlinetmp[0] != None:
#                    tempCL = self.Branches[i].Centerlinetmp[0]
#                    print i,'uso CL_basic'
#                else:
#                    print 'nn ho trovato la centerline!!!'
#                    return None
                m, tempCL = self.Branches[i].selectCenterline()
                for j in range(tempCL.Centerline.GetNumberOfPoints()):
                    ac.append(vtk_methods.CreateSphere(tempCL.Centerline.GetPoints().GetPoint(j), 0.1, [1, 1, 1]))
                    #coords
#                    ac.append(vtk_methods.CreateArrow(tempCL.Centerline, [1, 0, 0]))
                ac.append(vtk_methods.CreateActor(vtk_methods.ReadPolyData(tempCL.fn)))
        self.viewActors(ac)
    def selectCenterline(self):
        """ Seleziona la centerlina piu consona da utilizzare per i calcoli
        """
        if self.Centerlinetmp[2] != None:
            return 2, self.Centerlinetmp[2]
        elif self.Centerlinetmp[1] != None:
            return 1, self.Centerlinetmp[1]
        elif self.Centerlinetmp[0] != None:
            return 0, self.Centerlinetmp[0]
        else:
            return -1, None
    
    def printAllPoints(self):
        """ Scrive in una cartella, con nome <nome_file>_txt dei file .txt contenenti i punti (x,y,z) delle sezioni che compongono i vari branch
        """
        dir = self.dir+'/File_txt/'
        if os.path.exists(dir)==False:
            os.mkdir(dir)        
        
        for i in range(len(self.Branches)):
            if self.Branches[i].Blanking == False:
                dir1 = dir +'branch_'+str(i)+'/'
                if os.path.exists(dir1)==False:
                    os.mkdir(dir1)
                k1=1
                for j in range(len(self.Branches[i].Sections_resampled)):
                    if self.Branches[i].Sections_resampled[j] != None:
                        self.Branches[i].Sections_resampled[j].printPoints(dir1, k1)
                        k1+=1
                        
    def Execute(self, n_sez=None, n_pti_sez=None, resamplingStep=None):
        """ Esegue la pipeline di calcolo
        n_sez: numero delle sezioni, se None calcola in base al resemplingStep
        n_pti_sez : numero dei punti delle sezioni
        resamplingStep: passo di ricampionamento        
        """
        for i in range(len(self.Branches)):
            if self.Branches[i].Blanking == False:
                self.Branches[i].sectionizeResampling(n_sez, n_pti_sez, resamplingStep)
                for j in range(len(self.Branches[i].Sections_resampled)):
                    if self.Branches[i].Sections_resampled[j] != None:
                        self.Branches[i].Sections_resampled[j].resample(n_pti_sez)
                        
        print self.fn_base, '\n ------------ REPORT ------------'
        print '\n ------------ ERRORI ------------ \n'
        for i in range(len(errori)):
            print errori[i]
        print '\n ------------ WARNINGS ------------ \n'
        for i in range(len(warnings)):
            print warnings[i]
        print  '\n ------------ PRINTINGS ------------ \n'
        for i in range(len(printings)):
            print  printings[i]               
            
    def ReorderPoints(self, reorder = [None, None]):
        """ Riordina i punti
        reorder: [0] - metodologia di ricampionamento ( parallel, boundaries) [1] - (solo per boundaries, id del punto p)
        """
        
        if reorder[0] == 'first':
            for i in range(len(self.Branches)):
                if self.Branches[i].Blanking == False:
                    self.Branches[i].ReorderPoint(reorder)
                    
        elif reorder[0] == 'boundaries' and reorder[1] != None:
            for i in self.BifGroupId:
                print i , self.Branches[i].p[reorder[1]]
                # branch IN
                self.Branches[ self.Branches[i].adjacentIDs[0][0] ].ReorderPoint([reorder[0], self.Branches[i].p[reorder[1]]  ], 'in')
                # branch OUT
                for k in range(len(self.Branches[i].adjacentIDs[1])):
                    self.Branches[ self.Branches[i].adjacentIDs[1][k]].ReorderPoint([reorder[0], self.Branches[i].p[reorder[1]]  ], 'out')
            
    def computeBoundariesPoints(self, input=50):
        """ Calcola i punti di appoggio per creare le sezioni, da appendere ai branch, relative alla zona di biforcazione
        """
        
        #parametro di debug
        debug =0
        debug_id = 3
        ids=[] #contiene gli ids delle biforcazioni
        for i in range(len(self.GroupId)):
            if self.Branches[i].Blanking != False:
                ids.append(i)
        self.BifGroupId = ids
        boundaries=[]
        
        for i in ids:
            # prendo gli id dei branch in/out della biforcazione i-esima
            id_bc = []
            boundaries=[]
            id_bc.append(self.Branches[i].adjacentIDs[0][0])
            
            for k in range(len(self.Branches[i].adjacentIDs[1])):
                id_bc.append(self.Branches[i].adjacentIDs[1][k])
            #prendo il punto O
            O = self.Branches[i].midPoint
            #prendo le 3 boundaries
            #appendo quella in entrata
            boundaries.append(self.Branches[self.Branches[i].adjacentIDs[0][0]].BoundariesTmp[1])
            #appendo quelle in uscita
            for k in self.Branches[i].adjacentIDs[1]:
                boundaries.append( self.Branches[k].BoundariesTmp[0])
            # Calcolo i punti A,B,C
            # A e relativo al branch in entrata
            
            if debug == 1:
                print 'AdjacentIDs: bif:', i, self.Branches[i].adjacentIDs
                print 'A', len(self.Branches[ self.Branches[i].adjacentIDs[0][0] ].Sections_resampled), self.Branches[ self.Branches[i].adjacentIDs[0][0] ].sectioningID[1]-1
                print 'B', len(self.Branches[ self.Branches[i].adjacentIDs[1][0] ].Sections_resampled), self.Branches[ self.Branches[i].adjacentIDs[1][0] ].sectioningID[1]
                print 'C', len(self.Branches[ self.Branches[i].adjacentIDs[1][1] ].Sections_resampled), self.Branches[ self.Branches[i].adjacentIDs[1][1] ].sectioningID[1]
            
            A = self.Branches[ self.Branches[i].adjacentIDs[0][0] ].Sections_resampled[ self.Branches[ self.Branches[i].adjacentIDs[0][0] ].sectioningID[1]-1 ].clPoint
            B = self.Branches[ self.Branches[i].adjacentIDs[1][0] ].Sections_resampled[ self.Branches[ self.Branches[i].adjacentIDs[1][0] ].sectioningID[0] ].clPoint
            C = self.Branches[ self.Branches[i].adjacentIDs[1][1] ].Sections_resampled[ self.Branches[ self.Branches[i].adjacentIDs[1][1] ].sectioningID[0] ].clPoint            
            
            n=[(0, 0, 0) ]
            # CALCOLO I VERSORI n1,2,3
            
            n.append(vtk_methods.versors(vtk_methods.vettore(O, A)))
            n.append(vtk_methods.versors(vtk_methods.vettore(O, B)))
            n.append(vtk_methods.versors(vtk_methods.vettore(O, C)))
            
            # CALCOLO I VERSORI n4,n5,n6,n7,n8
            n.append(vtk_methods.versors(vtk_methods.somma_vettori(n[1], n[2])))
            n.append(vtk_methods.versors(vtk_methods.somma_vettori(n[1], n[3])))
            n.append(vtk_methods.versors(vtk_methods.somma_vettori(n[2], n[3])))
            n.append(vtk_methods.versors(vtk_methods.vector(n[1], n[2])))
            n.append(vtk_methods.versors(vtk_methods.vector(n[2], n[1])))
            
            # CALCOLO I PUNTI p
            
            p=[ [] , [] , [] ] # i: 1,2,3 boundaries ; j: indice dell n-esimo punto 
            
            for k in range(6): #scorro tutti i p

                for m in range(len(boundaries)): # scorro le boundaries

                    p[m].append(0)
                    zero = math.acos(vtk_methods.scalar( n[k+3],  vtk_methods.normalize(vtk_methods.vettore(O, boundaries[m][0]))))
                    p[m][k]=0
                    for j in range(len(boundaries[m])): #scorro i punti
                        a = math.acos(vtk_methods.scalar( n[k+3],  vtk_methods.normalize(vtk_methods.vettore(O, boundaries[m][j]))))
                        if a < zero:
                            p[m][k] = j
                            zero = a
            if debug == 1 and i == debug_id:
                p_test = p
                
            # SELEZIONO I PUNTI p DA TENERE, in p_new ci sono solo gli indici delle boundaries
            p_new = [ [p[0][1], None, p[0][2], p[0][3], p[0][4]], 
                             [p[1][1], p[1][3], None, p[1][5], p[1][4]], 
                             [None , p[2][3],p[2][2], p[2][1], p[2][4]]]
                     
            # IN p MANTENGO I PUNTI (come coordinate) CHE HANNO ANGOLO MINORE
            p=[]
            for j in range(5):
                for k in range(3):
                    if p_new[k][j] != None:
                        a = vtk_methods.normalize( vtk_methods.vettore(O, boundaries[k][p_new[k][j]] ) )
                        for k1 in range(3):
                            if p_new[k1][j] != None:
                                b = vtk_methods.normalize( vtk_methods.vettore(O, boundaries[k1][p_new[k1][j]] ) )
                                if k!=k1:
                                    min = vtk_methods.min( a, b )
                                    if min == a:
                                        min_idx = [k, j]
                                        min_pt = boundaries[k][ p_new[k][j]]
                                    elif min == b:
                                        min_idx = [k1, j]
                                        min_pt = boundaries[k1][ p_new[k1][j] ]
                p.append(min_pt)
           
            #assegno i punti alle Branches di biforcazione
            
#            self.Branches[ self.Branches[i].adjacentIDs[0][0] ].p = p
            self.Branches[i].p = p
#            self.Branches[ self.Branches[i].adjacentIDs[1][0] ].p = p
#            self.Branches[ self.Branches[i].adjacentIDs[1][1] ].p = p
            
            
            if debug ==1 and i == debug_id:
                # VISUALIZZAZIONE, solo per debug ----------------------
                ac=[]
                # punti p0-5
                for i in range(len(p)):
                    ac.append(vtk_methods.CreateSphere(p[i], 0.1, colors[i] ))
                print p
                self.viewActors(ac)
                #punti p_test
                print p_test[0]
                for i in range(6):

                    ac.append(vtk_methods.CreateSphere( boundaries[2][p_test[2][i]], 0.3, colors[i]))
                self.viewActors(ac)
                # punti O,A,B,C
                ac.append(vtk_methods.CreateSphere(O, 0.1, (0, 0, 0) ))
                ac.append(vtk_methods.CreateSphere(A, 0.1, colors[1] ))
                ac.append(vtk_methods.CreateSphere(B, 0.1, colors[2] ))
                ac.append(vtk_methods.CreateSphere(C, 0.1, colors[3] ))
                print O, A, B, C
                self.viewActors(ac)
                for i in range(len(boundaries)):
                    for j in range(len(p[i])):
                        try:
                            ac.append(vtk_methods.CreateSphere(boundaries[i][ p[i][j] ],  0.11,  colors[j] ))
                        except:
                            pass
                            
                ac1=[]
                # punti delle boundaries
                for i in range(len(boundaries)):
                    for j in range(len(boundaries[i])):
                        ac1.append( vtk_methods.CreateSphere(boundaries[i][j], 0.05, (1, 1, 1) ))
                        for m in range(len(ac1)):
                            ac.append(ac1[m])
                
                # versori
                ac1=[]
                for i in range(9):
                    ac1 = vtk_methods.CreateVersor(O, n[i], colors[i])
                    for m in range(len(ac1)):
                        ac.append(ac1[m])
                #surface
                actor_surface = vtk_methods.CreateActor(self.Branches[3].Surface)
                prop_surface = actor_surface.GetProperty()   
                prop_surface.SetColor(1, 1, 1)
                prop_surface.SetOpacity(0.1)
                ac.append(actor_surface)
                self.viewActors(ac) #DEBUGGING
                #------------------------------------------------
            
            # RICAMPIONATURA DELLE BOUNDARIES
            
            # CREO LE SEMIELLISSI, con delle spline che pasano attraverso i 3 punti
            ellipses=[]
            ellipses_pts=[[], [], []]
            ac=[]
            for j in range(3):
                points = vtk.vtkPoints()
                points.SetNumberOfPoints(3)
                Lines = vtk.vtkCellArray()
                Lines.InsertNextCell(3)
                for i in range(3):
                    Lines.InsertCellPoint(i)                
                points.SetPoint(0,p[4])
                points.SetPoint(1,p[j])
                points.SetPoint(2,p[3])
                
                Polygon = vtk.vtkPolyData()
                Polygon.SetPoints(points)
                Polygon.SetLines(Lines)
        
                cleaner = vtk.vtkCleanPolyData()
                cleaner.SetInput(Polygon)
                cleaner.Update()
                splineFilter = vtk.vtkSplineFilter()
                splineFilter.SetInput(cleaner.GetOutput())
                splineFilter.SetSubdivideToSpecified()
                splineFilter.SetNumberOfSubdivisions(input/2)
                splineFilter.Update()
                ellipses.append( splineFilter.GetOutput())
                
                for m in range(ellipses[j].GetNumberOfPoints()):
                    ellipses_pts[j].append(ellipses[j].GetPoints().GetPoint(m))
                    if debug == 1 and i == debug_id:
                        #DEBUGGING
                        if m == 1:
                            ac.append(vtk_methods.CreateSphere(ellipses[j].GetPoints().GetPoint(m), 0.15, [0.3, 0.3, 0.3]))
                        else:
                            ac.append(vtk_methods.CreateSphere(ellipses[j].GetPoints().GetPoint(m), 0.15, colors[j]))
                if debug == 1:
                    ac.append(vtk_methods.CreateActor(ellipses[j])) #DEBUGGING
            if debug == 1 and i == debug_id:
                self.viewActors(ac) #DEBUGGING
            
            boundaries_new=[]
            
            # boundaries del branch 0
            points  = vtk.vtkPoints()
            points.SetNumberOfPoints(input+1)
            Lines = vtk.vtkCellArray()
            Lines.InsertNextCell(input+1)
            n2=ellipses[0].GetNumberOfPoints()
            for m in range(n2):
                Lines.InsertCellPoint(m)
                points.SetPoint(m, ellipses_pts[0][m])
            n1=n2+ellipses[2].GetNumberOfPoints()-1
            max = ellipses[2].GetNumberOfPoints()-1
            for m in range(n2, n1):
                max += -1
                Lines.InsertCellPoint(m)
                points.SetPoint(m, ellipses_pts[2][max])
            Polygon = vtk.vtkPolyData()
            Polygon.SetPoints(points)
            Polygon.SetLines(Lines)
            cleaner = vtk.vtkCleanPolyData()
            
            cleaner.SetInput(Polygon)
            cleaner.Update()
            boundaries_new.append(cleaner.GetOutput())
            #boundaries del branch 2
            points  = vtk.vtkPoints()
            points.SetNumberOfPoints(input+1)
            Lines = vtk.vtkCellArray()
            Lines.InsertNextCell(input+1)
            n2=ellipses[0].GetNumberOfPoints()
            for m in range(n2):
                Lines.InsertCellPoint(m)
                points.SetPoint(m, ellipses_pts[0][m])
            n1=n2+ellipses[1].GetNumberOfPoints()-1
            max = ellipses[1].GetNumberOfPoints()-1
            for m in range(n2, n1):
                max += -1
                Lines.InsertCellPoint(m)
                points.SetPoint(m, ellipses_pts[1][max])
            Polygon = vtk.vtkPolyData()
            Polygon.SetPoints(points)
            Polygon.SetLines(Lines)
            cleaner = vtk.vtkCleanPolyData()
            
            cleaner.SetInput(Polygon)
            cleaner.Update()
            boundaries_new.append(cleaner.GetOutput())
            #boundaries del branch 3
            points  = vtk.vtkPoints()
            points.SetNumberOfPoints(input+1)
            Lines = vtk.vtkCellArray()
            Lines.InsertNextCell(input+1)
            n2=ellipses[1].GetNumberOfPoints()
            for m in range(n2):
                Lines.InsertCellPoint(m)
                points.SetPoint(m, ellipses_pts[1][m])
            n1=n2+ellipses[2].GetNumberOfPoints()-1
            max = ellipses[2].GetNumberOfPoints()-1
            for m in range(n2, n1):
                max += -1
                Lines.InsertCellPoint(m)
                points.SetPoint(m, ellipses_pts[2][max])
            Polygon = vtk.vtkPolyData()
            Polygon.SetPoints(points)
            Polygon.SetLines(Lines)
            cleaner = vtk.vtkCleanPolyData()
            
            cleaner.SetInput(Polygon)
            cleaner.Update()
            boundaries_new.append(cleaner.GetOutput())
                
            # CREO LE SUPERFICI
            for m in range(3):
                surf = vtk_methods.CreateSurface(boundaries_new[m], O)
                if m == 0:
                    fn = self.Branches[id_bc[0]].Sections_resampled[len( self.Branches[id_bc[0]].Sections_resampled)-1 ].fn
                    fn = fn[0:len(fn)-6]+str(self.Branches[id_bc[0]].sectioningID[1])+'.vtp'
                    vtk_methods.WritePolyData(surf, fn)
                    self.Branches[id_bc[0]].Sections_resampled.append( section(surf,fn, O, self.Branches[id_bc[0]].sectioningID[1]+1 ))
                    self.Branches[id_bc[0]].Sections_resampled[ len(self.Branches[id_bc[0]].Sections_resampled)-1 ].Update(surf)
                    self.Branches[id_bc[0]].Sections_resampled[ len(self.Branches[id_bc[0]].Sections_resampled)-1 ].Boundaries_flag = 1
                else:
                    fn = self.Branches[id_bc[m]].Sections_resampled[ self.Branches[id_bc[m]].sectioningID[1]-1 ].fn
                    fn = fn[0:len(fn)-6]+'0.vtp'
                    vtk_methods.WritePolyData(surf, fn)
                    self.Branches[id_bc[m]].Sections_resampled[0] = section (surf, fn,O,  0)
                    self.Branches[id_bc[m]].Sections_resampled[0].Update(surf)
                    self.Branches[id_bc[m]].Sections_resampled[0].Boundaries_flag = 1
            

        
        if debug == 1 and i == debug_id:
            boundaries=[[], [], []]
           # bound 0, formata da ellissi 0,2
            for j in range(len(ellipses_pts[0])):
                boundaries[0].append(ellipses_pts[0][j])
            j=0
            for k in range(len(ellipses_pts[2])):
                boundaries[0].append(ellipses_pts[2][len(ellipses_pts[2])-k-1])
            k=0
           # bound 1, formata da ellissi 0,1
            for j in range(len(ellipses_pts[0])):
                boundaries[1].append(ellipses_pts[0][j])
            j=0
            for k in range(len(ellipses_pts[1])):
                boundaries[1].append(ellipses_pts[1][len(ellipses_pts[1])-k-1])
            k=0
           # bound 2, formata da ellissi 2,1  
            for j in range(len(ellipses_pts[1])):
                boundaries[2].append(ellipses_pts[1][j])
            for k in range(len(ellipses_pts[2])):
                boundaries[2].append(ellipses_pts[2][len(ellipses_pts[2])-k-1])
            k=0     
            self.Viewer(surf)
   #                print surf
                    
            self.viewActors(ac)

    def viewAllBoundaries(self):
        ac=[]
        prop=[]
        for ids in self.GroupId:
            if self.Branches[ids].Blanking == False:
                # PUNTI:
                #-punto midPoint
                if ids == 0:
                    # BOUNDARIES
                    for i in range(len(self.Branches[ids].BoundariesTmp[1])):
                        ac.append(vtk_methods.CreateSphere(self.Branches[ids].BoundariesTmp[1][i], 0.1, [1, 1, 1]))
                    # MIDPOINTBOUNDARIES
                    ac.append(vtk_methods.CreateSphere(self.Branches[ids].midPointBoundaries[1], 0.3, colors[ids]))
                    # SEZIONE
                    ac1 = self.Branches[ids].Sections_resampled[len(self.Branches[ids].Sections_resampled)-1].createActor()
                    for i in range(len(ac1)):
                        ac.append(ac1[i])
                else:
                    # BOUNDARIES
                    for i in range(len(self.Branches[ids].BoundariesTmp[0][1])):
                        ac.append(vtk_methods.CreateSphere(self.Branches[ids].BoundariesTmp[0][i], 0.1, [1, 1, 1]))
                    # MIDPOINTBOUNDARIES
                    ac.append(vtk_methods.CreateSphere(self.Branches[ids].midPointBoundaries[0], 0.3, colors[ids]))
                    # SEZIONE
                    ac1 = self.Branches[ids].Sections_resampled[0].createActor()
                    for i in range(len(ac1)):
                        ac.append(ac1[i])
            else:
                ac.append(vtk_methods.CreateSphere(self.Branches[ids].midPoint, 0.5, colors[ids]))
                    
        self.viewActors(ac)
    def CreateSectionsVTP(self):
        appendData  = vtk.vtkAppendPolyData()
        for i in range(len(self.Branches)):
            if self.Branches[i].Blanking == False:                
                for j in range(len(self.Branches[i].Sections_resampled)):
                    if self.Branches[i].Sections_resampled[j] != None:
                        appendData .AddInputConnection(self.Branches[i].Sections_resampled[j].selectSurface())
                        print self.Branches[i].Sections_resampled[j].selectSurface()
        printings.append(self.dir+self.fn_base+'_all_sections.vtp')
        print appendData.GetOutput()
        vtk_methods.WritePolyData(appendData.GetOutput(), self.dir+self.fn_base+'_all_sections.vtp')
                    
class centerline(vessel):
    def __init__(self, file=None,fn_cl=None, blanking=0): #file: superficie ; fn_cl: centerline
        if file!=None or fn_cl!=None:
            if file != None and os.path.exists(file): #case1: ho il nome del file devo calcolare la centerline 
                self.fn = file[0:len(file)-4]+'_cl.vtp' #file dove ho la centerline
                #prima provo a caricare il file, se esistente
                if os.path.exists(self.fn):
                    self.fn_attributes = self.fn[0:len(self.fn)-4]+'_attributes.vtp'
                    self.Centerline = vtk_methods.ReadPolyData(self.fn)
                    self.fn_resampled = None
                    
                    printings.append('CL,fn: '+self.fn+' caricata!')
                else:
                    self.__ComputeCenterline(file)
            elif fn_cl != None and os.path.exists(fn_cl): #case2: ho il nome del file della centerline
                self.fn=fn_cl
                self.fn_attributes = self.fn[0:len(self.fn)-4]+'_attributes.vtp'
                self.Centerline = vtk_methods.ReadPolyData(self.fn)
                self.fn_resampled = None
        else:
            print 'ERRORE(centerline,init) parametri mancanti'
            return None  
        self.__setParams()
        if blanking == 0:            
            self.__attributes()            
            self.Geometry()   
        
    def Update(self, fn_cl=None):
        try:
            if fn_cl!=None:
                self.fn = fn_cl
                self.Centerline = vtk_methods.ReadPolyData(fn_cl)            
            self.Geometry()
            self.__setParams()
            self.__attributes()
        except:
            errori.append('CL: ', self.fn, ' - Update() error!!!')
        
    def __ComputeCenterline(self, file):
#        str0='vmtkcenterlines -ifile '+file+' -endpoints 1 -seedselector carotidprofiles'
        str0= 'vmtkcenterlines -ifile '+file+' -endpoints 1'
        str1=' -costfunction 1/R^2 -resampling 1 -resamplingstep 2.5 -ofile '+self.fn
        myArguments=str0+str1
        myPype=pypes.PypeRun(myArguments)
        self.Centerline = myPype.GetScriptObject('vmtkcenterlines', '0').Centerlines
        
    def __attributes(self):
        try:
            str0= 'vmtkcenterlineattributes -ifile '+self.fn+' -abscissasarray Abscissas -normalsarray ParallelTransportNormals'
            myArguments = str0
            myPype=pypes.PypeRun(myArguments)
#            print myPype.GetScriptObject('vmtkcenterlineattributes', '0').Centerlines.GetPointData().GetArray('Abscissas')
            self.AbscissasArray = myPype.GetScriptObject('vmtkcenterlineattributes', '0').Centerlines.GetPointData().GetArray('Abscissas')
            self.NormalsArray = myPype.GetScriptObject('vmtkcenterlineattributes', '0').Centerlines.GetPointData().GetArray('ParallelTransportNormals')
        except:
            errori.append('Centerline.__attributes: FAILED!!')
            print 'Centerline.__attributes: FAILED!!'
        
    def __setParams(self):
        self.LenghtArray = self.Centerline.GetPointData().GetArray('Length')
        self.CurvatureArray = self.Centerline.GetPointData().GetArray('Curvature')
        self.TorsionArray = self.Centerline.GetPointData().GetArray('Torsion')
        self.TortuosityArray = self.Centerline.GetPointData().GetArray('Tortuosity')
        self.FrenetNormalArray = self.Centerline.GetPointData().GetArray('FrenetNormal')
        self.FrenetTangentArray = self.Centerline.GetPointData().GetArray('FrenetTangent')
        self.FrenetBinormalArray = self.Centerline.GetPointData().GetArray('FrenetBinormal')
        self.CenterlineIds = self.Centerline.GetCellData().GetArray('CenterlineIds')
        self.GroupIds = self.Centerline.GetCellData().GetArray('GroupIds')
        self.TractIds = self.Centerline.GetCellData().GetArray('TractIds')
        self.Blanking = self.Centerline.GetCellData().GetArray('Blanking')
        
    def Geometry(self):
        centerlineGeometry = vtkvmtk.vtkvmtkCenterlineGeometry()        
        centerlineGeometry.SetInput(self.Centerline)
        centerlineGeometry.SetLengthArrayName('Length')
        centerlineGeometry.SetCurvatureArrayName('Curvature')
        centerlineGeometry.SetTorsionArrayName('Torsion')
        centerlineGeometry.SetTortuosityArrayName('Tortuosity')
        centerlineGeometry.SetFrenetTangentArrayName('FrenetTangent')
        centerlineGeometry.SetFrenetNormalArrayName('FrenetNormal')
        centerlineGeometry.SetFrenetBinormalArrayName('FrenetBinormal')
        centerlineGeometry.SetLineSmoothing(1)
        centerlineGeometry.SetSmoothingFactor(0.1)
        centerlineGeometry.Update()
        self.Centerline = centerlineGeometry.GetOutput()
        vtk_methods.WritePolyData(self.Centerline, self.fn)

    def Viewer(self):
        actor_cl = vtk_methods.CreateActor(self.Centerline)
        prop_cl = actor_cl.GetProperty()
        prop_cl.SetColor(1, 1, 1)
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor_cl)        
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()
        renderWindow.Render()
        renderWindowInteractor.Start()      
        
    def Viewer_ID(self, ID='GroupIds'): #visualizza la centerline colorata secondo ID: CenterlineIds, GroupIds, TractIds, Blanking
        cellCenters = vtk.vtkCellCenters()
        cellCenters.SetInput(self.Centerline)
        cellCenters.Update()        
        cellCenters.GetOutput().GetPointData().SetActiveScalars(ID)
        
        labelsMapper = vtk.vtkLabeledDataMapper()
        labelsMapper.SetInput(cellCenters.GetOutput())
        labelsMapper.SetLabelModeToLabelScalars()        
        
        labelsActor = vtk.vtkActor2D()
        labelsActor.SetMapper(labelsMapper)
        
        centerlineMapper = vtk.vtkPolyDataMapper()
        centerlineMapper.SetInput(self.Centerline)
        
        centerlineMapper.ScalarVisibilityOn()
        centerlineMapper.SetScalarModeToUseCellData()
        self.Centerlinetmp[1].GetCellData().SetActiveScalars(ID)
        centerlineMapper.SetScalarRange(self.Centerlinetmp[1].GetCellData().GetScalars().GetRange(0))
        
        centerlineActor = vtk.vtkActor()
        centerlineActor.SetMapper(centerlineMapper)
        
        ScalarBarActor = vtk.vtkScalarBarActor()
        ScalarBarActor.SetLookupTable(centerlineActor.GetMapper().GetLookupTable())
        ScalarBarActor.GetLabelTextProperty().ItalicOff()
        ScalarBarActor.GetLabelTextProperty().BoldOff()
        ScalarBarActor.GetLabelTextProperty().ShadowOff()
        ##             self.ScalarBarActor.GetLabelTextProperty().SetColor(0.0,0.0,0.0)
        ScalarBarActor.SetLabelFormat('%1.0f')
        ScalarBarActor.SetTitle('GroupIds')

        renderer = vtk.vtkRenderer()
        renderer.AddActor(labelsActor)         
        renderer.AddActor(centerlineActor)
        renderer.AddActor(ScalarBarActor)
        renderer.SetInteractive(1)
        renderer.TwoSidedLightingOn()
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()        
        renderWindow.Render()
        renderWindowInteractor.Start()
            
    def resampling(self, length=0.0):
        self.fn_resampled = self.fn[0:len(self.fn)-4]+'_resampled.vtp'
        str0='vmtkcenterlineresampling -ifile '+self.fn
        str1=' -length '+str(length)+' -ofile '+self.fn_resampled 
        myArguments=str0+str1
        myPype=pypes.PypeRun(myArguments)        
        return myPype.GetScriptObject('vmtkcenterlineresampling', '0').Centerlines
        
    def  createActor(self, centerline=None): 
        self.actor = []        
        if centerline == None:
            centerline = self.Centerline
        self.actor.append(vtk_methods.CreateActor(centerline))
        for i in range(0, centerline.GetNumberOfPoints()):
            actor_sphere = vtk_methods.CreateSphere(centerline.GetPoint(i))
            prop = actor_sphere.GetProperty()
            prop.SetColor(1, 0, 0) #red
            self.actor.append(actor_sphere)        
        return self.actor
    
class branch(vessel):
    def __init__(self, gID):
        # DATI
        self.Surface = None # superficie del branch
        self.Surface_capped = None #superficie con il cap
        self.Surface_adjusted = None #superficie con le boundaries corrette
        self.Centerlinetmp = [None, None, None] #0)base, 1) centerline, 2) merged
        self.GroupId = gID # ID del branch
        self.CenterlineIds = [] # centerlineIDs
        self.Blanking = False #se true r una biforcazione
        self.Sections = [] # lista degli oggetti sezione
        self.Sections_resampled =[] # lista delle sezioni con bordo ricampionato
        self.adjacentIDs=[[], []] #primo vettore: ID in entrata, secondo vettore: ID in uscita
        self.p=[] #p1,2,3,4,5 delle boundaries
        
        self.Boundaries = None #Cell(0): punti iniziali, cell(1): punti finali
        self.BoundariesTmp =[[], []]
        
        self.midPoint = None #punto relativo al centro di massa dei vasi confinanti. IN,OUT
        
        self.midPointBoundaries=[[0, 0, 0], [0, 0, 0]]  #punto relativo al centro di massa delle boundaries
        self.meanRadius=[None, None] #raggio medio boundaries [0]: parte iniziale, [1]: parte finale
        self.meanRadiusBc = 0 #raggio medio del vaso
        
        self.root = 0 # albero: radice
        self.leaf = 0# albero: foglia
        self.mid = 0 # albero: ramo

        self.sectioningID = [0, None] # range degli id per la sezionatura
                
        # FILE
        self.fn = None        # superficie        
        self.fn_cap = None #superficie cappata        
        self.fn_base = None
#        self.sectionize()
    def Execute(self, flag=None):
        if self.Blanking != True:
            self.__findCenterlineIDs()
            self.__computeBoundaries()
            self.boundaryMidPoint()
            self.__checkBoundaries()
            self.__meanRadius()                     
                
    def __findCenterlineIDs(self):
        try:
            tmp= self.Centerlinetmp[0].CenterlineIds.GetNumberOfTuples()
            
        except:
            errori.append('branch.__findCenterlineIDs: FAILED!!!')
            print 'branch.__findCenterlineIDs: FAILED!!!'
            
    def __meanRadius(self):
        """ Calcola il raggio medio delle boundaries
        """
        if self.Blanking != 1:
            
            radius = 0
            for i in range(len(self.BoundariesTmp[0])):
                radius = radius + vtk_methods.norm(self.midPointBoundaries[0], self.BoundariesTmp[0][i])
            self.meanRadius[0] = radius/len(self.BoundariesTmp[0])
            radius=0
            for i in range(len(self.BoundariesTmp[1])):
                radius = radius + vtk_methods.norm(self.midPointBoundaries[1], self.BoundariesTmp[1][i])
            self.meanRadius[1] = radius/len(self.BoundariesTmp[1])
            
        
    def __computeBoundaries(self):
        """  Calcolo i punti delle boundaries
        """
        #trovo i punti dei bordi
        try:
            boundaryExtractor = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
            boundaryExtractor.SetInput(self.Surface)
            boundaryExtractor.Update()
            self.Boundaries = boundaryExtractor.GetOutput()
        except:
            errori.append('(branch.Execute): Branch ID: '+str(self.GroupId)+'_computeBoundaries failed!')
            print '(branch.Execute): Branch ID: '+str(self.GroupId)+'_computeBoundaries failed!'
        #copio in BoundariesTmp
        for i in range(self.Boundaries.GetCell(0).GetNumberOfPoints()):
            self.BoundariesTmp[0].append(self.Boundaries.GetCell(0).GetPoints().GetPoint(i))
        for i in range(self.Boundaries.GetCell(1).GetNumberOfPoints()):
            self.BoundariesTmp[1].append(self.Boundaries.GetCell(1).GetPoints().GetPoint(i))

    def __checkBoundaries(self):
        """ Controlla che le Boundaries siano ordinate correttamente, 0: quelle in entrata; 1: quelle in uscita
        """
        # controllo che siano salvate correttamente
        # calcolo la distanza del loro centro di massa con il primo e ultimo punto della centerline in modo che siano congruenti
        CL = self.selectCenterline()[1]
        first = CL.Centerline.GetPoints().GetPoint(0)
        last = CL.Centerline.GetPoints().GetPoint(CL.Centerline.GetNumberOfPoints()-1)
        OF = vtk_methods.norm(self.midPointBoundaries[0], first)
        OL = vtk_methods.norm(self.midPointBoundaries[0], last)
        if OF > OL:
            a = self.BoundariesTmp[0]
            self.BoundariesTmp[0] = self.BoundariesTmp[1]
            self.BoundariesTmp[1] = a
            self.boundaryMidPoint()
            self.__meanRadius()
        else:
            pass

    def viewBoundaries(self):
        """ Visualizza le Boundaries del Branch
        """
        ac=[]
        for i in range(len(self.BoundariesTmp[0])):
            ac.append(vtk_methods.CreateSphere(self.BoundariesTmp[0][i], 0.1, [1, 0, 0]))
        for i in range(len(self.BoundariesTmp[1])):
            ac.append(vtk_methods.CreateSphere(self.BoundariesTmp[1][i], 0.1, [0, 1, 0]))
        ac.append(vtk_methods.CreateActor(self.Surface))
        ac.append(vtk_methods.CreateSphere(self.midPointBoundaries[1], 0.5, [1, 0, 0]))#rosso
        self.viewActors(ac)

    def boundaryMidPoint(self):
        """ calcola il punto centro di massa dei bordi [0]: in, 1 out
        """ 
        if self.Boundaries.GetCell(0).GetPoints().GetPoint(0) == self.Boundaries.GetCell(1).GetPoints().GetPoint(1):
            errori.append('Le Boundaries sono incongruenti')
            print 'Le Boundaries sono incongruenti'
        try:
            self.midPointBoundaries=[[0, 0, 0], [0, 0, 0]]
            for j in range(0, 2): #boundaries 0,1
                for i in range(len(self.BoundariesTmp[j])): #numero punti
                    for k in range(0, 3): #x,y,z
                        self.midPointBoundaries[j][k] = self.midPointBoundaries[j][k] + self.BoundariesTmp[j][i][k]
                for k in range(0, 3):
                    self.midPointBoundaries[j][k] = self.midPointBoundaries[j][k] / (len(self.BoundariesTmp[j]))
        except:
            errori.append('(branch.Execute): Branch ID: '+str(self.GroupId)+'boundaryMidPoint failed!')
            print '(branch.Execute): Branch ID: '+str(self.GroupId)+'boundaryMidPoint failed!'

    def sectionizeResampling(self, input=None, input2=None, resamplingStep=None):
        """ Sezionatura del Branch con resampling
        
        input   = Numero di punti per la ricampionatura della centerline del branch
        input2 = Numero di punti per la ricampionatura delle superfici
        resamplingStep = passo di campionatura, espresso in termini di lunghezza"""
        k, tempCL = self.selectCenterline()
        if tempCL == -1:
            errori.append('(branch.sectionizeResampling): Branch ID '+str(GroupId)+' Centerline non trovata!!')
            return None
        print input, input2, resamplingStep
        if (input==None and resamplingStep == None):
            input = int(raw_input('Quanti punti usare per la centerline['+str(k)+']? ('+' '+str(tempCL.Centerline.GetNumberOfPoints())+') '))
        elif input != None:
            input = int(input)
        
        if input2==None:
            input2 = int(raw_input('Quanti punti usare per le sezioni? '))
        elif input2 != None:
            input2 = int(input2)
        
        flag = 0
        if input != None and input != tempCL.Centerline.GetNumberOfPoints() and resamplingStep != None:
            tempCL.resampling(tempCL.Centerline.GetCellData().GetArray('Length').GetTuple(0)[0]/input)            
        elif input == None and resamplingStep != None:
            tempCL.resampling(resamplingStep)
        else:
            flag = 1
            
        if tempCL.fn_resampled != None or flag == 1:
            if k==0:
                self.Centerlinetmp[1] = centerline(None,  tempCL.fn_resampled)
                self.Centerlinetmp[1].Update()
                print '(branch.sectionizeResampling): Branch ID '+str(self.GroupId)+' Centerline[1] creata.'
                k=1
            elif k == 1:
                self.Centerlinetmp[1].Update(tempCL.fn_resampled)
                print '(branch.sectionizeResampling): Branch ID '+str(self.GroupId)+' Centerline[1] aggiornata.'
            elif k == 2:
                self.Centerlinetmp[2].Update(tempCL.fn_resampled)
                print '(branch.sectionizeResampling): Branch ID '+str(self.GroupId)+' Centerline[2] aggiornata.'
        else:
            errori.append('Branch.sectionizeResampling: ID'+str(self.GroupId)+' errore nel ricampionamento!')
            
        #creo la cartella contenente tutti i file        
        dir = os.path.dirname(self.fn)+'/'
        if os.path.exists(dir+'section_resampled')==False:
            os.mkdir(os.path.join(dir, 'section_resampled/'))
        
        self.check = None
        try:
            self.check = self.__findCenterlineIDforSectioning()        
        except:
            errori.append('(branch.Execute): Branch ID: '+str(self.GroupId)+'__findCenterlineIDforSectioning failed!'+str(self.sectioningID))
  
        if self.check == None:
            self.sectioningID[0] = 1
            self.sectioningID[1] = self.Centerlinetmp[k].Centerline.GetNumberOfPoints()-2
        else:        
            if self.sectioningID[0] == 0:
                self.sectioningID[0] = 1
            if self.sectioningID[1] == 0:
                self.sectioningID[1] = self.Centerlinetmp[k].Centerline.GetNumberOfPoints()-2
                warnings.append('Il Bc: '+str(self.GroupId)+' ha non ha sfruttato i sectioningID')
            
        imax = range(self.sectioningID[0], self.sectioningID[1]+1)
        #creo le superfici e le centerline splittate
        j=0
        for i in range(self.sectioningID[1]):
            if i in imax:
                if (self.Centerlinetmp[k].FrenetTangentArray == None 
                        or self.Centerlinetmp[k].FrenetBinormalArray == None 
                        or self.Centerlinetmp[k].FrenetNormalArray == None
                        or self.Centerlinetmp[k].NormalsArray == None
                        or self.Centerlinetmp[k].AbscissasArray == None):
                        self.Centerlinetmp[k].Update()
                
                
                fn = os.path.split(self.fn)[0]
                if not os.path.exists(fn+'/section_resampled/section_'+str(i)+'.vtp'):
                   
                    self.Sections_resampled.append(section(
                                                            self.__cutWithPlane(
                                                                                         self.Surface,
                                                                                         self.Centerlinetmp[k].Centerline.GetPoint(i),
                                                                                         self.Centerlinetmp[k].FrenetTangentArray.GetTuple3(i)
                                                                                          )
                                                                                         ,fn+'/section_resampled/section_'+str(i)+'.vtp'
                                                                                         , self.Centerlinetmp[k].Centerline.GetPoint(i)
                                                                                         , i
                                                                                    )
                                                                        )
                else:
                    self.Sections_resampled.append(section( vtk_methods.ReadPolyData(fn+'/section_resampled/section_'+str(i)+'.vtp')
                                                           ,fn+'/section_resampled/section_'+str(i)+'.vtp'
                                                           , self.Centerlinetmp[k].Centerline.GetPoint(i)
                                                           , i
                                                           )
                                                           )
                    j+=1
                self.Sections_resampled[i].setCoords(
                                                self.Centerlinetmp[k].FrenetBinormalArray.GetTuple3(i),
                                                self.Centerlinetmp[k].FrenetNormalArray.GetTuple3(i),
                                                self.Centerlinetmp[k].FrenetTangentArray.GetTuple3(i), 
                                                self.Centerlinetmp[k].NormalsArray.GetTuple3(i), 
                                                self.Centerlinetmp[k].AbscissasArray.GetTuple3(i)                                                    
                                                )
            else:
                self.Sections_resampled.append(None)
        self.__BranchMeanRadius()
        if j!=0:
            printings.append('Bc: '+str(self.GroupId)+' ha caricato '+str(j)+'/'+str(len(self.Sections_resampled))+' sezioni')
                                                           
    def sectionBoundaries(self,CL, input, dir):
        # il problema qui e che devo creare un nuovo vtkPolyData dai punti delle Boundaries
        try:
            ID = 0
            sectionBoundaries= []
            for i in range(0, 2):            
                if i == 0:
    #                Points = self.Boundaries.GetCell(0).GetPoints()
                    Points = vtk.vtkPoints()
                    for j in range(len(self.BoundariesTmp[0])):
                        Points.InsertNextPoint(self.BoundariesTmp[0][j])
                    origin = CL.Centerline.GetPoint(0)             
                else:
    #                Points = self.Boundaries.GetCell(1).GetPoints()
                    Points = vtk.vtkPoints()
                    for j in range(len(self.BoundariesTmp[1])):
                        Points.InsertNextPoint(self.BoundariesTmp[1][j])
                    origin = CL.Centerline.GetPoint(self.sectioningID[1])
                
                Lines = vtk.vtkCellArray()
    #            Lines.InsertNextCell(Points.GetNumberOfPoints()+1)  
                Lines.InsertNextCell(Points.GetNumberOfPoints()+1)
    #            for i in range(0, Points.GetNumberOfPoints()):
                for i in range(Points.GetNumberOfPoints()):
                    Lines.InsertCellPoint(i)
                Lines.InsertCellPoint(0)
                Polygon = vtk.vtkPolyData()
                Polygon.SetPoints(Points)
                Polygon.SetLines(Lines) 
                
                connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
                connectivityFilter.SetInput(Polygon)
                connectivityFilter.SetClosestPoint(origin)
                connectivityFilter.SetExtractionModeToClosestPointRegion()
                connectivityFilter.Update()        
                contour = connectivityFilter.GetOutput()
                numberOfContourPoints = contour.GetNumberOfPoints()
                section1 = vtk.vtkPolyData()
                sectionPoints = vtk.vtkPoints()
                sectionCellArray = vtk.vtkCellArray() 
                sectionCellArray.InsertNextCell(numberOfContourPoints)        
                for j in range(numberOfContourPoints):
                    point = contour.GetPoint(j)
                    sectionPoints.InsertNextPoint(point)
                    sectionCellArray.InsertCellPoint(j)
                sectionCellArray.InsertCellPoint(0)            
                section1.SetPoints(sectionPoints)
                section1.SetPolys(sectionCellArray)
                sectionBoundaries.append(section1)
                
            return sectionBoundaries      
        except:
            errori.append('branch.sectionBoundaries_'+str(self.GroupId)+' :FAILED!!!')
            print 'branch.sectionBoundaries_'+str(self.GroupId)+' :FAILED!!!'
    
    def viewSections(self, a=None):
        renderer = vtk.vtkRenderer()        
        if a!=None:
            for i in range(0, self.Centerlinetmp[2].Centerline.GetNumberOfPoints()):
                #sezione
                mapper = vtk.vtkPolyDataMapper()
                temp = vtk_methods.ReadPolyData(self.Sections_resampled[i].fn)
                mapper.SetInput(temp)
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)   
                prop=actor.GetProperty()
                prop.SetColor(1, 1, 1)
                renderer.AddActor(actor)
                #sfere centerline
                sphere = vtk.vtkSphereSource()
                actor_sphere=vtk.vtkActor()
                mapper_sphere=vtk.vtkPolyDataMapper()
                sphere.SetCenter(self.Centerlinetmp[2].Centerline.GetPoint(i))
                sphere.SetRadius(0.1) 
                sphere.SetThetaResolution(18)
                sphere.SetPhiResolution(18)                
                mapper_sphere.SetInput(sphere.GetOutput())                                
                prop_shere = actor_sphere.GetProperty()
                prop_shere.SetColor(1.0, 0, 0)
                actor_sphere.SetMapper(mapper_sphere)
                renderer.AddActor(actor_sphere)
                for j in range(0, temp.GetNumberOfPoints()): #sfere sul bordo della superficie
                    sphere = vtk.vtkSphereSource()
                    actor_sphere=vtk.vtkActor()
                    mapper_sphere=vtk.vtkPolyDataMapper()
                    sphere.SetCenter(self.Sections_resampled[i].Surface.GetPoint(j))
                    sphere.SetRadius(0.05) 
                    sphere.SetThetaResolution(18)
                    sphere.SetPhiResolution(18)
                    mapper_sphere.SetInput(sphere.GetOutput())
                    prop_shere = actor_sphere.GetProperty()
                    prop_shere.SetColor(1.0, 1.0, 0)
                    actor_sphere.SetMapper(mapper_sphere)
                    renderer.AddActor(actor_sphere)
            actor_cl = vtk_methods.CreateActor(self.Centerlinetmp[2].Centerline)            
        else:            
            for k in range(0, self.Centerlinetmp[1].Centerline.GetNumberOfPoints()):
                mapper = vtk.vtkPolyDataMapper()
                temp = vtk_methods.ReadPolyData(self.Sections[k].fn)
                mapper.SetInput(temp)
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)  
                prop = actor.GetProperty()
                prop.SetColor(1, 1, 1)
                renderer.AddActor(actor)
            actor_cl = vtk_methods.CreateActor(self.Centerlinetmp[1].Centerline)
        prop_cl = actor_cl.GetProperty()
        prop_cl.SetColor(1, 1, 1)
        renderer.AddActor(actor_cl)        
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()
        renderWindow.Render()
        renderWindowInteractor.Start()

    def __cutWithPlane(self, surface,origin,normal):    
        plane = vtk.vtkPlane()
        plane.SetOrigin(origin)
        plane.SetNormal(normal)        
        cutter = vtk.vtkCutter()
        cutter.SetInput(surface)
        cutter.SetCutFunction(plane)
        cutter.Update()        
        stripper = vtk.vtkStripper()
        stripper.SetInput(cutter.GetOutput())
        stripper.Update()        
        connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
        connectivityFilter.SetInput(stripper.GetOutput())
        connectivityFilter.SetClosestPoint(origin)
        connectivityFilter.SetExtractionModeToClosestPointRegion()
        connectivityFilter.Update()        
        contour = connectivityFilter.GetOutput()
        numberOfContourPoints = contour.GetNumberOfPoints()
        section = vtk.vtkPolyData()
        sectionPoints = vtk.vtkPoints()
        sectionCellArray = vtk.vtkCellArray() 
        sectionCellArray.InsertNextCell(numberOfContourPoints)        
        for j in range(numberOfContourPoints):
            point = contour.GetPoint(j)
            sectionPoints.InsertNextPoint(point)
            sectionCellArray.InsertCellPoint(j)
        sectionCellArray.InsertCellPoint(0)
        section.SetPoints(sectionPoints)
        section.SetPolys(sectionCellArray)
        return section
    
    def createActor(self):
        self.actor=[]
        self.actor.append(vtk_methods.CreateActor(self.Surface))
        prop = self.actor[0].GetProperty()
        prop.SetOpacity(0.5)
        try:
            for i in range(0, len(self.Centerlinetmp[1].createActor())):
                self.actor.append(self.Centerlinetmp[1].createActor()[i])
                    
        except:
            pass
        return self.actor
        
    def createActorSections(self, surf=None):
        self.actor_sections=[]
        if surf==None:
            surf=self.Sections
        for i in range(0, len(self.Centerlinetmp[1].createActor())):
            self.actor_sections.append(self.Centerlinetmp[1].createActor()[i])
        for i in range(4, 9):#len(self.Sections)):
            for j in range(0, len(surf[i].createActor())):
                self.actor_sections.append(surf[i].createActor()[j])
#                print i, j
        return self.actor_sections
        
    def __findCenterlineIDforSectioning(self):
        id = -1
        id0 = -1               
        k, tempCL = self.selectCenterline()
        
        #controllo che il branch non sia troppo piccolo, in quel caso assegno dei valori fittizzi
        rapporto = tempCL.Centerline.GetCellData().GetArray('Length').GetTuple(0)[0] / self.meanRadius[1]
        if rapporto > 1:
            warnings.append('Il Bc: '+str(self.GroupId)+' ha rapporto <1')
            
            if k == -1:
                self.errori.append('(branch.__findCenterlineIDforSectioning): Branch ID '+str(GroupId)+' centerline mancante')
            imax = tempCL.Centerline.GetNumberOfPoints()-1
            
            #top
            dist_1i = vtk_methods.norm(self.midPointBoundaries[1],  tempCL.Centerline.GetPoints().GetPoint(imax))
            while dist_1i < self.meanRadius[1]/2:
                imax+=-1
                dist_1i = vtk_methods.norm(self.midPointBoundaries[1],  tempCL.Centerline.GetPoints().GetPoint(imax))
            errori.append(str(dist_1i / self.meanRadius[1])+' imax:'+str(imax)+'/'+str(tempCL.Centerline.GetNumberOfPoints()-1))
    
            self.sectioningID[1] = imax+1
    
            #bottom
            i = 0
            dist_0i = vtk_methods.norm(self.midPointBoundaries[0], tempCL.Centerline.GetPoints().GetPoint(i))
            while dist_0i < self.meanRadius[0]/2:
                i+=1
                dist_0i = vtk_methods.norm(self.midPointBoundaries[0], tempCL.Centerline.GetPoints().GetPoint(i))
                
            self.sectioningID[0] = i
            
            if imax < i or i==0 or imax == 0:
                errori.append('Il Bc: '+str(self.GroupId)+' non ha ritornato degli IDforSectioning idonei!!!')
                return None
            else:
                return 1
        else:
            self.sectioningID[0] = 0
            self.sectioningID[1] = tempCL.Centerline.GetNumberOfPoints()
            warnings.append('Il Bc: '+str(self.GroupId)+' ha rapporto <1')
            return 1
            
    def __BranchMeanRadius(self):
        k=0
        for i in range(len(self.Sections_resampled)):
            if self.Sections_resampled[i] != None:
                self.meanRadiusBc  += self.Sections_resampled[i].meanRadius
                k+=1
        try:
            self.meanRadiusBc = self.meanRadiusBc / k
        except:
            errori.append( 'BC: '+str(self.GroupId)+' __BranchMeanRadius ERROR: k:'+str(k)+' meanRadiusBc:'+str(self.meanRadiusBc))

            
    def findFirstBoundaries(self, point = None, last=None):
        """Calcola il first point delle sezioni in base alla distanza dal punto dato
        last 0 se partire dalla sezione[0], 1 se partire dalla sezione finale
        """
        if  point == None or last == None:
            errori.append('BC: '+str(self.GroupId)+' parametri mancanti')
        else:
            pto = point
            if last == None or last == 0:
                ids = [] # contiene l'id del punto di ogni sezione
                for i in range(len(self.Sections_resampled)):            
                        ids.append(0)
#                        print '\n-----'+str(i)+'/'+str(len(self.Sections_resampled))+' '+str(len(ids))+' '+str(pto)
                        if self.Sections_resampled[i]!= None and self.Sections_resampled[i].selectSurface().GetNumberOfPoints() > 0:
                            surf = self.Sections_resampled[i].selectSurface()
#                            print surf
                            try:
                                dist = abs(vtk_methods.norm(surf.GetPoints().GetPoint(surf.GetNumberOfPoints()-1), pto))
                                for j in range(surf.GetNumberOfPoints()-1):
                                    if dist > abs(vtk_methods.norm(surf.GetPoints().GetPoint(j), pto)):
                                        ids[i] = j
                                        dist = abs(vtk_methods.norm(surf.GetPoints().GetPoint(j), pto))
                                pto = surf.GetPoints().GetPoint(ids[i])
                            except:
                                pass
            elif last == 1:
                imax = len(self.Sections_resampled)-1
                ids = [0 for i in range(imax)]
                rng = [imax-i-1 for i in range(imax)]
                k=0
                for i in rng:
                    
#                    print '\n-----'+str(i)+'/'+str(imax-1)+' '+str(len(ids))+' '+str(imax-i-1)+str(pto), self.Sections_resampled[i]
                    if self.Sections_resampled[i]!= None and self.Sections_resampled[i].selectSurface().GetNumberOfPoints() > 0:
                            surf = self.Sections_resampled[i].selectSurface()
                            dist = abs(vtk_methods.norm(surf.GetPoints().GetPoint(surf.GetNumberOfPoints()-1), pto) )                
                            for j in range(surf.GetNumberOfPoints()):
                                if dist > abs(vtk_methods.norm(surf.GetPoints().GetPoint(j), pto)):
                                    dist = abs(vtk_methods.norm(surf.GetPoints().GetPoint(j), pto))
                                    ids[k] = j
#                                print 'j: '+str(j)+' dist:'+str(dist)
                            k+=1
                            pto = surf.GetPoints().GetPoint(ids[imax-i-1])
#                            print ids[k], pto
                    else:
                        ids[k] = -1
                        k+=1
                
            return ids
                        
    def ReorderPoint(self, reorder = [None, None], out= None):
        """ Computa il riordino dei punti delle sezioni, secondo due differenti metodi
        reorder [0] == 
         - first: ordina calcolando il punto numero 0 come quello che ha angolo minore rispetto al ParallelTransportNormals
         - boundaries: ordina calcolando il punto che ha distanza minore rispetto al punto dato in reorder [1]
        """
        if reorder[0] == None and reorder[1] == None:
            errori.append('BC: '+str(self.GroupId)+' ReorderPoint, parametri mancanti')
            print 'BC: '+str(self.GroupId)+' ReorderPoint, parametri mancanti'
            
        elif reorder[0] == 'first':
            for i in range(len(self.Sections_resampled)-1):
                if self.Sections_resampled[i] != None and self.Sections_resampled[i].selectSurface().GetNumberOfPoints() > 0:
                    if self.Sections_resampled[i].Boundaries_flag == None:
                        try:
                            self.Sections_resampled[i].findFirstParallel()
                            print 'ReorderPoint: BC:'+str(self.GroupId)+' i:'+str(i)+'/'+str(len(self.Sections_resampled)-1)+' Done'
                        except:
                            print 'ReorderPoint: BC:'+str(self.GroupId)+' i:'+str(i)+'/'+str(len(self.Sections_resampled)-1)+' ERROR'
                            print self.Sections_resampled[i]
                        
        elif reorder[0] == 'boundaries' and reorder[1] != None:
            if out == 'in':
                ids = self.findFirstBoundaries(reorder[1], 1)
#                print ids
                for i in range(len(self.Sections_resampled)-1):
                    if self.Sections_resampled[i] != None and self.Sections_resampled[i].Boundaries_flag == None and self.Sections_resampled[i].selectSurface().GetNumberOfPoints() > 0:
#                        print str(i)+'/'+ str(len(ids))+' '+str(ids[i])
                        if ids[i] != -1 :
                            self.Sections_resampled[i].reorderFromPointGiven( ids[i])
                        else:
                            self.Sections_resampled[i].reorderFromPointGiven(ids[i-1])
            elif out == 'out':
                ids = self.findFirstBoundaries(reorder[1], 0)
#                print ids
                for i in range(len(self.Sections_resampled)-1):
                    if self.Sections_resampled[i] != None and self.Sections_resampled[i].Boundaries_flag == None and self.Sections_resampled[i].selectSurface().GetNumberOfPoints() > 0:
#                        print i, len(ids)
                        if ids[i] != -1 :
                            self.Sections_resampled[i].reorderFromPointGiven( ids[i])
                
            
    
        
class section(vessel): 
    def __init__(self, surf, fn,clPoint,  ID=None):
        #DATI
        self.Surface = surf        
        self.Surface_resampled=None
        self.ID = ID
        self.Normal = None
        self.first = None
        self.clPoint=clPoint
        self.midPoint = None
        self.actor=[]
        self.FrenetBinormalArray = None
        self.FrenetNormalArray = None
        self.FrenetTangentArray = None
        self.NormalsArray = None
        self.AbscissasArray = None
        self.Boundaries_flag = None
        #FILE
        self.fn = vtk_methods.WritePolyData(surf, fn)
        self.fn_resampled = self.fn#[0:len(self.fn)-4]+'_resampled.vtp'        
        #ELABORAZIONE
        self.Origin_points = self.Surface.GetPoints()
        self.midPoint= self.__centralPoint()     
        self.__meanRadius()
        
    def Update(self, surf):
        """ Aggiorna la superficie ricampionata a surf
        
        surf : nuova superficie"""
        self.Surface_resampled = surf
        self.midPoint= self.__centralPoint()       
        
    def __lineResampling(self, input=None):  
        """ Ricampionamento di una linea con una spline
        
        input: numero di punti da utilizzare per il ricampionamento"""      
        if input == None:
            errori.append("ERRORE(section.__lineResampling): parametro mancante!!!")
            pass
        elif self.Origin_points == None:
            return None
        elif input==self.Origin_points.GetNumberOfPoints():
            errori.append("(section.__lineResampling): si e richiesto di ricampionare senza cambiare il numero di punti, aborted!")
            pass
        elif self.Origin_points.GetNumberOfPoints() > 0:
            Lines = vtk.vtkCellArray()
            Lines.InsertNextCell(self.Origin_points.GetNumberOfPoints()+1)  
            for i in range(0, self.Origin_points.GetNumberOfPoints()):
                Lines.InsertCellPoint(i)
            Lines.InsertCellPoint(0)
            Polygon = vtk.vtkPolyData()
            Polygon.SetPoints(self.Origin_points)
            Polygon.SetLines(Lines)            
            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInput(Polygon)
            cleaner.Update()
            splineFilter = vtk.vtkSplineFilter()
            splineFilter.SetInput(cleaner.GetOutput())
            splineFilter.SetSubdivideToSpecified()
            splineFilter.SetNumberOfSubdivisions(input)            
            splineFilter.Update()
            return splineFilter.GetOutput()
        else:
            return None
                
        
    def resample(self, input=None):
        """Ricampionamento della superficie
        
        input: numero di punti da utilizzare"""
        if input==None:
            input = raw_input("Quanti punti vuoi usare? ("+ " ".join(str(self.Origin_points.GetNumberOfPoints()))+') ')                    
        a=self.__lineResampling(int(input))
        
        if a!=None and self.midPoint!= [None, None, None]:
            connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
            connectivityFilter.SetInput(a)
            connectivityFilter.SetClosestPoint(self.midPoint)
            connectivityFilter.SetExtractionModeToClosestPointRegion()
            connectivityFilter.Update()        
            contour = connectivityFilter.GetOutput()
            numberOfContourPoints = contour.GetNumberOfPoints()
            self.Surface_resampled = vtk.vtkPolyData()
            sectionPoints = vtk.vtkPoints()
            sectionCellArray = vtk.vtkCellArray() 
            sectionCellArray.InsertNextCell(numberOfContourPoints)        
            for j in range(numberOfContourPoints):
                point = contour.GetPoint(j)
                sectionPoints.InsertNextPoint(point)
                sectionCellArray.InsertCellPoint(j)
            sectionCellArray.InsertCellPoint(0)
            self.Surface_resampled.SetPoints(sectionPoints)
            self.Surface_resampled.SetPolys(sectionCellArray)
            vtk_methods.WritePolyData(self.Surface_resampled, self.fn_resampled)
            
        else:
            return None
    
    def __centralPoint(self):
        """ Calcola il centro di massa della sezione
        """
        self.midPoint = [0, 0, 0]
        surf = self.selectSurface()
        for i in range(surf.GetNumberOfPoints()):
            for j in range(3):
                self.midPoint[j]+=surf.GetPoint(i)[j]
        try:
            return [self.midPoint[0]/surf.GetNumberOfPoints(), self.midPoint[1]/surf.GetNumberOfPoints() , self.midPoint[2]/surf.GetNumberOfPoints()]
        except:
            return None
    def __meanRadius(self):
        surf = self.selectSurface()
        self.meanRadius = 0
        if surf.GetNumberOfPoints() > 0:
            for i in range(surf.GetNumberOfPoints()):
                self.meanRadius += vtk_methods.norm(self.midPoint, surf.GetPoint(i))
            self.meanRadius = self.meanRadius/surf.GetNumberOfPoints()
            
            
    def viewSection(self):
        surf = self.selectSurface()
        renderer = vtk.vtkRenderer()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(surf)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        prop=actor.GetProperty()
        prop.SetColor(1, 1, 1)
        renderer.AddActor(actor)
        #sfera midpoint
        sphere = vtk.vtkSphereSource()
        actor_sphere=vtk.vtkActor()
        mapper_sphere=vtk.vtkPolyDataMapper()
        sphere.SetCenter(self.midPoint)
        sphere.SetRadius(0.1) 
        sphere.SetThetaResolution(18)
        sphere.SetPhiResolution(18)
        mapper_sphere.SetInput(sphere.GetOutput())
        prop_shere = actor_sphere.GetProperty()
        prop_shere.SetColor(1.0, 1.0, 0)
        actor_sphere.SetMapper(mapper_sphere)
        renderer.AddActor(actor_sphere)
        for j in range(0, surf.GetNumberOfPoints()): #sfere sul bordo della superficie
            sphere = vtk.vtkSphereSource()
            actor_sphere=vtk.vtkActor()
            mapper_sphere=vtk.vtkPolyDataMapper()
            sphere.SetCenter(surf.GetPoint(j))
            sphere.SetRadius(0.1) 
            sphere.SetThetaResolution(18)
            sphere.SetPhiResolution(18)
            mapper_sphere.SetInput(sphere.GetOutput())
            prop_shere = actor_sphere.GetProperty()
            prop_shere.SetColor(1.0, 0, 0)
            actor_sphere.SetMapper(mapper_sphere)
            renderer.AddActor(actor_sphere)

        ac1=vtk_methods.CreateCoords(self.midPoint, self.FrenetBinormalArray, self.FrenetNormalArray, self.FrenetTangentArray)
        for i in range(len(ac1)):
            renderer.AddActor(ac1[i])
        
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)        
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)        
        renderWindow.SetInteractor(renderWindowInteractor)        
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        renderWindowInteractor.SetInteractorStyle(interactorStyle)
        renderWindowInteractor.Initialize()
        renderWindow.Render()
        renderWindowInteractor.Start()   
    
    def createActor(self, sphere=None):
        self.actor=[]
        #punti del contorno il numero 1 e colorato di giallo, in 2 arancione
        surf = self.selectSurface()
        for i in range(surf.GetNumberOfPoints()):
            if i==1:
                self.actor.append(vtk_methods.CreateSphere(surf.GetPoints().GetPoint(i) , 0.3, [1, 0, 1] ))
            elif i==2:
                self.actor.append(vtk_methods.CreateSphere(surf.GetPoints().GetPoint(i), 0.2, [1, 0.7, 0]))
#            else:
#                self.actor.append(vtk_methods.CreateSphere(self.Surface.GetPoints().GetPoint(i) , 0.1, [0, 1, 1]))
            Points = vtk.vtkPoints()
            Lines = vtk.vtkCellArray()
            Lines.InsertNextCell(surf.GetNumberOfPoints()+1)
            for j in range(surf.GetNumberOfPoints()):
                Points.InsertNextPoint(surf.GetPoints().GetPoint(j))
                Lines.InsertCellPoint(j)
            Lines.InsertCellPoint(0)
            Polygon = vtk.vtkPolyData()
            Polygon.SetPoints(Points)
            Polygon.SetLines(Lines)
            self.actor.append(vtk_methods.CreateActor(Polygon))
        return self.actor
        
    def setCoords(self, fb = None, fn = None, ft = None, na = None, ab = None):
        """ Setta le coordinate della CL, nella sezione
        """
        self.FrenetBinormalArray = fb
        self.FrenetNormalArray = fn
        self.FrenetTangentArray = ft
        self.NormalsArray = na
        self.AbscissasArray = ab
        
    def selectSurface(self):
        """ Seleziona la superficie piu consona
        """
        if self.Surface_resampled != None:
            if self.Surface_resampled.GetNumberOfPoints()>0:
                return self.Surface_resampled
            else:
                return self.Surface
        else:
                return self.Surface
    
    def findFirstParallel(self):
        """ Cerca il punto della sezione piu vicino al vettore ParallelTransportNormals
        """
        zero = math.acos(vtk_methods.scalar(
                                      self.NormalsArray, 
                                      vtk_methods.normalize(vtk_methods.vettore(self.clPoint, self.Surface_resampled.GetPoints().GetPoint(0))
                                      )))
        first = 0
        for i in range(self.Surface_resampled.GetNumberOfPoints()):
            a = math.acos(vtk_methods.scalar(
                                       self.NormalsArray,
                                       vtk_methods.normalize(vtk_methods.vettore(self.clPoint, self.Surface_resampled.GetPoints().GetPoint(i)))
                                       ))
            if a < zero:
                first = i
                zero = a
        self.first = first
        if first != 0:
            self.__reorderPoints()        
        return first
        
    def printPoints(self, dir, k):
        """
        dir = directory"""
        out_file = open(dir+'section_'+str(k)+'.txt',"w")
        surf = self.selectSurface()
        for i in range(surf.GetNumberOfPoints()):
            out_file.write(str(surf.GetPoints().GetPoint(i)[0])+' '+str(surf.GetPoints().GetPoint(i)[1])+' '+str(surf.GetPoints().GetPoint(i)[2])+'\n')
        out_file.close()
    def __reorderPoints(self):
        """ Riordina i punti congruentemente con il punto first
        """
        debug = 1
        
        if self.Boundaries_flag == None:
            Lines = vtk.vtkCellArray()
            surf = self.selectSurface()
            Lines.InsertNextCell(surf.GetNumberOfPoints()+1)  
            Points = vtk.vtkPoints()
            k=0
            if debug == 1:
                print 'Ranges, 1:'+str(range(self.first, surf.GetNumberOfPoints()))+' 2:'+str(range(0, self.first+1)) 
                print 'Surf: '+str(self.ID), surf!= None
            for i in range(self.first, surf.GetNumberOfPoints()):
                Lines.InsertCellPoint(k)
                Points.InsertNextPoint(surf.GetPoints().GetPoint(i))
                k+=1
            for i in range(0, self.first+1):
                Lines.InsertCellPoint(k)
                Points.InsertNextPoint(surf.GetPoints().GetPoint(i))
                k+=1
            Lines.InsertCellPoint(0)
            Polygon = vtk.vtkPolyData()
            Polygon.SetPoints(Points)
            Polygon.SetLines(Lines)
            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInput(Polygon)
            cleaner.Update()
            a = cleaner.GetOutput()
            connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
            connectivityFilter.SetInput(a)
            connectivityFilter.SetClosestPoint(self.midPoint)
            connectivityFilter.SetExtractionModeToClosestPointRegion()
            connectivityFilter.Update()        
            contour = connectivityFilter.GetOutput()
            numberOfContourPoints = contour.GetNumberOfPoints()
            self.Surface_resampled = vtk.vtkPolyData()
            sectionPoints = vtk.vtkPoints()
            sectionCellArray = vtk.vtkCellArray() 
            sectionCellArray.InsertNextCell(numberOfContourPoints)        
            for j in range(numberOfContourPoints):
                point = contour.GetPoint(j)
                sectionPoints.InsertNextPoint(point)
                sectionCellArray.InsertCellPoint(j)
            sectionCellArray.InsertCellPoint(0)
            self.Surface_resampled.SetPoints(sectionPoints)
            self.Surface_resampled.SetPolys(sectionCellArray)
            vtk_methods.WritePolyData(self.Surface_resampled, self.fn_resampled)
        
    def reorderFromPointGiven(self,id=None):
        """ Riordina i punti a partire dal punto passato come parametro, si appoggia al metodo __reorderPoints, modificando self.first
        """
        if id!= None:
            self.first = id
            self.__reorderPoints()
            
        
#root = Tkinter.Tk()
#print root
## Carotid_multimple_branches.vtp , PRE_INNER_WALL.stl
#FileIn = tkFileDialog.askopenfile(parent=root, mode='rb',title='Choose a file', initialfile='')
#root.destroy()
#vaso=vessel(FileIn.name)

vaso = vessel('/home/walter/Desktop/Tesi_23_3_14/tesi/File STL/Carotid_multimple_branches.vtp')
vaso.Execute(None, 10, 2)
print 'Execute Done'
vaso.computeBoundariesPoints()
print 'CopmputeBoundariesPoints Done'
#vaso.viewAllSections()
print vaso.BifGroupId
vaso.viewBranchSections(6)

vaso.ReorderPoints(['first', None])
#vaso.ReorderPoints(['boundaries', 4])
print 'ReorderPoints Done'

#vaso.CreateSectionsVTP()
vaso.viewAllSections()
#vaso.printAllPoints()


#files = ['PRE_INNER_WALL.stl', 'Carotid_multimple_branches.vtp', 'STLexpORIG_noCap.stl', 'Lumen_pre_noCap.stl', 'CA_SX_SMOOTH_bin_noCap.stl']
#files = 'Carotid_multimple_branches.vtp'
#dir_files = '/home/walter/Dropbox/Tesi/'

#vaso=[]
#for i in range(len(files)):
#    print '\n Eseguo:'+dir_files+files[i]
#    vaso.append(vessel(dir_files+files[i]))
#    vaso[i].Execute(20, 50)
#    vaso[i].computeBoundariesPoints()
#    time.sleep(5)
#    vaso[i].printAllPoints()
#    
#    print 'Elaborazione terminata'
#    vaso[i].viewAllSections()
#    errori=[]
#    warnings=[]
#    printings=[]
#    except:
#        errori.append(files[i]+' - Non funzionante')
    
#vaso=vessel('Carotid_multimple_branches.vtp')
#
##vaso = vessel('PRE_INNER_WALL.stl')
#
#print 'Creato'

#vaso.printAllPoints()



#vaso.viewAllBoundaries()







    

#vaso.viewAllSections()


