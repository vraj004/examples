!> \file
!> $Id: MoreComplexMeshExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program which sets up a field which uses a more complex mesh using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example Ca_Dynamics_3DExample.f90
!! Example program which sets up a singpe species diffusion problem.
!! \par Latest Builds:
!<

!> Main program
PROGRAM ELEMVOLSEXAMPLE

  USE OPENCMISS
  USE MPI
  USE FIELDML_API


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: LENGTH=10.0_CMISSDP,WIDTH=10.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=10.0_CMISSDP
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: ElemVolFieldUserNumber=19


  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,CONDITION,NUMBER_GLOBAL_Y_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS,NODE_NUMBER,NodeDomain,NUMBER_GLOBAL_ELEMENTS
  INTEGER(CMISSIntg),DIMENSION(50) :: BCNODES
  INTEGER(CMISSIntg),DIMENSION(9) :: RYRNODES
  INTEGER(CMISSIntg) :: MPI_IERROR,NUMBER_GLOBAL_Z_ELEMENTS,NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: node,STARTNODE,ENDNODE,ne,number_surrounding_elements
  REAL(CMISSDP) :: VALUE,elem_volume
  INTEGER(CMISSIntg) :: constantModelIndex,GeometricMeshComponent
  INTEGER(CMISSIntg), POINTER :: SURROUNDING_ELEMENTS(:)

  !INTEGER(INTG) :: first_global_dof,first_local_dof,first_local_rank,last_global_dof,last_local_dof,last_local_rank,rank_idx
  !INTEGER(INTG) :: EQUATIONS_SET_INDEX
  !TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  
    !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,SourceField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField,ElemVolField
  TYPE(CMISSMeshElementsType) :: MeshElements
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSFieldMLIOType) :: fieldmlInfo,outputInfo
  INTEGER(CMISSIntg) :: typeHandle

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = "."
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = outputDirectory//"/finalSoln.xml"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "Test_3D"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: dataFormat = "PLAIN_TEXT"


  LOGICAL :: EXPORT_FIELD

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: Err
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,CELLPATH,RyRModel
  INTEGER(CMISSIntg) :: NUMBER_OF_ATTRIBUTES,BOUNDARY_MARKER,ELE_ATTRIBUTES
  INTEGER :: st,i,NUMBER_OF_COORDS,NODES_PER_ELE
  INTEGER :: NUMBER_OF_NODES,element,NUMBER_OF_RYRS,NODE_BD_LABEL
  REAL(CMISSDP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums
  REAL(CMISSDP) :: nodex,nodey,nodez
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diganostic levels on for testing

  CALL MPI_BCAST(NUMBER_GLOBAL_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)
  

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  !CALL CMISSBasis_Initialise(Basis,Err)
  !CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a linear Lagrange basis
  !CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
  !Finish the creation of the basis
  !CALL CMISSBasis_CreateFinish(Basis,Err)

  !Start the creation of a basis (default is trilinear simplex)
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a linear simplex basis
  CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis,(/CMISS_Basis_Linear_Simplex_Interpolation, &
   &   CMISS_Basis_Linear_Simplex_Interpolation, CMISS_Basis_Linear_Simplex_Interpolation/),Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)

  !Start the creation of a mesh from a file in the region
  open(unit=10,file='test.node',status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeNums(NUMBER_OF_NODES,2))
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS))
    DO i = 1,NUMBER_OF_NODES
      READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeCoords(i,3),NodeNums(i,2)
      WRITE(*,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeCoords(i,3),NodeNums(i,2)
    ENDDO
  ENDIF
  CLOSE(10)
  !Read in elements
  OPEN(unit=11,file='test.ele',status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES
    ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,5))
    DO i = 1,NUMBER_OF_ELEMENTS
      READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5)
      WRITE(*,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5)
    ENDDO
  ENDIF 
  CLOSE(11)

  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  DO i = 1,NUMBER_OF_NODES
    CALL CMISSNodes_UserNumberSet(Nodes,i,NodeNums(i,1),Err)
    WRITE(*,*) i,NodeNums(i,1)
  ENDDO

  CALL CMISSNodes_CreateFinish(Nodes,Err)
  DO i = 1,NUMBER_OF_NODES
    CALL CMISSNodes_UserNumberGet(Nodes,i,NODE_NUMBER,Err)
    WRITE(*,*) i,NODE_NUMBER
  ENDDO
  
  WRITE(*,*) 'PAST NODES INITIALISATION'
  
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,1,Err)

  WRITE(*,*) 'PAST MESH  INITIALISATION'
  
  CALL CMISSMeshElements_Initialise(MeshElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  DO i = 1,NUMBER_OF_ELEMENTS
    element = ElemMap(i,1)
    !CALL CMISSMeshElements_UserNumberSet(MeshElements,i,ElemMap(i,1),Err)
    CALL CMISSMeshElements_NodesSet(MeshElements,i,(/ElemMap(i,2),ElemMap(i,3), &
     &   ElemMap(i,4),ElemMap(i,5)/),Err)

  ENDDO
  CALL CMISSMeshElements_CreateFinish(MeshElements,Err)
  CALL CMISSMesh_CreateFinish(Mesh,Err)  
  WRITE(*,*) 'PAST MESH ELEMENTS INITIALISATION'

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  WRITE(*,*) 'PAST DECOMPOSITION INITIALISATION'
  
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)
  WRITE(*,*) 'PAST GEOMETRIC FIELD INITIALISATION'  
  
  !Update the geometric field parameters
  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)
      nodey = NodeCoords(i,2)
      nodez = NodeCoords(i,3)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex,Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,3,nodez,Err)
     ENDIF
    ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)



  EXPORT_FIELD = .TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"Test_3DMesh","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Test_3DMesh","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF 

  !Calculate the volumes of each of the mesh elements within GEOMETRIC_PARAMETERS_CALCULATE inside field_routines.f90
  !store these volumes in an element field called volumes
  !Create a volume field
  CALL CMISSField_Initialise(ElemVolField,Err)
  CALL CMISSField_CreateStart(ElemVolFieldUserNumber,Region,ElemVolField,Err)
  CALL CMISSField_TypeSet(ElemVolField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(ElemVolField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(ElemVolField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(ElemVolField,1,Err)
  CALL CMISSField_VariableTypesSet(ElemVolField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMISSField_DataTypeSet(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
  CALL CMISSField_DimensionSet(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_NumberOfComponentsSet(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,"Element Volumes Field",Err)
  CALL CMISSField_ComponentMeshComponentGet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMISSField_ComponentMeshComponentSet(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be element constant
  CALL CMISSField_ComponentInterpolationSet(ELemVolField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,ERR)
  CALL CMISSField_CreateFinish(ElemVolField,Err)
  !Initialise Element Volumes to zero value
  CALL CMISSField_ComponentValuesInitialise(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  WRITE(*,*) 'PAST ELEM VOLS FIELD INITIALISATION'

  !assign element volumes calculated in the geometric parameters field to elemvolfield
  DO ne=1,NUMBER_OF_ELEMENTS
    CALL CMISSField_GeometricParametersElementVolumeGet(GeometricField,ElemMap(ne,1),elem_volume,Err)
    WRITE(*,*) "Elem number", ElemMap(ne,1) !print out element
    WRITE(*,*) "volume:", elem_volume !and it's volume
    CALL CMISSField_ParameterSetUpdateElement(ElemVolField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & ElemMap(ne,1),GeometricMeshComponent,elem_volume,Err)
  ENDDO
  

  !Find surrounding elements of a particular node.
  NODE_NUMBER=3_CMISSIntg
  CALL CMISSDecomposition_NodeNumberSurroundingElementsGet(Decomposition,NODE_NUMBER,1, &
    & number_surrounding_elements,Err)
  WRITE(*,*) "The Number of elements for node", NODE_NUMBER
  WRITE(*,*) "is", number_surrounding_elements
  CALL CMISSDecomposition_NodeSurroundingElementsGet(Decomposition,NODE_NUMBER,1, &
    & SURROUNDING_ELEMENTS,Err)
  DO ne=1,number_surrounding_elements
    WRITE(*,*) SURROUNDING_ELEMENTS(ne)
  ENDDO



  !FieldML Output
  !CALL CMISSFieldMLIO_Initialise(outputInfo,Err)
  !CALL CMISSFieldML_OutputCreate(Mesh,outputDirectory,basename,dataFormat,outputInfo,Err)
  !CALL CMISSFieldML_OutputAddImport(outputInfo,"coordinates.rc.3d",typeHandle,Err)
  !CALL CMISSFieldML_OutputAddField(outputInfo,baseName//".geometric",dataFormat,GeometricField, &
  !  & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSFieldML_OutputWrite(outputInfo,outputFilename,Err)
  !CALL CMISSFieldMLIO_Finalise(outputInfo,Err)


  

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP


END PROGRAM ELEMVOLSEXAMPLE
