!> \file
!> $Id: CantileverExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program to solve the cantilever problem with quadratic simplex elements and uses the experimental data from 
!> V Rajagopal et.al, 2007 IJNME for validating the mechanics code.
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
!> The Original Code is openCMISS
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

!> \example FiniteElasticity/SimplexElements/Cantilever/CantileverExample.f90
!! Example program to solve a mechanics equation for cantilever beam from V Rajagopal et.al IJNME, 2007 using openCMISS calls.
!! \htmlinclude FiniteElasticity/SimplexElements/Cantilever/history.html
!<

!> Main program
PROGRAM CANTILEVER_EXAMPLE


  USE OPENCMISS
  USE MPI
  USE FIELDML_API


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: SimplexBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: LinSimplexBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DisplacementFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MechPropertiesFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MechEquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: MechProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: GravityFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=21  
  INTEGER(CMISSIntg), PARAMETER :: MechEquationsSetFieldUserNumber=1337


  !CMISS variables
  TYPE(CMISSBasisType) :: SimplexBasis,LinSimplexBasis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: MechEquations
  TYPE(CMISSEquationsSetType) :: MechEquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DisplacementField,MechPropertiesField,GravityField,FibreField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSBoundaryConditionsType) :: MechBoundaryConditions
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: MechProblem
  TYPE(CMISSControlLoopType) :: LoadControlLoop
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  TYPE(CMISSSolverEquationsType) :: MechSolverEquations
  TYPE(CMISSFieldType) :: MechEquationsSetField
  TYPE(CMISSMeshElementsType) :: MeshElements
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh    
  !Program variables
  INTEGER(CMISSIntg) :: CONDITION,st,i,NODES_PER_ELE,NUMBER_OF_COORDS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS,NODE_NUMBER,NodeDomain,NUMBER_GLOBAL_ELEMENTS  
  INTEGER(CMISSIntg) :: NUMBER_OF_ATTRIBUTES,BOUNDARY_MARKER,ELE_ATTRIBUTES
  INTEGER :: NUMBER_OF_NODES,NUMBER_OF_ELEMENTS,element,NODE_BD_LABEL
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi
  INTEGER(CMISSIntg) :: MEMBRANE_MARKER
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums
  REAL(CMISSDP) :: nodex,nodey,nodez
  REAL(CMISSDP), PARAMETER :: LENGTH=50.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=25.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=25.0_CMISSDP
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER :: node,STARTNODE,ENDNODE
  INTEGER(CMISSIntg) :: constantModelIndex,component_idx,user_node,node_idx,NodeNumber
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ScalingType=CMISS_FIELD_NO_SCALING
  INTEGER(CMISSIntg), PARAMETER :: NumberOfLoadIncrements=3
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: MechEquationsSetIndex
  INTEGER(CMISSIntg) :: Err

  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP), PARAMETER :: Density=9.0E-4_CMISSDP !in g mm^-3
  REAL(CMISSDP), PARAMETER :: Gravity(3)=[0.0_CMISSDP,0.0_CMISSDP,-9.81_CMISSDP] !in m s^-2
  REAL(CMISSDP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  CHARACTER(250) :: MESHFILE,NODEFILE,ELEMFILE

  !FieldML output parameters
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = "."
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = outputDirectory//"/cellv1_6x4x4.xml"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "cellv1"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: dataFormat = "PLAIN_TEXT"
  !FieldML parsing variables
  LOGICAL :: EXPORT_FIELD=.TRUE.
  TYPE(CMISSFieldMLIOType) :: fieldmlInfo, outputInfo
  INTEGER(CMISSIntg) :: meshComponentCount  
  INTEGER(CMISSIntg) :: typeHandle
  INTEGER(CMISSIntg) :: coordinateCount
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_GLOBAL_X_ELEMENTS=4
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_GLOBAL_Y_ELEMENTS=4
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_GLOBAL_Z_ELEMENTS=4 

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  

  
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


!_________________________________________________________________________________________________

  EXPORT_FIELD=.TRUE.

! THIS IS AN EXAMPLE VALIDATING SIMPLEX BASIS FUNCTIONS FOR MECHANICS USING THE CANTILEVER EXPERIMENTAL DATA FROM V RAJAGOPAL, IJNME 2007
! THE MODEL DEFLECTION OF 35 mm AGREES WELL WITH THE EXPERIMENTAL DEFLECTION MEASUREMENT OF 34.6 mm. FOR TESTING PURPOSES THE SIMULATION
! USES A LOW MESH RESOLUTION. INCREASE THE NUMBER OF ELEMENTS IN THE MESH IF YOU WANT TO SEE THE MATCH WITH THE EXPERIMENTAL MEASUREMENT.
!_________________________________________________________________________________________________
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

!_________________________________________________________________________________________________  


  !Define basis functions - tri-quadratic simplex and linear hydrostatic pressure.

  CALL CMISSBasis_Initialise(SimplexBasis,Err)
  CALL CMISSBasis_CreateStart(SimplexBasisUserNumber,SimplexBasis,Err)
  CALL CMISSBasis_TypeSet(SimplexBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(SimplexBasis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(SimplexBasis,(/CMISS_Basis_Quadratic_Simplex_Interpolation, &
   &   CMISS_Basis_Quadratic_Simplex_Interpolation, CMISS_Basis_Quadratic_Simplex_Interpolation/),Err)
  CALL CMISSBasis_QuadratureOrderSet(SimplexBasis,3,Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(SimplexBasis,.TRUE.,Err)
  CALL CMISSBasis_CreateFinish(SimplexBasis,Err)

  !Define basis function - Simplex tri-linear (pressure)
  CALL CMISSBasis_Initialise(LinSimplexBasis,Err)
  CALL CMISSBasis_CreateStart(LinSimplexBasisUserNumber,LinSimplexBasis,Err) 
  CALL CMISSBasis_TypeSet(LinSimplexBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinSimplexBasis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(LinSimplexBasis,(/CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION/),Err)
  CALL CMISSBasis_CreateFinish(LinSimplexBasis,Err)

!_________________________________________________________________________________________________  
  !Time to create a mesh - wohoo!


  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[SimplexBasis,LinSimplexBasis],Err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/LENGTH,WIDTH,HEIGHT/),Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS, &
    NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !!Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !!Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !!Set the domain to be used by the field components. We have 3 field components in 1 mesh component
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !!Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  CALL CMISSMeshElements_Initialise(MeshElements,Err)
  CALL CMISSMesh_ElementsGet(Mesh,1,MeshElements,Err)
  CALL CMISSMesh_NumberOfElementsGet(Mesh,NUMBER_OF_ELEMENTS,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"SimplexCantilever","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"SimplexCantilever","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

!______________________________________________________________________________________________________________
  !Some additional fields
  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

!______________________________________________________________________________________________________________
  !Set up mechanics equations

  !Create the equations_set
  CALL CMISSField_Initialise(MechEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(MechEquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
    & MechEquationsSetFieldUserNumber, &
    & MechEquationsSetField, &
    & MechEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(MechEquationsSet,Err)


  !Create the displacement field
  CALL CMISSField_Initialise(DisplacementField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(MechEquationsSet,DisplacementFieldUserNumber,DisplacementField,Err)
  CALL CMISSField_VariableLabelSet(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,"Displacements",Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentMeshComponentSet(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(DisplacementField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO

  !Add hydrostatic pressure as nodally varying field component and match it to mesh component 2 for linear pressure interp.
  CALL CMISSField_ComponentMeshComponentSet(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,4,2,Err)
  CALL CMISSField_ComponentMeshComponentSet(DisplacementField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  !setting the 4th component of the displacement field to nodally based.
  CALL CMISSField_ComponentInterpolationSet(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,4, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION, &
    & Err)
  CALL CMISSField_ComponentInterpolationSet(DisplacementField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(MechEquationsSet,Err)

  !Create the mechanical properties field
  CALL CMISSField_Initialise(MechPropertiesField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(MechEquationsSet,MechPropertiesFieldUserNumber,MechPropertiesField,Err)
  CALL CMISSField_VariableLabelSet(MechPropertiesField,CMISS_FIELD_U_VARIABLE_TYPE,"Stiffness",Err)
  CALL CMISSField_VariableLabelSet(MechPropertiesField,CMISS_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(MechEquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL CMISSField_ComponentValuesInitialise(MechPropertiesField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,0.64_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MechPropertiesField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,2,0.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MechPropertiesField,CMISS_FIELD_V_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,Density,Err)

  !Create the source field with the gravity vector
  CALL CMISSField_Initialise(GravityField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(MechEquationsSet,GravityFieldUserNumber,GravityField,Err)
  CALL CMISSEquationsSet_SourceCreateFinish(MechEquationsSet,Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentValuesInitialise(GravityField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

  !Create the equations set equations
  CALL CMISSEquations_Initialise(MechEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(MechEquationsSet,MechEquations,Err)
  CALL CMISSEquations_SparsityTypeSet(MechEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(MechEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(MechEquationsSet,Err)
!______________________________________________________________________________________________________________
  
  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & -0.64_CMISSDP, &
    & Err)
  CALL CMISSField_ParameterSetUpdateStart(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!______________________________________________________________________________________________________________

  !Set up the problem
  CALL CMISSProblem_Initialise(MechProblem,Err)
  CALL CMISSProblem_CreateStart(MechProblemUserNumber,MechProblem,Err)
  CALL CMISSProblem_SpecificationSet(MechProblem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(MechProblem,Err)

!______________________________________________________________________________________________________________
  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(MechProblem,Err)
  CALL CMISSControlLoop_Initialise(LoadControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(MechProblem,CMISS_CONTROL_LOOP_NODE,LoadControlLoop,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(LoadControlLoop,NumberOfLoadIncrements,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(MechProblem,Err)

!______________________________________________________________________________________________________________
  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(MechProblem,Err)
  CALL CMISSProblem_SolverGet(MechProblem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_NewtonTypeSet(Solver,CMISS_SOLVER_NEWTON_LINESEARCH,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(Solver,1e-3_CMISSDP,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(MechProblem,Err)
!______________________________________________________________________________________________________________
  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(MechSolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(MechProblem,Err)
  CALL CMISSProblem_SolverGet(MechProblem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,MechSolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(MechSolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(MechSolverEquations,MechEquationsSet,MechEquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(MechProblem,Err)
!______________________________________________________________________________________________________________

  !Set up loading and boundary conditions
  !Prescribe boundary conditions (relative nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(MechBoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(MechSolverEquations,MechBoundaryConditions,Err)


  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes in x, y and z
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,3
        CALL CMISSBoundaryConditions_AddNode(MechBoundaryConditions,DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & component_idx,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDDO
    ENDIF
  ENDDO

   CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(MechSolverEquations,Err)

!______________________________________________________________________________________________________________
  !Solve problem

  WRITE(*,*) "About to solve problem"
  CALL CMISSProblem_Solve(MechProblem,Err)
!______________________________________________________________________________________________________________
  !Output solutions and geometries

  IF(EXPORT_FIELD) THEN

    CALL CMISSField_ParametersToFieldParametersComponentCopy(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE, & 
      & CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
    CALL CMISSField_ParametersToFieldParametersComponentCopy(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE, & 
      & CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
    CALL CMISSField_ParametersToFieldParametersComponentCopy(DisplacementField,CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)

    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"Cantilever","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Cantilever","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)

  ENDIF
  

  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM CANTILEVER_EXAMPLE
