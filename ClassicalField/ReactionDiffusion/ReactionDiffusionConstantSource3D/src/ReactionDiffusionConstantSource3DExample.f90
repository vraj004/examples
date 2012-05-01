!> \file
!> $Id: DiffusionExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program to solve a diffusion equation using openCMISS calls.
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

!> \example ClassicalField/ReactionDiffusion/ReactionDiffusionConstantSource3D/src/ReactionDiffusionConstantSource3DExample.f90
!! Example program to solve a diffusion equation using openCMISS calls.
!! \htmlinclude ClassicalField/ReactionDiffusion/ReactionDiffusionConstantSource3D/history.html
!<

!> Main program
PROGRAM REACTIONDIFFUSIONCONSTANTSOURCE3DEXAMPLE


  USE OPENCMISS
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: LENGTH=100.0_CMISSDP,WIDTH=20.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=20.0_CMISSDP
  
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


  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,CONDITION,NUMBER_GLOBAL_Y_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS,NODE_NUMBER,NodeDomain,NUMBER_GLOBAL_ELEMENTS
  INTEGER(CMISSIntg),DIMENSION(50) :: BCNODES
  INTEGER(CMISSIntg),DIMENSION(9) :: RYRNODES
  INTEGER(CMISSIntg) :: MPI_IERROR,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER :: node,STARTNODE,ENDNODE
  REAL(CMISSDP) :: VALUE
  INTEGER(CMISSIntg) :: constantModelIndex

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

  NUMBER_GLOBAL_X_ELEMENTS=6
  NUMBER_GLOBAL_Y_ELEMENTS = 4
  NUMBER_GLOBAL_Z_ELEMENTS = 4
  NUMBER_GLOBAL_ELEMENTS = NUMBER_GLOBAL_X_ELEMENTS*NUMBER_GLOBAL_Y_ELEMENTS*NUMBER_GLOBAL_Z_ELEMENTS
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
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a linear Lagrange basis
  CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular 1D mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/LENGTH,WIDTH,HEIGHT/),Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh, &
    & (/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  
  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
  
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
  
  
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  
  !Create the cellml reaction with split reaction diffusion equations_set 
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region, &
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  CALL CMISSField_VariableLabelSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 1,0.39_CMISSDP,Err) !diff coeff in x
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 2,0.39_CMISSDP,Err) !diff coeff in y
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 3,0.39_CMISSDP,Err) !diff coeff in z
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 3,1.0_CMISSDP,Err) ! storage coefficient

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  CALL CMISSField_VariableLabelSet(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,"Source Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(SourceField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)

  !set initial value of the dependent field/state variable, ca.
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)
 
  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Create the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a strang split reaction diffusion problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMISS_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,0.0_CMISSDP,0.1_CMISSDP,0.01_CMISSDP,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)


  !Set up the problem solver for pde solving
  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  !Solver is the dynamic solver for solving the parabolic equation
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL CMISSSolver_DynamicThetaSet(Solver,1.0_CMISSDP,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL CMISSSolver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LibraryTypeSet(LinearSolver,CMISS_SOLVER_LAPACK_LIBRARY,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSSolver_LinearDirectTypeSet(LinearSolver,CMISS_SOLVER_DIRECT_LU,Err)
  !CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverCMISSLibrary,Err)
  !CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverMUMPSLibrary,Err)
  !CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,10000,Err)

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Create the equations set boundary conditions

  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  BCNODES = (/1,8,15,22,29,36,43,50,57,64,71,78,85,92,99,106,113,120,127,134,141,148,155,162,169, &
  & 7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175/)
  DO node=1,50
    NODE_NUMBER = BCNODES(node)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CONDITION = CMISS_BOUNDARY_CONDITION_FIXED
      VALUE=5.0_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField, &
       & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
       & NODE_NUMBER,1,CONDITION,VALUE,Err)
    ENDIF
  ENDDO


  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"ReactionDiffusion3D","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"ReactionDiffusion3D","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)

  ENDIF
  

  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM REACTIONDIFFUSIONCONSTANTSOURCE3DEXAMPLE
