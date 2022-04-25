
#include <iostream>
#include "vtkPolyDataToitkMesh.h"

#ifndef vtkDoubleType
#define vtkDoubleType double
#endif

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

vtkPolyDataToitkMesh
::vtkPolyDataToitkMesh()
{
  m_itkMesh = TriangleMeshType::New();
  m_PolyData = vtkSmartPointer<vtkPolyData>::New();

}

vtkPolyDataToitkMesh
::~vtkPolyDataToitkMesh()
{

}

void
vtkPolyDataToitkMesh
::SetInput(vtkSmartPointer<vtkPolyData> polydata)
{
  m_PolyData = polydata;
  this->ConvertvtkToitk();
}

vtkPolyDataToitkMesh::TriangleMeshType *
vtkPolyDataToitkMesh
::GetOutput()
{
  return m_itkMesh;
}

void
vtkPolyDataToitkMesh
::ConvertvtkToitk()
{
  //
  // Transfer the points from the vtkPolyData into the itk::Mesh
  //
  const unsigned int numberOfPoints = m_PolyData->GetNumberOfPoints();
  vtkSmartPointer<vtkPoints>        vtkpoints =  m_PolyData->GetPoints();

  m_itkMesh->GetPoints()->Reserve( numberOfPoints );
  for( unsigned int p = 0; p < numberOfPoints; p++ )
    {

    #if VTK_MAJOR_VERSION > 5
    double * apoint = vtkpoints->GetPoint( p );
    #else
    vtkFloatingPointType* apoint = vtkpoints->GetPoint( p );
    #endif
    m_itkMesh->SetPoint( p, TriangleMeshType::PointType( apoint ) );

    // Need to convert the point to PoinType
    TriangleMeshType::PointType pt;
    for( unsigned int i = 0; i < 3; i++ )
      {
      pt[i] = apoint[i];
      }
    m_itkMesh->SetPoint( p, pt);

    }
  //
  // Transfer the cells from the vtkPolyData into the itk::Mesh
  //
  vtkSmartPointer<vtkCellArray> triangleStrips = m_PolyData->GetStrips();


  vtkIdType numberOfCellPoints;

  //
  // First count the total number of triangles from all the triangle strips.
  //
  unsigned int numberOfTriangles = 0;

  /*triangleStrips->InitTraversal();
  while( triangleStrips->GetNextCell( numberOfCellPoints, cellPoints ) )
    {
    numberOfTriangles += numberOfCellPoints - 2;
  }*/

  vtkSmartPointer<vtkCellArray> polygons = m_PolyData->GetPolys();

  /*polygons->InitTraversal();

  while( polygons->GetNextCell( numberOfCellPoints, cellPoints ) )
    {
    if( numberOfCellPoints == 3 )
      {
      numberOfTriangles++;
      }
    }*/

    vtkIdType cellSizes = 0;
    for (vtkIdType i=0;i<triangleStrips->GetNumberOfCells(); i++)
    {
      numberOfTriangles+=triangleStrips->GetCellSize(i);
      if (triangleStrips->GetCellSize(i) > cellSizes)
        cellSizes = triangleStrips->GetCellSize(i);
    }

    for (vtkIdType i=0;i<polygons->GetNumberOfCells(); i++)
    {
      if (polygons->GetCellSize(i) == 3)
        numberOfTriangles++;

      if (polygons->GetCellSize(i) > cellSizes)
        cellSizes = polygons->GetCellSize(i);
    }


  //
  // Reserve memory in the itk::Mesh for all those triangles
  //
  m_itkMesh->GetCells()->Reserve( numberOfTriangles );

  //
  // Copy the triangles from vtkPolyData into the itk::Mesh
  //
  //

  typedef TriangleMeshType::CellType CellType;

  typedef itk::TriangleCell<CellType> TriangleCellType;

  vtkNew<vtkIdList> cellPoints;

  // first copy the triangle strips
  int cellId = 0;
  triangleStrips->InitTraversal();
  while( numberOfCellPoints = triangleStrips->GetNextCell( cellPoints ) )
    {
    unsigned int numberOfTrianglesInStrip = numberOfCellPoints - 2;

    uint64_t pointIds[3];
    pointIds[0] = cellPoints->GetId(0);
    pointIds[1] = cellPoints->GetId(1);
    pointIds[2] = cellPoints->GetId(2);
    for( unsigned int t = 0; t < numberOfTrianglesInStrip; t++ )
      {
      TriangleMeshType::CellAutoPointer c;
      TriangleCellType * tcell = new TriangleCellType;
      TriangleCellType::PointIdentifier itkPts[3];
      for (int ii = 0; ii < 3; ++ii)
        {
        itkPts[ii] = static_cast<TriangleCellType::PointIdentifier>(pointIds[ii]);
        }
      tcell->SetPointIds( itkPts );
      c.TakeOwnership( tcell );
      m_itkMesh->SetCell( cellId, c );
      cellId++;
      pointIds[0] = pointIds[1];
      pointIds[1] = pointIds[2];
      pointIds[2] = cellPoints->GetId(t+3);
      }

    }

  // then copy the triangles
  polygons->InitTraversal();
  while( numberOfCellPoints = polygons->GetNextCell( cellPoints ) )
    {
    if( numberOfCellPoints != 3 ) // skip any non-triangle.
      {
      continue;
      }
    TriangleMeshType::CellAutoPointer c;
    TriangleCellType * t = new TriangleCellType;
    TriangleCellType::PointIdentifier itkPts[3];
    for (int ii = 0; ii < numberOfCellPoints; ++ii)
      {
      itkPts[ii] = static_cast<TriangleCellType::PointIdentifier>(cellPoints->GetId(ii));
      }
    t->SetPointIds( itkPts );
    c.TakeOwnership( t );
    m_itkMesh->SetCell( cellId, c );
    cellId++;
    }
}
