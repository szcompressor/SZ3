set(classes
  vtkSZ3Reader)

vtk_module_add_module(SZ3Reader::vtkSZ3Reader
  CLASSES ${classes})
vtk_module_link(SZ3Reader::vtkSZ3Reader
  PUBLIC
    SZ3)

paraview_add_server_manager_xmls(
  XMLS 
    SZ3Reader.xml)
