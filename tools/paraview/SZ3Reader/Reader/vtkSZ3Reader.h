#ifndef vtkSZ3Reader_h
#define vtkSZ3Reader_h

#include "vtkSZ3ReaderModule.h" // for export macro
#include "vtkDataObjectAlgorithm.h"
#include "vtkImageAlgorithm.h"

#include <SZ3/api/sz.hpp>
#include <string>

class VTKSZ3READER_EXPORT vtkSZ3Reader : public vtkImageAlgorithm
{
public:
  static vtkSZ3Reader* New();
  vtkTypeMacro(vtkSZ3Reader, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  void SetDomainDimensions(int x, int y, int z);
  void GetDomainDimensions(int& x, int& y, int& z);

  void SetDoublePrecision(int flag) { this->UseDoublePrecision = flag; this->Modified(); }
  int GetDoublePrecision() const { return this->UseDoublePrecision; }

protected:
  vtkSZ3Reader();
  ~vtkSZ3Reader();
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  char* FileName;
  int DomainDimensions[3];
  int UseDoublePrecision;

private:
  vtkSZ3Reader(const vtkSZ3Reader&);
  void operator=(const vtkSZ3Reader&);

  template <typename T>
  void Decompress(vtkImageData* output, SZ3::Config& conf, std::vector<char>& compressedBuffer);
};

#endif
