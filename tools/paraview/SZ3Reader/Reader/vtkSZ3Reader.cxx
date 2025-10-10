#include "vtkSZ3Reader.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkAOSDataArrayTemplate.h"
#include "vtkImageAlgorithm.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPVInformationKeys.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <iostream>

vtkStandardNewMacro(vtkSZ3Reader);

vtkSZ3Reader::vtkSZ3Reader()
{
  this->FileName = nullptr;
  this->DomainDimensions[0] = 0;
  this->DomainDimensions[1] = 0;
  this->DomainDimensions[2] = 0;
  this->UseDoublePrecision = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

vtkSZ3Reader::~vtkSZ3Reader()
{
  this->SetFileName(nullptr);
}

void vtkSZ3Reader::SetDomainDimensions(int x, int y, int z)
{
  this->DomainDimensions[0] = x;
  this->DomainDimensions[1] = y;
  this->DomainDimensions[2] = z;
  this->Modified();
}

void vtkSZ3Reader::GetDomainDimensions(int& x, int& y, int& z)
{
  x = this->DomainDimensions[0];
  y = this->DomainDimensions[1];
  z = this->DomainDimensions[2];
}

int vtkSZ3Reader::RequestData(
  vtkInformation* /*request*/, vtkInformationVector** /*inputVector*/, vtkInformationVector* outputVector)
{

  std::cout << "FileName: " << (this->FileName ? this->FileName : "(none)") << std::endl;
  std::cout << "DomainDimensions: " << this->DomainDimensions[0] << ", " << this->DomainDimensions[1] << ", " << this->DomainDimensions[2] << std::endl;
  std::cout << "DoublePrecision: " << this->UseDoublePrecision << std::endl;

  if (!this->FileName) {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
  }

  if (this->DomainDimensions[0] <= 0 || this->DomainDimensions[1] <= 0 || this->DomainDimensions[2] <= 0) {
    vtkErrorMacro("Dimensions must be greater than zero.");
    return 0;
  }

  std::ifstream file(this->FileName, std::ios::binary);
  if(!file) {
    vtkErrorMacro("Could not open file: " << this->FileName);
    return 0;
  }

  file.seekg(0, std::ios::end);
  size_t compressedSize = file.tellg();
  file.seekg(0, std::ios::beg);
  std::vector<char> compressedBuffer(compressedSize);
  file.read(compressedBuffer.data(), compressedSize);
  file.close();

  auto conf = SZ3::Config(DomainDimensions[0], DomainDimensions[1], DomainDimensions[2]);

  vtkImageData* output = vtkImageData::GetData(outputVector);
  output->SetDimensions(this->DomainDimensions);
  
  if (this->UseDoublePrecision) {
    this->Decompress<double>(output, conf, compressedBuffer);
  }
  else {
    this->Decompress<float>(output, conf, compressedBuffer);
  }

  return 1;
}

template <typename T>
void vtkSZ3Reader::Decompress(
  vtkImageData* output, SZ3::Config& conf, std::vector<char>& compressedBuffer)
{
  using VtkArrayT = vtkAOSDataArrayTemplate<T>;
  vtkNew<VtkArrayT> dataArray;

  dataArray->SetNumberOfComponents(1);
  dataArray->SetNumberOfTuples(conf.num);
  dataArray->SetName("scalar");

  // Decompress directly into the VTK array's buffer
  T* decompressedPtr = static_cast<T*>(dataArray->GetVoidPointer(0));

  try {
    SZ_decompress<T>(
      conf, reinterpret_cast<const char*>(compressedBuffer.data()), compressedBuffer.size(), decompressedPtr);
  }
  catch (std::exception& e) {
    vtkErrorMacro("Decompression failed: " << e.what());
    return;
  }

  output->GetPointData()->AddArray(dataArray);
  output->GetPointData()->SetScalars(dataArray);
}

int vtkSZ3Reader::RequestInformation(
    vtkInformation* /*request*/,
    vtkInformationVector** /*inputVector*/,
    vtkInformationVector* outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // Set the whole extent
  int extent[6] = {0, 0, 0, 0, 0, 0};
  extent[1] = this->DomainDimensions[0] > 0 ? this->DomainDimensions[0] - 1 : 0;
  extent[3] = this->DomainDimensions[1] > 0 ? this->DomainDimensions[1] - 1 : 0;
  extent[5] = this->DomainDimensions[2] > 0 ? this->DomainDimensions[2] - 1 : 0;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

  // Set default origin and spacing
  double origin[3] = {0.0, 0.0, 0.0};
  double spacing[3] = {1.0, 1.0, 1.0};
  outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), spacing, 3);

  return 1;
}

void vtkSZ3Reader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " << (this->FileName ? this->FileName : "(none)") << "\n";
  os << indent << "DomainDimensions: " << this->DomainDimensions[0] << ", " << this->DomainDimensions[1] << ", " << this->DomainDimensions[2] << "\n";
  os << indent << "UseDoublePrecision: " << this->UseDoublePrecision << "\n";
}
