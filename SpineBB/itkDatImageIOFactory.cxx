#include "itkDatImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkDatImageIO.h"
#include "itkVersion.h"

  
namespace itk
{
DatImageIOFactory::DatImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkDatImageIO",
                         "Dat Image IO",
                         1,
                         CreateObjectFunction<DatImageIO>::New());
}
  
DatImageIOFactory::~DatImageIOFactory()
{
}

const char* 
DatImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}

const char* 
DatImageIOFactory::GetDescription() const
{
  return "Dat ImageIO Factory, allows the loading of Dat images into insight";
}

} // end namespace itk
