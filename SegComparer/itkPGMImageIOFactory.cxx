#include "itkPGMImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkPGMImageIO.h"
#include "itkVersion.h"

  
namespace itk
{
PGMImageIOFactory::PGMImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkPGMImageIO",
                         "PGM Image IO",
                         1,
                         CreateObjectFunction<PGMImageIO>::New());
}
  
PGMImageIOFactory::~PGMImageIOFactory()
{
}

const char* 
PGMImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}

const char* 
PGMImageIOFactory::GetDescription() const
{
  return "PGM ImageIO Factory, allows the loading of Netpbm PGM images into ITK";
}

} // end namespace itk
